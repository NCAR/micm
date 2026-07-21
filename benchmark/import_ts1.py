#!/usr/bin/env python
"""Import the MOZART-TS1 mechanism from a CAM chem_mech.in into a simple
line-based table read by the ts1_dae benchmark.

Source: CAM (ESCOMP/CAM, branch cam_cesm2_2_rel)
  src/chemistry/pp_trop_strat_mam4_vbsext/chem_mech.in
  — tag "TS1-fullVBS for CESM2.2" (MZ264_TS1_fullVBS_20190611), the
  MOZART-TS1 chemistry of Emmons et al. (2020).

Import policy (documented approximations, all counted in the header):
  * Arrhenius rates (1 or 2 parameters): exact — k = A * exp(C/T).
  * Troe rates (5 parameters, MOZART convention k0 = a*(300/T)^n):
    exact — micm k0_B = -n, kinf_B = -n_inf; M is implicit in the
    falloff kernel, so M is dropped from those reactions' species lists.
  * Photolysis: emitted as user-defined rates (label = the j-tag);
    the benchmark assigns fixed clear-sky midday values to the major
    channels and a small default to the rest (steady-J benchmark).
  * usr_* / het_* rates (no rate parameters in the file): EXCLUDED —
    these encode CAM-specific functional forms and heterogeneous
    chemistry that are out of scope for a gas-phase solver benchmark.
  * Species left with no remaining reactions (mostly aerosol tracers
    coupled only through excluded reactions) are dropped.

Output format (benchmark/data/ts1_mechanism.txt), one record per line:
  S <species>
  R ARRH <label> <A> <C> ; <reactants> ; <coef*product ...>
  R TROE <label> <k0_A> <k0_B> <kinf_A> <kinf_B> <Fc> ; ... ; ...
  R PHOT <label> ; <reactant> ; <coef*product ...>
"""

from __future__ import annotations

import re
import sys
from pathlib import Path


def parse_species(lines: list[str]) -> list[str]:
    species = []
    for entry in ",".join(lines).split(","):
        entry = entry.strip()
        if not entry:
            continue
        name = entry.split("->")[0].strip()
        if name:
            species.append(name)
    return species


def parse_side(text: str, drop: set[str]) -> list[tuple[float, str]]:
    terms = []
    for token in text.split("+"):
        token = token.strip()
        if not token:
            continue
        if "*" in token:
            coef_text, name = token.split("*", 1)
            coef = float(coef_text)
        else:
            coef, name = 1.0, token
        name = name.strip()
        if name in drop:
            continue
        terms.append((coef, name))
    return terms


def main() -> None:
    src = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).parent / "data" / "ts1_chem_mech.in"
    out = Path(sys.argv[2]) if len(sys.argv) > 2 else Path(__file__).parent / "data" / "ts1_mechanism.txt"
    text = src.read_text().splitlines()

    def section(start_tag: str, end_tag: str) -> list[str]:
        start = next(i for i, l in enumerate(text) if l.strip() == start_tag)
        end = next(i for i, l in enumerate(text) if l.strip() == end_tag)
        return text[start + 1 : end]

    species = parse_species([l for l in section("Solution", "End Solution") if not l.strip().startswith("*")])
    fixed = parse_species([l for l in section("Fixed", "End Fixed") if not l.strip().startswith("*")])
    solution = set(species) | set(fixed)

    records: list[str] = []
    excluded: list[str] = []
    used_species: set[str] = set()

    tag_re = re.compile(r"^\[([^\]]+)\]\s*(.*)$")

    def emit(kind: str, label: str, params: list[float], lhs: str, rhs: str, drop_m: bool) -> None:
        drop = {"hv"} | ({"M"} if drop_m else set())
        reactants = parse_side(lhs, drop)
        products = parse_side(rhs, drop)
        names = [n for _, n in reactants] + [n for _, n in products]
        unknown = [n for n in names if n not in solution and n != "M"]
        if unknown:
            excluded.append(f"{label} (unknown species {unknown})")
            return
        used_species.update(n for n in names)
        used_species.update({"M"} if not drop_m and "M" in lhs else set())
        reactant_text = " ".join(n for c, n in reactants for _ in range(int(round(c))))
        product_text = " ".join(f"{c:.6g}*{n}" for c, n in products)
        param_text = (" " + " ".join(f"{p:.10g}" for p in params)) if params else ""
        records.append(f"R {kind} {label}{param_text} ; {reactant_text} ; {product_text}")

    # Photolysis: [tag] SPECIES + hv -> products
    for line in section("Photolysis", "End Photolysis"):
        line = line.strip()
        if not line or line.startswith("*"):
            continue
        match = tag_re.match(line)
        if not match:
            continue
        label = match.group(1).split(",")[0].split("=")[0].strip()
        lhs, rhs = match.group(2).split("->", 1)
        emit("PHOT", label, [], lhs, rhs, drop_m=False)

    # Kinetic reactions.
    for line in section("Reactions", "End Reactions"):
        line = line.strip()
        if not line or line.startswith("*"):
            continue
        match = tag_re.match(line)
        if not match:
            continue
        label = match.group(1).split(",")[0].strip()
        body = match.group(2)
        if ";" in body:
            equation, rate_text = body.split(";", 1)
            params = [float(p) for p in rate_text.replace(",", " ").split()]
        else:
            equation, params = body, []
        # Three usr_ rates are chemically load-bearing (ozone formation,
        # CO + OH, HO2 self-reaction) and have standard JPL forms; they are
        # imported with those forms (HO2 + HO2 simplified to its low-pressure
        # channel, dropping the M and H2O enhancements). All other usr_/het_
        # rates stay excluded.
        special = {
            # O + O2 + M -> O3 + M: k = 6.0e-34 * (T/300)^-2.4 * [M] (M as reactant)
            "usr_O_O2": ("ARRH", [6.0e-34, 0.0, -2.4], False),
            # CO + OH: k = 1.5e-13 * (1 + 0.6 * P_atm); micm E in Pa^-1
            "usr_CO_OH_b": ("ARRH", [1.5e-13, 0.0, 0.0, 0.6 / 101325.0], False),
            # HO2 + HO2 -> H2O2 + O2 (low-pressure channel only)
            "usr_HO2_HO2": ("ARRH", [3.0e-13, 460.0], False),
        }
        if label in special:
            kind, sparams, drop_m = special[label]
            lhs, rhs = equation.split("->", 1)
            emit(kind, label, sparams, lhs, rhs, drop_m=drop_m)
            continue
        if label.startswith("usr_") or label.startswith("het") or not params:
            excluded.append(label)
            continue
        lhs, rhs = equation.split("->", 1)
        if len(params) == 1:
            emit("ARRH", label, [params[0], 0.0], lhs, rhs, drop_m=False)
        elif len(params) == 2:
            emit("ARRH", label, [params[0], params[1]], lhs, rhs, drop_m=False)
        elif len(params) == 5:
            k0_a, k0_n, kinf_a, kinf_n, fc = params
            # MOZART: k0 = a*(300/T)^n  ->  micm (T/300)^B with B = -n.
            emit("TROE", label, [k0_a, -k0_n, kinf_a, -kinf_n, fc], lhs, rhs, drop_m=True)
        else:
            excluded.append(f"{label} (unhandled rate with {len(params)} parameters)")

    kept_species = [s for s in species if s in used_species]
    with open(out, "w") as f:
        f.write(f"# MOZART-TS1 (TS1-fullVBS for CESM2.2, MZ264_TS1_fullVBS_20190611)\n")
        f.write(f"# imported by import_ts1.py; species kept {len(kept_species)}/{len(species)} plus fixed {fixed}; ")
        f.write(f"reactions kept {len(records)}; excluded {len(excluded)} (usr_/het_/unrated):\n")
        f.write("#   " + " ".join(sorted(set(e.split(' ')[0] for e in excluded))) + "\n")
        for s in fixed:
            f.write(f"F {s}\n")
        for s in kept_species:
            f.write(f"S {s}\n")
        for r in records:
            f.write(r + "\n")
    print(f"wrote {out}: {len(kept_species)} species + {len(fixed)} fixed, {len(records)} reactions, {len(excluded)} excluded")


if __name__ == "__main__":
    main()
