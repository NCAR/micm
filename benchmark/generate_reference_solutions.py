#!/usr/bin/env python
"""External reference solutions for the MICM DAE/ODE benchmarks.

Integrates the Robertson and tropospheric O3-NOx-HOx benchmark systems with an
independent solver (scipy's Radau implicit Runge-Kutta) at tight tolerance and
writes the trajectories at the benchmarks' output times to CSV. The benchmark
programs measure their errors against these files instead of against a
tight-tolerance run of the same solver, removing the possibility that a shared
systematic error cancels in the comparison (phase 2 of
docs/superpowers/plans/2026-07-20-dae-performance-and-rigor.md).

The chemistry here must match benchmark/robertson_system.hpp and
benchmark/tropospheric_system.hpp exactly: MICM Arrhenius rate constants
k = A * exp(C/T) * (T/300)^B at the fixed box temperature, the same fixed
photolysis rates, the same stoichiometry (including reservoir species), and
the same QSSA-manifold initial radicals (damped fixed point, 500 iterations,
0.5 damping).

Usage:
    python generate_reference_solutions.py [output_dir]   # default: ./data
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp


def output_times(t_first: float, t_last: float, n: int) -> np.ndarray:
    """Replicates OutputTimes() in the benchmark drivers."""
    lo, hi = np.log10(t_first), np.log10(t_last)
    return 10.0 ** (lo + (hi - lo) * np.arange(n) / (n - 1))


def write_csv(path: Path, times: np.ndarray, names: list[str], solution: np.ndarray) -> None:
    with open(path, "w") as f:
        f.write("time," + ",".join(names) + "\n")
        for i, t in enumerate(times):
            f.write(f"{t:.17e}," + ",".join(f"{v:.17e}" for v in solution[:, i]) + "\n")
    print(f"wrote {path}")


# ---------------------------------------------------------------------------
# Robertson (robertson_system.hpp): A -> B, B+B -> B+C, B+C -> A+C
# ---------------------------------------------------------------------------

def robertson_reference(out_dir: Path) -> None:
    k1, k2, k3 = 0.04, 3.0e7, 1.0e4  # benchmark defaults; k2 is the default stiffness

    def rhs(_t, y):
        a, b, c = y
        return [-k1 * a + k3 * b * c, k1 * a - k2 * b * b - k3 * b * c, k2 * b * b]

    def jac(_t, y):
        _a, b, c = y
        return [
            [-k1, k3 * c, k3 * b],
            [k1, -2.0 * k2 * b - k3 * c, -k3 * b],
            [0.0, 2.0 * k2 * b, 0.0],
        ]

    times = output_times(1.0e-3, 1.0e6, 19)
    result = solve_ivp(
        rhs,
        (0.0, times[-1]),
        [1.0, 0.0, 0.0],
        method="Radau",
        jac=jac,
        rtol=1.0e-12,
        atol=1.0e-14,
        t_eval=times,
    )
    if not result.success:
        raise RuntimeError(f"Robertson reference failed: {result.message}")
    write_csv(out_dir / "robertson_reference.csv", times, ["A", "B", "C"], result.y)

    # Long-time reference for the conservation-DAE benchmark (IVP test-set
    # style, t -> 1e11): B falls to ~8e-14, so its absolute tolerance must sit
    # well below that for the reference to resolve it relatively.
    long_times = output_times(1.0e-3, 1.0e11, 29)
    long_result = solve_ivp(
        rhs,
        (0.0, long_times[-1]),
        [1.0, 0.0, 0.0],
        method="Radau",
        jac=jac,
        rtol=1.0e-13,
        atol=[1.0e-19, 1.0e-22, 1.0e-14],
        t_eval=long_times,
    )
    if not long_result.success:
        raise RuntimeError(f"Robertson long-time reference failed: {long_result.message}")
    write_csv(out_dir / "robertson_longtime_reference.csv", long_times, ["A", "B", "C"], long_result.y)


# ---------------------------------------------------------------------------
# Tropospheric O3-NOx-HOx (tropospheric_system.hpp)
# ---------------------------------------------------------------------------

T_BOX = 298.15
N_M = 2.5e19
N_N2 = 0.78 * N_M
N_O2 = 0.21 * N_M
N_H2O = 0.01 * N_M
N_O3_0 = 7.5e11
N_NO_0 = 2.5e10
N_NO2_0 = 2.5e10
N_CO_0 = 2.5e12
N_CH4_0 = 4.4e13
N_CO2_0 = 1.0e16

JNO2 = 8.0e-3
JO1D = 3.0e-5
JO3P = 4.0e-4

SPECIES = ["O1D", "O", "OH", "HO2", "O3", "NO", "NO2", "CO", "CH4", "H2O2", "HNO3", "CO2", "O2", "N2", "M", "H2O"]


def arrhenius(A: float, C: float = 0.0, B: float = 0.0, T: float = T_BOX) -> float:
    """MICM CalculateArrhenius with default D = 300 K and E = 0."""
    return A * np.exp(C / T) * (T / 300.0) ** B


K = {
    "O1D_N2": arrhenius(2.15e-11, C=110),
    "O1D_O2": arrhenius(3.30e-11, C=55),
    "O1D_H2O": arrhenius(1.63e-10, C=60),
    "O_O2_M": arrhenius(6.0e-34, B=-2.4),
    "O_O3": arrhenius(8.0e-12, C=-2060),
    "O3_NO": arrhenius(3.0e-12, C=-1500),
    "OH_CO": arrhenius(1.5e-13),
    "OH_CH4": arrhenius(2.45e-12, C=-1775),
    "HO2_NO": arrhenius(3.3e-12, C=270),
    "OH_NO2": arrhenius(1.0e-11),
    "HO2_HO2": arrhenius(3.0e-13, C=460),
}


def project_radicals(O3, NO, NO2, CO, CH4, N2, O2, M, H2O):
    """Replicates tropospheric::ProjectRadicals (damped fixed point)."""
    O1D = (JO1D * O3) / (K["O1D_N2"] * N2 + K["O1D_O2"] * O2 + K["O1D_H2O"] * H2O)
    O = (JNO2 * NO2 + JO3P * O3 + (K["O1D_N2"] * N2 + K["O1D_O2"] * O2) * O1D) / (K["O_O2_M"] * O2 * M + K["O_O3"] * O3)

    oh_loss = K["OH_CO"] * CO + K["OH_CH4"] * CH4 + K["OH_NO2"] * NO2
    oh_to_ho2 = K["OH_CO"] * CO + K["OH_CH4"] * CH4
    oh_src = 2.0 * K["O1D_H2O"] * H2O * O1D
    OH = oh_src / oh_loss
    HO2 = (oh_to_ho2 * OH) / (K["HO2_NO"] * NO + 1.0e-30)
    for _ in range(500):
        a = 2.0 * K["HO2_HO2"]
        b = K["HO2_NO"] * NO
        cc = oh_to_ho2 * OH
        HO2_new = (2.0 * cc) / (b + np.sqrt(b * b + 4.0 * a * cc))
        OH_new = (oh_src + K["HO2_NO"] * NO * HO2_new) / oh_loss
        OH = 0.5 * OH + 0.5 * OH_new
        HO2 = 0.5 * HO2 + 0.5 * HO2_new
    return O1D, O, OH, HO2


def tropospheric_rhs(_t, y):
    O1D, O, OH, HO2, O3, NO, NO2, CO, CH4, H2O2, HNO3, _CO2, O2, N2, M, H2O = y

    r_p1 = JNO2 * NO2
    r_p2 = JO1D * O3
    r_p3 = JO3P * O3
    r1 = K["O1D_N2"] * O1D * N2
    r2 = K["O1D_O2"] * O1D * O2
    r3 = K["O1D_H2O"] * O1D * H2O
    r4 = K["O_O2_M"] * O * O2 * M
    r5 = K["O_O3"] * O * O3
    r6 = K["O3_NO"] * O3 * NO
    r7 = K["OH_CO"] * OH * CO
    r8 = K["OH_CH4"] * OH * CH4
    r9 = K["HO2_NO"] * HO2 * NO
    r10 = K["OH_NO2"] * OH * NO2
    r11 = K["HO2_HO2"] * HO2 * HO2

    d = np.empty(16)
    d[0] = r_p2 - r1 - r2 - r3                       # O1D
    d[1] = r_p1 + r_p3 + r1 + r2 - r4 - r5           # O
    d[2] = 2.0 * r3 - r7 - r8 + r9 - r10             # OH
    d[3] = r7 + r8 - r9 - 2.0 * r11                  # HO2
    d[4] = -r_p2 - r_p3 + r4 - r5 - r6               # O3
    d[5] = r_p1 - r6 - r9                            # NO
    d[6] = -r_p1 + r6 + r9 - r10                     # NO2
    d[7] = -r7                                       # CO
    d[8] = -r8                                       # CH4
    d[9] = r11                                       # H2O2
    d[10] = r10                                      # HNO3
    d[11] = r7                                       # CO2
    d[12] = r_p2 + r_p3 - r4 + 2.0 * r5 + r6 + r11   # O2 (reservoir; net from stoichiometry)
    d[13] = 0.0                                      # N2 (reactant and product in r1)
    d[14] = 0.0                                      # M (reactant and product in r4)
    d[15] = -r3                                      # H2O (consumed by r3, not reproduced)
    return d


def tropospheric_reference(out_dir: Path) -> None:
    O1D0, O0, OH0, HO20 = project_radicals(N_O3_0, N_NO_0, N_NO2_0, N_CO_0, N_CH4_0, N_N2, N_O2, N_M, N_H2O)
    y0 = np.array(
        [O1D0, O0, OH0, HO20, N_O3_0, N_NO_0, N_NO2_0, N_CO_0, N_CH4_0, 0.0, 0.0, N_CO2_0, N_O2, N_N2, N_M, N_H2O]
    )

    times = output_times(1.0e-3, 1.0e5, 25)
    atol = np.maximum(1.0e-10 * np.abs(y0), 1.0e-8)
    result = solve_ivp(
        tropospheric_rhs,
        (0.0, times[-1]),
        y0,
        method="Radau",
        rtol=1.0e-12,
        atol=atol,
        t_eval=times,
    )
    if not result.success:
        raise RuntimeError(f"Tropospheric reference failed: {result.message}")
    write_csv(out_dir / "tropospheric_reference.csv", times, SPECIES, result.y)


# ---------------------------------------------------------------------------
# Chemical Akzo Nobel (CWI/Bari IVP test set): index-1 DAE of dimension 6,
#   M = diag(1,1,1,1,1,0), 0 = Ks*y1*y4 - y6.
# y6 is explicitly slaved, so the exact equivalent reduced ODE substitutes
# y6 = Ks*y1*y4 into r5 — integrable by scipy Radau (no mass-matrix support)
# with y6 reconstructed afterwards. Constants from the test-set description.
# ---------------------------------------------------------------------------

AKZO_K1, AKZO_K2, AKZO_K3, AKZO_K4 = 18.7, 0.58, 0.09, 0.42
AKZO_BIGK, AKZO_KLA, AKZO_KS = 34.4, 3.3, 115.83
AKZO_PCO2, AKZO_H = 0.9, 737.0


def akzo_reduced_rhs(_t, y):
    y1, y2, y3, y4, y5 = y
    y6 = AKZO_KS * y1 * y4
    sq2 = np.sqrt(max(y2, 0.0))
    r1 = AKZO_K1 * y1**4 * sq2
    r2 = AKZO_K2 * y3 * y4
    r3 = (AKZO_K2 / AKZO_BIGK) * y1 * y5
    r4 = AKZO_K3 * y1 * y4**2
    r5 = AKZO_K4 * y6**2 * sq2
    f_in = AKZO_KLA * (AKZO_PCO2 / AKZO_H - y2)
    return [
        -2.0 * r1 + r2 - r3 - r4,
        -0.5 * r1 - r4 - 0.5 * r5 + f_in,
        r1 - r2 + r3,
        -r2 + r3 - 2.0 * r4,
        r2 - r3 + r5,
    ]


def akzo_reference(out_dir: Path) -> None:
    y0 = [0.444, 0.00123, 0.0, 0.007, 0.0]
    times = output_times(1.0e-1, 180.0, 21)
    result = solve_ivp(
        akzo_reduced_rhs,
        (0.0, times[-1]),
        y0,
        method="Radau",
        rtol=1.0e-12,
        atol=1.0e-16,
        t_eval=times,
    )
    if not result.success:
        raise RuntimeError(f"Akzo Nobel reference failed: {result.message}")
    y6 = AKZO_KS * result.y[0] * result.y[3]
    solution = np.vstack([result.y, y6])
    write_csv(out_dir / "akzo_reference.csv", times, ["Y1", "Y2", "Y3", "Y4", "Y5", "Y6"], solution)
    print("  Akzo terminal state:", ", ".join(f"{v:.10e}" for v in solution[:, -1]))


# ---------------------------------------------------------------------------
# HIRES (IVP test set): stiff ODE of dimension 8, t_end = 321.8122.
# ---------------------------------------------------------------------------

def hires_rhs(_t, y):
    y1, y2, y3, y4, y5, y6, y7, y8 = y
    return [
        -1.71 * y1 + 0.43 * y2 + 8.32 * y3 + 0.0007,
        1.71 * y1 - 8.75 * y2,
        -10.03 * y3 + 0.43 * y4 + 0.035 * y5,
        8.32 * y2 + 1.71 * y3 - 1.12 * y4,
        -1.745 * y5 + 0.43 * y6 + 0.43 * y7,
        -280.0 * y6 * y8 + 0.69 * y4 + 1.71 * y5 - 0.43 * y6 + 0.69 * y7,
        280.0 * y6 * y8 - 1.81 * y7,
        -280.0 * y6 * y8 + 1.81 * y7,
    ]


def hires_reference(out_dir: Path) -> None:
    y0 = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0057]
    times = output_times(1.0e-1, 321.8122, 21)
    result = solve_ivp(
        hires_rhs,
        (0.0, times[-1]),
        y0,
        method="Radau",
        rtol=1.0e-12,
        atol=1.0e-16,
        t_eval=times,
    )
    if not result.success:
        raise RuntimeError(f"HIRES reference failed: {result.message}")
    write_csv(out_dir / "hires_reference.csv", times, [f"H{i}" for i in range(1, 9)], result.y)
    print("  HIRES terminal state:", ", ".join(f"{v:.10e}" for v in result.y[:, -1]))


if __name__ == "__main__":
    out_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parent / "data"
    out_dir.mkdir(parents=True, exist_ok=True)
    robertson_reference(out_dir)
    tropospheric_reference(out_dir)
    akzo_reference(out_dir)
    hires_reference(out_dir)
