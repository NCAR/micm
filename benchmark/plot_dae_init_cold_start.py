#!/usr/bin/env python3
"""Turn dae_init_cold_start's CSV output into the figures for D1, D2 and D4.

Usage:
    python3 plot_dae_init_cold_start.py <input_dir> <output_dir>

<input_dir> holds dae_init_basin.csv and dae_init_cost.csv as written by the
benchmark, and optionally dae_init_trace.csv. Figures are written to
<output_dir> as PNG at 150 dpi.

Deliberately dependency-light: matplotlib and the standard library only, no
pandas. This is the first plotting script in the repository, so it sets the
convention rather than following one.
"""

import csv
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import ListedColormap  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

# Variant order is the story order: what main did, what PR #1083 did, what the
# line search does. "undamped" is PR #1083 with the search disabled, which the
# opt-out test pins as bit-identical to PR #1083 itself.
VARIANTS = ["main", "undamped", "damped"]
VARIANT_TITLES = {
    "main": "origin/main",
    "undamped": "PR #1083 (undamped)",
    "damped": "PR #1083 + line search",
}

# Converged green, everything else a warm colour: the eye should find failures.
STATUS_COLORS = {
    "Converged": "#2a9d5c",
    "Constraint Initialization Failed": "#c1443c",
    "NaN Detected": "#d97b29",
    "Inf Detected": "#8c5bb5",
}
UNKNOWN_COLOR = "#7a7a7a"


def read_csv(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle))


def plot_basin(rows, out_path):
    present = [v for v in VARIANTS if any(r["variant"] == v for r in rows)]
    if not present:
        return None

    temperatures = sorted({float(r["temperature"]) for r in rows})
    guesses = sorted({float(r["initial_xm"]) for r in rows})
    statuses = sorted({r["status"] for r in rows})
    palette = [STATUS_COLORS.get(s, UNKNOWN_COLOR) for s in statuses]
    index = {s: i for i, s in enumerate(statuses)}

    fig, axes = plt.subplots(
        1, len(present), figsize=(4.6 * len(present), 4.4), sharex=True, sharey=True, squeeze=False
    )
    for ax, variant in zip(axes[0], present):
        cell = {
            (float(r["temperature"]), float(r["initial_xm"])): index[r["status"]]
            for r in rows
            if r["variant"] == variant
        }
        grid = [[cell.get((t, g), -1) for t in temperatures] for g in guesses]
        ax.imshow(
            grid,
            aspect="auto",
            origin="lower",
            interpolation="nearest",
            cmap=ListedColormap(palette),
            vmin=-0.5,
            vmax=len(statuses) - 0.5,
        )
        ax.set_xticks(range(len(temperatures)))
        ax.set_xticklabels([f"{t:.0f}" for t in temperatures], fontsize=7, rotation=90)
        ax.set_yticks(range(len(guesses)))
        ax.set_yticklabels([f"{g:g}" for g in guesses], fontsize=7)
        ax.set_title(VARIANT_TITLES.get(variant, variant), fontsize=10)
        ax.set_xlabel("temperature (K)", fontsize=9)
    axes[0][0].set_ylabel("initial [XM] guess", fontsize=9)

    fig.legend(
        handles=[Patch(facecolor=STATUS_COLORS.get(s, UNKNOWN_COLOR), label=s) for s in statuses],
        loc="lower center",
        ncol=min(len(statuses), 4),
        fontsize=8,
        frameon=False,
    )
    fig.suptitle("DAE constraint-initialization basin: consistent-IC projection outcome", fontsize=11)
    fig.tight_layout(rect=(0, 0.10, 1, 0.96))
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    return out_path


def plot_cost(rows, out_path):
    starts = sorted({r["start"] for r in rows})
    variants = [v for v in VARIANTS if any(r["variant"] == v for r in rows)]
    metrics = [
        ("constraint_init_iterations", "Newton updates"),
        ("solves", "triangular solves"),
        ("median_us", "median wall time (us)"),
    ]

    fig, axes = plt.subplots(1, len(metrics), figsize=(4.2 * len(metrics), 3.8), squeeze=False)
    width = 0.8 / max(len(variants), 1)
    for ax, (field, label) in zip(axes[0], metrics):
        for i, variant in enumerate(variants):
            values = []
            for start in starts:
                match = [r for r in rows if r["variant"] == variant and r["start"] == start]
                values.append(float(match[0][field]) if match else 0.0)
            offset = (i - (len(variants) - 1) / 2) * width
            bars = ax.bar([x + offset for x in range(len(starts))], values, width, label=VARIANT_TITLES.get(variant, variant))
            ax.bar_label(bars, fmt="%.3g", fontsize=7, padding=1)
        ax.set_xticks(range(len(starts)))
        ax.set_xticklabels(starts, fontsize=9)
        ax.set_title(label, fontsize=10)
        ax.margins(y=0.18)
    axes[0][0].set_ylabel("per projection", fontsize=9)
    axes[0][-1].legend(fontsize=8, frameon=False)
    fig.suptitle("Constraint-initialization cost: the line search is free on the consistent path", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    return out_path


def plot_trace(rows, out_path):
    if not rows:
        return None
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    for variant in VARIANTS:
        series = [r for r in rows if r["variant"] == variant]
        if not series:
            continue
        series.sort(key=lambda r: int(r["update"]))
        updates = [int(r["update"]) for r in series]
        norms = [float(r["weighted_correction_norm"]) for r in series]
        ax.semilogy(updates, norms, marker="o", markersize=4, label=VARIANT_TITLES.get(variant, variant))
        for r in series:
            step = float(r["accepted_step"])
            # The last update of a run accepts no step - it either converged or exhausted the
            # budget - and is recorded as NaN. Annotate only genuinely damped steps.
            if step == step and step != 1.0:
                ax.annotate(
                    f"lambda={step:g}",
                    (int(r["update"]), float(r["weighted_correction_norm"])),
                    textcoords="offset points",
                    xytext=(4, 6),
                    fontsize=7,
                )
    ax.set_xlabel("Newton update", fontsize=9)
    ax.set_ylabel(r"$\|\delta\|_w$ (weighted Newton correction)", fontsize=9)
    ax.set_title("Cold-start convergence trace, T = 284 K, initial [XM] = 1", fontsize=11)
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(fontsize=8, frameon=False)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    return out_path


def main():
    if len(sys.argv) != 3:
        print(__doc__)
        return 1
    in_dir, out_dir = sys.argv[1], sys.argv[2]
    os.makedirs(out_dir, exist_ok=True)

    written = []
    basin = os.path.join(in_dir, "dae_init_basin.csv")
    if os.path.exists(basin):
        written.append(plot_basin(read_csv(basin), os.path.join(out_dir, "dae_init_basin.png")))
    cost = os.path.join(in_dir, "dae_init_cost.csv")
    if os.path.exists(cost):
        written.append(plot_cost(read_csv(cost), os.path.join(out_dir, "dae_init_cost.png")))
    trace = os.path.join(in_dir, "dae_init_trace.csv")
    if os.path.exists(trace):
        written.append(plot_trace(read_csv(trace), os.path.join(out_dir, "dae_init_trace.png")))
    else:
        print("note: dae_init_trace.csv absent, skipping the trace figure")

    for path in [w for w in written if w]:
        print(f"wrote {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
