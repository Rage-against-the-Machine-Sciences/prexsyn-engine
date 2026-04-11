"""
plot_mean_curves.py

Plot mean analog count vs threshold for all testset subdirs on one figure.

Usage:
    python plot_mean_curves.py results/diversity_ceiling_100k/
    python plot_mean_curves.py results/diversity_ceiling_1M/
"""

import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker


def main(root_dir: Path) -> None:
    subdirs = sorted([d for d in root_dir.iterdir() if d.is_dir()])
    if not subdirs:
        print(f"No subdirectories found in {root_dir}")
        return

    fig, ax = plt.subplots(figsize=(10, 5))
    fig.patch.set_facecolor("#f5f0e8")
    ax.set_facecolor("#fdfaf4")

    colors = plt.cm.tab10.colors

    for i, subdir in enumerate(subdirs):
        npz_path = subdir / "count_matrix.npz"
        if not npz_path.exists():
            print(f"  skipping {subdir.name} — no count_matrix.npz")
            continue

        data = np.load(npz_path)
        count_matrix = data["count_matrix"]   # (n_targets, n_thresholds)
        thresholds   = data["thresholds"].tolist()
        means        = count_matrix.mean(axis=0)

        label = subdir.name.replace("1k_", "").replace("_", " ")
        color = colors[i % len(colors)]

        ax.plot(thresholds, means, "o-", color=color, linewidth=2,
                markersize=5, label=label)

    ax.set_xlabel("Tanimoto threshold", fontsize=12)
    ax.set_ylabel("Mean # analogs per target", fontsize=12)
    ax.set_title(
        f"Mean analog count vs threshold — {root_dir.name}",
        fontsize=13, fontweight="bold",
    )
    ax.legend(fontsize=9, framealpha=0.8, loc="upper right")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ymax = ax.get_ylim()[1]
    if ymax < 1:
        ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{x:.3f}"))
    elif ymax < 10:
        ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{x:.1f}"))
    else:
        ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{int(x):,}"))

    ax.set_ylim(bottom=0)
    ax.set_xticks([round(t, 2) for t in thresholds])

    plt.tight_layout()
    out_path = root_dir / "mean_analog_curves_all_testsets.png"
    fig.savefig(out_path, dpi=150, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)
    print(f"Saved → {out_path}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python plot_mean_curves.py <results_dir>")
        sys.exit(1)
    main(Path(sys.argv[1]))
