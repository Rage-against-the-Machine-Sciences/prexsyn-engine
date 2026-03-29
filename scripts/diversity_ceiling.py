"""
diversity_ceiling.py

For each molecule in a testset, compute how many generated molecules (from a
processed .jsonl) fall within each Tanimoto similarity threshold.  Reports:
  - A heatmap: rows = testset molecules, cols = thresholds, values = analog count
  - At threshold=0.80: distribution stats + histogram

Usage:
    python diversity_ceiling.py \
        --jsonl   results/processed.jsonl \
        --testset data/testset.csv \
        --out-dir results/diversity_ceiling/ \
        [--n-generated 100000] \
        [--thresh-min 0.50] \
        [--thresh-max 0.95] \
        [--thresh-step 0.05] \
        [--primary-thresh 0.80] \
        [--n-workers 8]

Testset CSV must have a "smiles" column (or "target" column).
"""

from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path
from typing import Iterator

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.colors import LogNorm
from scipy import stats
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
from tqdm import tqdm


# ---------------------------------------------------------------------------
# Fingerprint helpers
# ---------------------------------------------------------------------------

def smiles_to_fp(smi: str, radius: int = 2, nbits: int = 2048):
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None
    return AllChem.GetMorganFingerprintAsBitVect(mol, radius, nbits)


def bulk_tanimoto(query_fp, fps: list) -> np.ndarray:
    """Return Tanimoto similarity of query_fp against each fp in fps."""
    sims = DataStructs.BulkTanimotoSimilarity(query_fp, fps)
    return np.array(sims, dtype=np.float32)


# ---------------------------------------------------------------------------
# Loaders
# ---------------------------------------------------------------------------

def load_testset(path: Path) -> list[str]:
    suffix = path.suffix.lower()
    if suffix in (".txt", ".smi", ".smiles"):
        raw = [l.strip().split()[0] for l in path.read_text().splitlines()
               if l.strip() and not l.startswith("#")]
    else:
        # CSV — try "smiles" then "target" column
        df = pd.read_csv(path)
        col = "smiles" if "smiles" in df.columns else "target"
        raw = df[col].dropna().tolist()

    valid = []
    for s in raw:
        mol = Chem.MolFromSmiles(s)
        if mol is not None:
            valid.append(Chem.MolToSmiles(mol, canonical=True))
    print(f"Testset: {len(valid)} valid molecules from {path}")
    return valid


def iter_jsonl(path: Path, n_max: int | None) -> Iterator[str]:
    """Yield product_smiles from processed JSONL, up to n_max records."""
    count = 0
    with open(path) as f:
        for line in f:
            if n_max is not None and count >= n_max:
                break
            line = line.strip()
            if not line:
                continue
            try:
                rec = json.loads(line)
                smi = rec.get("product_smiles", "")
                if smi:
                    yield smi
                    count += 1
            except json.JSONDecodeError:
                continue


def load_generated_fps(jsonl_path: Path, n_max: int | None) -> tuple[list[str], list]:
    """Load up to n_max generated SMILES and compute their fingerprints."""
    smiles_list, fps = [], []
    for smi in tqdm(iter_jsonl(jsonl_path, n_max), desc="loading generated", unit="mol"):
        fp = smiles_to_fp(smi)
        if fp is not None:
            smiles_list.append(smi)
            fps.append(fp)
    print(f"Loaded {len(fps)} valid generated fingerprints.")
    return smiles_list, fps


# ---------------------------------------------------------------------------
# Core computation
# ---------------------------------------------------------------------------

def compute_analog_counts(
    testset_smiles: list[str],
    generated_fps: list,
    thresholds: list[float],
) -> np.ndarray:
    """
    Returns count_matrix of shape (n_testset, n_thresholds).
    count_matrix[i, j] = number of generated molecules with
    Tanimoto >= thresholds[j] to testset molecule i.
    """
    n_test = len(testset_smiles)
    n_thresh = len(thresholds)
    thresholds_arr = np.array(thresholds, dtype=np.float32)

    count_matrix = np.zeros((n_test, n_thresh), dtype=np.int32)

    for i, smi in enumerate(tqdm(testset_smiles, desc="scoring testset", unit="mol")):
        query_fp = smiles_to_fp(smi)
        if query_fp is None:
            continue
        sims = bulk_tanimoto(query_fp, generated_fps)  # (n_generated,)
        for j, thresh in enumerate(thresholds_arr):
            count_matrix[i, j] = int(np.sum(sims >= thresh))

    return count_matrix


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def plot_heatmap(
    count_matrix: np.ndarray,
    thresholds: list[float],
    out_path: Path,
    n_generated: int,
) -> None:
    n_test, n_thresh = count_matrix.shape

    fig, ax = plt.subplots(figsize=(max(8, n_thresh * 1.2), min(20, n_test * 0.08 + 2)))
    fig.patch.set_facecolor("#f5f0e8")
    ax.set_facecolor("#f5f0e8")

    # Use log scale so sparse high-threshold columns aren't invisible
    data = count_matrix.astype(float)
    data_clipped = np.clip(data, 1, None)  # avoid log(0)

    vmin = 1
    vmax = max(data.max(), 2)

    im = ax.imshow(
        data_clipped,
        aspect="auto",
        cmap="inferno",
        norm=LogNorm(vmin=vmin, vmax=vmax),
        interpolation="nearest",
    )

    thresh_labels = [f"{t:.2f}" for t in thresholds]
    ax.set_xticks(range(n_thresh))
    ax.set_xticklabels(thresh_labels, fontsize=9)
    ax.set_xlabel("Tanimoto threshold", fontsize=11, labelpad=8)
    ax.set_ylabel("Testset molecule index", fontsize=11, labelpad=8)
    ax.set_title(
        f"Analog count heatmap  |  {n_test} targets × {n_generated:,} generated",
        # color="#ffffff",
        fontsize=13,
        pad=12,
        fontweight="bold",
    )
    # ax.tick_params(colors="#888888", length=3)
    # for spine in ax.spines.values():
    #     spine.set_edgecolor("#333333")

    cbar = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    # cbar.ax.yaxis.set_tick_params(color="#888888")
    plt.setp(cbar.ax.yaxis.get_ticklabels(), fontsize=8)
    cbar.set_label("# analogs (log scale)", fontsize=9)

    plt.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)
    print(f"Saved heatmap → {out_path}")


def plot_distribution(
    counts: np.ndarray,
    thresh: float,
    n_generated: int,
    out_path: Path,
) -> dict:
    """Plot histogram + stats for one threshold. Returns stats dict."""
    mean   = float(np.mean(counts))
    median = float(np.median(counts))
    mode_r = stats.mode(counts, keepdims=True)
    mode_v = int(mode_r.mode[0])
    mn     = int(counts.min())
    mx     = int(counts.max())
    p25    = float(np.percentile(counts, 25))
    p75    = float(np.percentile(counts, 75))
    frac_zero = float(np.mean(counts == 0))

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    fig.patch.set_facecolor("#f5f0e8")
    for ax in axes:
        ax.set_facecolor("#fdfaf4")

    # --- histogram ---
    ax = axes[0]
    nonzero = counts[counts > 0]
    bins = min(60, max(10, int(math.sqrt(len(nonzero)))))
    ax.hist(nonzero, bins=bins, color="#c0392b", edgecolor="#8b1a1a", linewidth=0.4, alpha=0.85)
    ax.axvline(mean,   color="#2c3e50", linewidth=2,   linestyle="-",  label=f"mean  {mean:.1f}")
    ax.axvline(median, color="#16a085", linewidth=2,   linestyle="--", label=f"median {median:.1f}")
    ax.axvline(mode_v, color="#8e44ad", linewidth=1.5, linestyle=":",  label=f"mode  {mode_v}")
    ax.set_xlabel("# analogs per testset molecule", fontsize=11)
    ax.set_ylabel("# testset molecules", fontsize=11)
    ax.set_title(
        f"Analog count distribution  (thresh={thresh:.2f})\n"
        f"{frac_zero*100:.1f}% targets with 0 analogs (excluded from histogram)",
        fontsize=10,
    )
    ax.legend(fontsize=9, framealpha=0.7)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # --- ECDF of non-zero counts (right panel) ---
    ax2 = axes[1]
    nonzero_counts = counts[counts > 0]
    if len(nonzero_counts) > 0:
        sorted_vals = np.sort(nonzero_counts)
        ecdf_y = np.arange(1, len(sorted_vals) + 1) / len(sorted_vals)
        ax2.step(sorted_vals, ecdf_y, color="#c0392b", linewidth=1.8, where="post")
        ax2.axvline(np.mean(nonzero_counts),   color="#2c3e50", linewidth=1.5,
                    linestyle="-",  label=f"mean  {np.mean(nonzero_counts):.1f}")
        ax2.axvline(np.median(nonzero_counts), color="#16a085", linewidth=1.5,
                    linestyle="--", label=f"median {np.median(nonzero_counts):.0f}")
        ax2.legend(fontsize=8, framealpha=0.7)
    else:
        ax2.text(0.5, 0.5, "No targets with\n>0 analogs",
                 ha="center", va="center", transform=ax2.transAxes,
                 fontsize=11, color="#888")
    ax2.set_xlabel("# analogs (non-zero targets only)", fontsize=10)
    ax2.set_ylabel("Cumulative fraction", fontsize=10)
    ax2.set_title(f"ECDF of analog counts\n(non-zero only, n={len(nonzero_counts)})", fontsize=10)
    ax2.set_ylim(0, 1.05)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)

    # annotation box — placed on histogram panel (right side), always visible
    stats_text = (
        f"n targets : {len(counts)}\n"
        f"n generated: {n_generated:,}\n"
        f"thresh     : {thresh:.2f}\n"
        f"mean       : {mean:.2f}\n"
        f"median     : {median:.1f}\n"
        f"mode       : {mode_v}\n"
        f"min / max  : {mn} / {mx}\n"
        f"Q1 / Q3    : {p25:.1f} / {p75:.1f}\n"
        f"zero-analog: {frac_zero*100:.1f}%"
    )
    axes[0].text(
        0.97, 0.97, stats_text,
        transform=axes[0].transAxes,
        fontsize=8.5,
        verticalalignment="top",
        horizontalalignment="right",
        fontfamily="monospace",
        bbox=dict(boxstyle="round,pad=0.5", facecolor="#fdfaf4", edgecolor="#aaa", alpha=0.9),
    )

    plt.suptitle(
        f"Diversity ceiling analysis  |  threshold = {thresh:.2f}",
        fontsize=13, fontweight="bold", y=1.01,
    )
    plt.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)
    print(f"Saved distribution plot → {out_path}")

    return {
        "threshold": thresh,
        "n_targets": len(counts),
        "n_generated": n_generated,
        "mean": mean,
        "median": median,
        "mode": mode_v,
        "min": mn,
        "max": mx,
        "p25": p25,
        "p75": p75,
        "frac_zero_analog": frac_zero,
    }


def plot_mean_curve(
    thresholds: list[float],
    count_matrix: np.ndarray,
    out_path: Path,
) -> None:
    """Mean analog count vs threshold — quick sanity curve."""
    means = count_matrix.mean(axis=0)
    medians = np.median(count_matrix, axis=0)

    fig, ax = plt.subplots(figsize=(8, 4))
    fig.patch.set_facecolor("#f5f0e8")
    ax.set_facecolor("#fdfaf4")

    ax.plot(thresholds, means,   "o-", color="#c0392b", linewidth=2, label="mean",   markersize=5)
    ax.plot(thresholds, medians, "s--", color="#2980b9", linewidth=1.5, label="median", markersize=4)
    ax.fill_between(
        thresholds,
        np.percentile(count_matrix, 25, axis=0),
        np.percentile(count_matrix, 75, axis=0),
        alpha=0.18, color="#c0392b", label="IQR",
    )
    ax.set_xlabel("Tanimoto threshold", fontsize=11)
    ax.set_ylabel("# analogs per target", fontsize=11)
    ax.set_title("Mean/median analog count vs threshold", fontsize=12, fontweight="bold")
    ax.legend(fontsize=9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    # Use float formatting when values are < 1 (sparse high-threshold regime)
    ymax = max(means.max(), 1e-9)
    if ymax < 1:
        ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{x:.3f}"))
    elif ymax < 10:
        ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{x:.1f}"))
    else:
        ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{int(x):,}"))
    ax.set_ylim(bottom=0)

    plt.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)
    print(f"Saved mean curve → {out_path}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(description="Diversity ceiling analysis")
    p.add_argument("--jsonl",          required=True,  help="Processed JSONL from process_batches.py")
    p.add_argument("--testset",        required=True,  help="CSV with testset SMILES (col: 'smiles' or 'target')")
    p.add_argument("--out-dir",        default="results/diversity_ceiling/", help="Output directory")
    p.add_argument("--n-generated",    type=int, default=None, help="Max generated molecules to load (default: all)")
    p.add_argument("--thresh-min",     type=float, default=0.50)
    p.add_argument("--thresh-max",     type=float, default=0.95)
    p.add_argument("--thresh-step",    type=float, default=0.05)
    p.add_argument("--primary-thresh", type=float, default=0.80, help="Threshold for detailed distribution plot")
    return p.parse_args()


def main():
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Build threshold list
    thresholds = []
    t = args.thresh_min
    while t <= args.thresh_max + 1e-9:
        thresholds.append(round(t, 4))
        t += args.thresh_step

    print(f"Thresholds: {thresholds}")

    # Load data
    testset_smiles = load_testset(Path(args.testset))
    gen_smiles, gen_fps = load_generated_fps(Path(args.jsonl), args.n_generated)
    n_generated = len(gen_fps)

    if n_generated == 0:
        print("ERROR: no valid generated fingerprints loaded.")
        return
    if not testset_smiles:
        print("ERROR: no valid testset molecules.")
        return

    # Core computation
    print(f"\nComputing analog counts: {len(testset_smiles)} targets × "
          f"{n_generated:,} generated × {len(thresholds)} thresholds ...")
    count_matrix = compute_analog_counts(testset_smiles, gen_fps, thresholds)

    # Save raw matrix
    np_out = out_dir / "count_matrix.npz"
    np.savez_compressed(
        np_out,
        count_matrix=count_matrix,
        thresholds=np.array(thresholds),
        testset_smiles=np.array(testset_smiles),
    )
    print(f"Saved count matrix → {np_out}")

    # Save summary CSV (per target, per threshold)
    thresh_cols = {f"thresh_{t:.2f}": count_matrix[:, j] for j, t in enumerate(thresholds)}
    summary_df = pd.DataFrame({"target_smiles": testset_smiles, **thresh_cols})
    summary_csv = out_dir / "analog_counts_per_target.csv"
    summary_df.to_csv(summary_csv, index=False)
    print(f"Saved per-target counts → {summary_csv}")

    # Heatmap
    plot_heatmap(
        count_matrix, thresholds,
        out_dir / "heatmap_analog_counts.png",
        n_generated,
    )

    # Mean curve
    plot_mean_curve(
        thresholds, count_matrix,
        out_dir / "mean_analog_count_vs_threshold.png",
    )

    # Detailed distribution at primary threshold
    primary_idx = min(
        range(len(thresholds)),
        key=lambda j: abs(thresholds[j] - args.primary_thresh),
    )
    primary_thresh = thresholds[primary_idx]
    primary_counts = count_matrix[:, primary_idx]

    dist_stats = plot_distribution(
        primary_counts, primary_thresh, n_generated,
        out_dir / f"distribution_thresh_{primary_thresh:.2f}.png",
    )

    # Print stats table
    print("\n" + "=" * 50)
    print(f"  STATS AT PRIMARY THRESHOLD = {primary_thresh:.2f}")
    print("=" * 50)
    for k, v in dist_stats.items():
        if isinstance(v, float):
            print(f"  {k:<22}: {v:.4f}")
        else:
            print(f"  {k:<22}: {v}")

    # Per-threshold summary table
    print("\nPer-threshold summary:")
    rows = []
    for j, t in enumerate(thresholds):
        col = count_matrix[:, j]
        rows.append({
            "threshold": t,
            "mean":      round(float(col.mean()), 2),
            "median":    round(float(np.median(col)), 2),
            "max":       int(col.max()),
            "frac_zero": round(float(np.mean(col == 0)), 4),
        })
    per_thresh_df = pd.DataFrame(rows)
    print(per_thresh_df.to_string(index=False))
    per_thresh_df.to_csv(out_dir / "per_threshold_summary.csv", index=False)

    print(f"\nAll outputs saved to {out_dir}/")


if __name__ == "__main__":
    main()
