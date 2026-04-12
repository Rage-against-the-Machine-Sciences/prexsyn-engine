import numpy as np
import json
from collections import Counter
from pathlib import Path
from tqdm import tqdm
from rdkit import Chem
import argparse
import matplotlib.pyplot as plt
import pickle
import prexsyn_engine as pe


def gini(x):
    x = np.array(x, dtype=np.float64)
    if np.amin(x) < 0:
        x -= np.amin(x)
    x += 1e-9
    x = np.sort(x)
    n = x.size
    return (np.sum((2 * np.arange(1, n + 1) - n - 1) * x)) / (n * np.sum(x))


def canonical(smi):
    # canonical form
    # drop stereochemistry (chirality and double bond geometry)
    mol = Chem.MolFromSmiles(smi)
    return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=False) if mol else None


def plot_frequency(vals, out_path="freq_plot.png", log: bool = False):
    vals_sorted = np.sort(vals)[::-1]  # descending
    ranks = np.arange(1, len(vals_sorted) + 1)

    plt.figure()
    plt.plot(ranks, vals_sorted)
    plt.xlabel("BB rank (sorted by frequency)")
    plt.ylabel("Frequency (count)")
    plt.title("BB Usage Frequency Distribution")
    if log:
        plt.yscale("log")
        plt.xscale("log")
    plt.tight_layout()
    plt.savefig(out_path)
    plt.close()


def plot_lorenz(vals, out_path="lorenz_curve.png"):
    vals = np.sort(vals)
    cumvals = np.cumsum(vals)
    cumvals = cumvals / cumvals[-1]

    x = np.linspace(0, 1, len(vals))

    plt.figure()
    plt.plot(x, cumvals, label="Lorenz curve")
    plt.plot([0, 1], [0, 1], linestyle="--", label="Uniform")
    plt.xlabel("Fraction of BBs")
    plt.ylabel("Fraction of total usage")
    plt.title("Lorenz Curve of BB Usage")
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_path)
    plt.close()


def plot_hist(vals, out_path="hist.png", log: bool = False):
    plt.figure()
    plt.hist(vals, bins=100, log=log)
    plt.xlabel("Frequency")
    plt.ylabel("Number of BBs")
    plt.title("Histogram of BB Frequencies")
    plt.tight_layout()
    plt.savefig(out_path)
    plt.close()


def get_counts_and_vis(input_path: Path, num_samples: int, lib_size: int) -> Counter:
    counts = Counter()

    with open(input_path) as f:
        for line in tqdm(f, desc="Processing"):
            rec = json.loads(line)
            for smi in rec["bb_smiles"]:
                c = canonical(smi)
                if c:
                    counts[c] += 1

    print(f"Unique BBs seen: {len(counts)}")
    print(f"Top 10: {counts.most_common(10)}")
    print(f"Total occurrences: {sum(counts.values())}")

    vals = np.array(list(counts.values()))
    print(f"% of library seen: {len(counts)/lib_size*100:.1f}%")
    print(f"Gini coeff: {gini(vals):.3f}")

    out_dir = input_path.parent / f"train_distribution_{num_samples}"
    out_dir.mkdir(parents=True, exist_ok=True)

    plot_frequency(vals, out_path=out_dir / "freq_nonlog.png", log=False)
    plot_frequency(vals, out_path=out_dir / "freq_log.png", log=True)
    plot_lorenz(vals, out_path=out_dir / "lorenz_curve.png")
    plot_hist(vals, out_path=out_dir / "hist_nonlog.png", log=False)
    plot_hist(vals, out_path=out_dir / "hist_log.png", log=True)

    return counts


def build_bb_distribution(
    counts: Counter, cs: pe.chemspace.ChemicalSpace, outfile: Path
):
    # iterate bb_lib to build smiles -> index map
    bb_lib = cs.bb_lib()
    n = len(bb_lib)

    smiles_to_idx = {}
    for i in range(n):
        item = bb_lib[i]
        smi = canonical(item.molecule.smiles)
        if smi:
            smiles_to_idx[smi] = i

    print(f"Built smiles -> index map: {len(smiles_to_idx)} entries")

    # Build distribution vector
    weights = np.zeros(n, dtype=np.float64)
    matched = 0
    for smi, count in counts.items():
        idx = smiles_to_idx.get(smi)
        if idx is not None:
            weights[idx] = float(count)
            matched += 1

    print(f"Matched {matched}/{len(counts)} corpus BBs to library indices")
    print(f"Non-zero entries: {np.count_nonzero(weights)}/{n}")

    np.save(outfile, weights)
    print(f"Saved weights to {outfile}")


def main():
    parser = argparse.ArgumentParser(description="Analyze BB distribution from JSONL")
    parser.add_argument(
        "--input",
        type=str,
        required=True,
        help="Path to input JSONL file",
    )
    parser.add_argument(
        "--nsamples",
        type=int,
        help="Number of samples in jsonl",
    )
    parser.add_argument(
        "--chemspace",
        required=True,
        help="Path to .bin chemspace file",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output .npy path",
    )

    args = parser.parse_args()
    input_path = Path(args.input)
    outfile = Path(args.output)
    num_samples = int(args.nsamples)
    cs = pe.chemspace.ChemicalSpace.deserialize(args.chemspace)
    lib_size = cs.bb_lib().size()

    counts = get_counts_and_vis(input_path, num_samples, lib_size)
    build_bb_distribution(counts, cs, outfile)


if __name__ == "__main__":
    main()
