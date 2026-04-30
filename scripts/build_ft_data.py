from pathlib import Path
import json
import argparse
from collections import defaultdict
from tqdm import tqdm

from rdkit import Chem
from rdkit.Chem import DataStructs
from rdkit.Chem.rdFingerprintGenerator import GetMorganGenerator

morgan_gen = GetMorganGenerator(radius=2, fpSize=2048)


def canonicalize_smiles(smi: str) -> str | None:
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None
    return Chem.MolToSmiles(mol, canonical=True)


def smiles_to_fp(smi: str):
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None
    return morgan_gen.GetFingerprint(mol)


def load_product_groups(jsonl_path: Path, n: int, max_paths_per_product: int):
    groups = defaultdict(list)

    with jsonl_path.open() as f:
        for line in tqdm(f, desc="loading + grouping", unit="lines"):
            line = line.strip()
            if not line:
                continue

            try:
                rec = json.loads(line)
            except json.JSONDecodeError:
                continue

            smi = canonicalize_smiles(rec.get("product_smiles", ""))
            if smi is None:
                continue

            pathway = rec.get("json_string")
            if not isinstance(pathway, str):
                continue
            assert pathway.strip().startswith("{")

            if pathway is None:
                continue

            # cap pathways per product
            if len(groups[smi]) < max_paths_per_product:
                groups[smi].append(pathway)

            # stop once we have enough unique products
            if len(groups) >= n:
                break

    return groups


def build_analog_groups(groups, threshold: float, min_analogs: int):
    smiles = list(groups.keys())
    fps = [smiles_to_fp(s) for s in smiles]

    neighbors = {s: set() for s in smiles}

    for i, fp_i in enumerate(tqdm(fps, desc="building analog graph", unit="mol")):
        if fp_i is None:
            continue

        sims = DataStructs.BulkTanimotoSimilarity(fp_i, fps[i + 1 :])
        for offset, sim in enumerate(sims, start=1):
            if sim >= threshold:
                j = i + offset
                si = smiles[i]
                sj = smiles[j]
                neighbors[si].add(sj)
                neighbors[sj].add(si)

    # debug stats
    degrees = [len(v) for v in neighbors.values()]
    print("max degree:", max(degrees))
    print("num nodes with >=1:", sum(d >= 1 for d in degrees))
    print("num nodes with >=2:", sum(d >= 2 for d in degrees))
    print("num nodes with >=3:", sum(d >= 3 for d in degrees))

    dataset = {}
    for smi in smiles:
        analogs = sorted(neighbors[smi])
        if len(analogs) < min_analogs:
            continue

        dataset[smi] = {
            # multiple pathways
            "pathways": groups[smi],
            "analogs": [
                {
                    "analog_smiles": a,
                    "pathways": groups[a],
                }
                for a in analogs
            ],
        }

    return dataset, neighbors


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--n", type=int, required=True)
    p.add_argument("--threshold", type=float, default=0.8)
    p.add_argument("--min-analogs", type=int, default=3)
    p.add_argument("--max-paths-per-product", type=int, default=3)
    return p.parse_args()


def main():
    args = parse_args()

    in_path = Path(args.input)
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    groups = load_product_groups(
        in_path,
        args.n,
        args.max_paths_per_product,
    )

    if not groups:
        raise SystemExit("No valid data.")

    dataset, neighbors = build_analog_groups(
        groups,
        args.threshold,
        args.min_analogs,
    )

    payload = {
        "meta": {
            "input": str(in_path),
            "n_products": len(groups),
            "threshold": args.threshold,
            "min_analogs": args.min_analogs,
            "max_paths_per_product": args.max_paths_per_product,
            "n_targets_kept": len(dataset),
        },
        "data": dataset,
    }

    with out_path.open("w") as f:
        json.dump(payload, f, indent=2)

    n_edges = sum(len(v) for v in neighbors.values()) // 2
    print(f"unique products        : {len(groups)}")
    print(f"analog pairs           : {n_edges}")
    print(f"targets kept           : {len(dataset)}")
    print(f"written                : {out_path}")


if __name__ == "__main__":
    main()
