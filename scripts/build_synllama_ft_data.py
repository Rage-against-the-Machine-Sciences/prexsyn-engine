import argparse
import json
import random

SYNLLAMA_TEMPLATE = {
    "instruction": "You are an expert synthetic organic chemist. Your task is to design a synthesis pathway for a given target molecule using common and reliable reaction templates and building blocks. Follow these instructions:\n\n1. **Input the SMILES String:** Read in the SMILES string of the target molecule and identify common reaction templates that can be applied.\n\n2. **Decompose the Target Molecule:** Use the identified reaction templates to decompose the target molecule into different intermediates.\n\n3. **Check for Building Blocks:** For each intermediate:\n   - Identify if it is a building block. If it is, wrap it in <bb> and </bb> tags and save it for later use.\n   - If it is not a building block, apply additional reaction templates to further decompose it into building blocks.\n\n4. **Document Reactions:** For each reaction documented in the output, wrap the reaction template in <rxn> and </rxn> tags.\n\n5. **Repeat the Process:** Continue this process until all intermediates are decomposed into building blocks, and document each step clearly in a structured JSON format.",
    "input": "Provide a synthetic pathway for this SMILES string: SMILES_STRING",
    "output": '{"reactions": [REACTIONS], "building_blocks": [BUILDING_BLOCKS]}',
}


def load_dataset(path):
    with open(path, "r") as f:
        obj = json.load(f)
    return obj["data"] if "data" in obj else obj


def collect_candidates(target, info):
    candidates = []

    # target
    for p in info.get("pathways", []):
        candidates.append((target, p))

    # analogs
    for a in info.get("analogs", []):
        smi = a["analog_smiles"]
        for p in a.get("pathways", []):
            candidates.append((smi, p))

    return candidates


def build_input(base_input, prev_smiles):
    if not prev_smiles:
        return base_input

    prev_str = ", ".join(prev_smiles)
    return (
        base_input
        + f"\nPreviously generated pathways target these molecules: {prev_str}. "
        + "Generate a pathway targeting a DIFFERENT, structurally distinct analog."
    )


def sample_subset(rng, pool):
    if not pool:
        return []
    k = rng.randint(0, len(pool))
    return rng.sample(pool, k) if k > 0 else []


def make_examples(rng, target, info, M=None):
    candidates = collect_candidates(target, info)
    if not candidates:
        return []
    if M is None:
        M = len(candidates)

    examples = []

    # full reconstruction sample
    recon = [c for c in candidates if c[0] == target]
    chosen = rng.choice(recon if recon else candidates)
    base_input = SYNLLAMA_TEMPLATE["input"].replace("SMILES_STRING", target)

    examples.append(
        {
            "instruction": SYNLLAMA_TEMPLATE["instruction"],
            "input": base_input,
            "output": chosen[1],
        }
    )

    # analog hints
    for _ in range(M - 1):
        out_smi, out_path = rng.choice(candidates)
        prev_pool = [s for s, _ in candidates if s != out_smi]
        prev = sample_subset(rng, prev_pool)
        inp = build_input(base_input, prev)

        examples.append(
            {
                "instruction": SYNLLAMA_TEMPLATE["instruction"],
                "input": inp,
                "output": out_path,
            }
        )

    return examples


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_json", required=True)
    parser.add_argument("--output_jsonl", required=True)
    parser.add_argument("--samples_per_target", type=int, default=None)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    rng = random.Random(args.seed)
    data = load_dataset(args.input_json)

    with open(args.output_jsonl, "w") as out:
        for target, info in data.items():
            examples = make_examples(
                rng,
                target,
                info,
                M=args.samples_per_target,
            )
            for ex in examples:
                out.write(json.dumps(ex) + "\n")


if __name__ == "__main__":
    main()
