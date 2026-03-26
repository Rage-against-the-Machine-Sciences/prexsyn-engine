from pathlib import Path
from typing import Iterator
import json

import click
import numpy as np
from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem import AllChem

import prexsyn_engine
from prexsyn_engine.chemspace import ChemicalSpace
from prexsyn_engine.chemistry import Molecule, Reaction

TOK_PAD = 0
TOK_END = 1
TOK_START = 2
TOK_BB = 3
TOK_RXN = 4


def decode_pathway(
    sample_id: int,
    synthesis_tokens: np.ndarray,
    chemspace: ChemicalSpace,
    rxn_smarts: list[str],
) -> dict | None:
    """
    Decode a single (16, 3) synthesis token sequence into all representations.

    Returns a dict with keys:
        product_smiles      : canonical SMILES of the final product
        bb_smiles           : list of BB SMILES in order of appearance (BU order)
        rxn_indices         : list of reaction indices in order of appearance
        pathway_flat_bu     : "BB1;BB2;R{idx};..." semicolon-joined, BU order
                              (compatible with diversity_eval.py)
        pathway_flat_td     : same tokens reversed (TD order, for ReaSyn TD)
        json_string         : SynLlama JSON format
        synthesis_indices   : raw token triplets "type,idx,aux;..." for PrexSyn FT

    Returns None if the pathway is invalid (reaction fails, empty stack, etc.).
    """
    bb_lib = chemspace.bb_lib()
    rxn_lib = chemspace.rxn_lib()

    # ---- stack-based execution (BU order) ----
    stack: list[str] = []  # SMILES strings of current intermediates
    bb_smiles: list[str] = []  # BBs in appearance order
    rxn_indices: list[int] = []  # reaction indices in appearance order
    rxn_int_smiles: list[str] = []  # intermediate SMILES in order of rxn appearance
    flat_tokens: list[str] = []  # for pathway_flat_bu
    raw_triplets: list[str] = []  # for synthesis_indices
    rxn_steps = []  # for synllama

    for frame in synthesis_tokens:
        tok_type, idx, aux = int(frame[0]), int(frame[1]), int(frame[2])

        if tok_type == TOK_PAD:
            break
        if tok_type == TOK_END:
            break
        if tok_type == TOK_START:
            continue

        elif tok_type == TOK_BB:
            bb_item = bb_lib.get(idx)
            smi = bb_item.molecule.smiles()
            stack.append(smi)
            bb_smiles.append(smi)
            flat_tokens.append(smi)
            raw_triplets.append(f"{tok_type},{idx},{aux}")

        elif tok_type == TOK_RXN:
            rxn_item = rxn_lib.get(aux)  # aux column holds rxn index
            rxn_sma = rxn_smarts[aux]
            n_reactants = rxn_item.reaction.num_reactants()
            rxn_indices.append(aux)
            flat_tokens.append(rxn_sma)
            raw_triplets.append(f"{tok_type},{idx},{aux}")

            if len(stack) < n_reactants:
                print(
                    f"Malformed sequence: len(stack) = {len(stack)} < n_reactants = {n_reactants}"
                )
                return None

            reactant_smiles = list(stack[-n_reactants:])
            stack = stack[:-n_reactants]

            # prexsyn api to apply reaction and get heaviest product
            try:
                rxn = Reaction.from_smarts(rxn_sma)
                reactants = list(map(Molecule.from_smiles, reactant_smiles))
                outcomes = rxn.apply(reactants, ignore_errors=True)
                if not outcomes:
                    return None  # invalid pathway

                # TODO: should we consider all outcomes?
                prod = outcomes[0].main_product()
                if prod is None:
                    return None
                prod_smi = prod.smiles()

                rxn_int_smiles.append(prod_smi)
                stack.append(prod_smi)
                rxn_steps.append(
                    {
                        "reaction_template": rxn_sma,
                        "rxn_index": aux,
                        "reactants": reactant_smiles,
                        "product": prod_smi,
                    }
                )
            except Exception:
                print(f"Error while applying rxn `{rxn_sma}` to: {reactant_smiles}")
                return None

    if len(stack) != 1:
        return None

    product_smiles = stack[0]

    pathway_flat_bu = ";".join(flat_tokens)
    pathway_flat_td = ";".join(reversed(flat_tokens))
    synthesis_indices = ";".join(raw_triplets)

    json_string = build_synllama_json(rxn_steps)

    return {
        "sample_id": sample_id,
        "product_smiles": product_smiles,
        "bb_smiles": bb_smiles,
        "rxn_indices": rxn_indices,
        "rxn_int_smiles": rxn_int_smiles,
        "pathway_flat_bu": pathway_flat_bu,
        "pathway_flat_td": pathway_flat_td,
        "json_string": json_string,
        "synthesis_indices": synthesis_indices,
    }


def build_synllama_json(
    rxn_steps: list[dict],
) -> str:
    """
    Build SynLlama-style JSON.

    Expected rxn_steps order:
        bottom-up / postfix order
        [
          {
            "reaction_template": str,   # SMARTS
            "rxn_index": int,
            "reactants": list[str],     # direct reactant SMILES, in application order
            "product": str,             # chosen single product SMILES
          },
          ...
        ]

    Output:
        {
          "reactions": [
            {
              "reaction_number": 1,
              "reaction_template": "...",
              "reactants": [...],
              "product": "..."
            },
            ...
          ],
          "building_blocks": [...]
        }
    """

    # final step first, leafward steps later.
    reactions_out = []
    for reaction_number, step in enumerate(reversed(rxn_steps), start=1):
        reactions_out.append(
            {
                "reaction_number": reaction_number,
                "reaction_template": step["reaction_template"],
                "reactants": step["reactants"],
                "product": step["product"],
            }
        )

    # BBs: reactants that are never products of any other step.
    produced_products = {step["product"] for step in rxn_steps if step["product"]}
    building_blocks = []
    seen_bb = set()

    for step in rxn_steps:
        for smi in step["reactants"]:
            if smi not in produced_products and smi not in seen_bb:
                seen_bb.add(smi)
                building_blocks.append(smi)

    payload = {
        "reactions": reactions_out,
        "building_blocks": building_blocks,
    }
    return json.dumps(payload, ensure_ascii=False)


def iter_batches(batch_dir: Path) -> Iterator[tuple[int, np.ndarray]]:
    npz_files = sorted(batch_dir.glob("batch_*.npz"))
    if not npz_files:
        raise FileNotFoundError(f"No batch_*.npz files found in {batch_dir}")
    for npz_path in npz_files:
        batch_idx = int(npz_path.stem.split("_")[1])
        data = np.load(npz_path)
        yield batch_idx, data["synthesis"]  # (1024, 16, 3) int64


@click.command()
@click.option(
    "--cs-path",
    type=str,
    default="./data/enamine_rxn115_chemspace.bin",
    help="Path to chemical space .bin file",
)
@click.option(
    "--smarts-path",
    type=str,
    default="./data/rxn115.txt",
    help="Path to the file containing the SMARTS for the 115 reaction templates set.",
)
@click.option(
    "--batch-dir",
    type=str,
    required=True,
    help="Directory containing batch_*.npz files from datapipe_save_batches.py",
)
@click.option(
    "--output",
    type=str,
    required=True,
    help="Output JSONL file path",
)
@click.option(
    "--max-batches",
    type=int,
    default=None,
    help="Stop after this many batches (default: process all)",
)
def main(cs_path, smarts_path, batch_dir, output, max_batches):
    cs_path = Path(cs_path)
    smarts_path = Path(smarts_path)
    batch_dir = Path(batch_dir)
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)

    print(f"Loading chemical space from {cs_path} ...")
    chemspace = prexsyn_engine.chemspace.ChemicalSpace.deserialize(str(cs_path))
    print(
        f"  BBs: {len(chemspace.bb_lib())}  "
        f"RXNs: {len(chemspace.rxn_lib())}  "
        f"Ints: {len(chemspace.int_lib())}"
    )

    rxn_smarts = smarts_path.read_text().splitlines()
    print(f"Loaded {len(rxn_smarts)} reaction template SMARTS from {smarts_path}.")

    total_samples = 0
    total_decoded = 0
    total_written = 0
    total_failed = 0

    batch_iter = iter_batches(batch_dir)

    with open(output, "w") as out_f:
        for batch_count, (batch_idx, synthesis) in enumerate(
            tqdm(batch_iter, desc="batches")
        ):
            if max_batches is not None and batch_count >= max_batches:
                break

            batch_size = synthesis.shape[0]

            for sample_offset in range(batch_size):
                sample_id = batch_idx * 1024 + sample_offset
                tokens = synthesis[sample_offset]  # (16, 3)

                total_samples += 1

                decoded = decode_pathway(sample_id, tokens, chemspace, rxn_smarts)
                if decoded is None:
                    total_failed += 1
                    continue

                total_decoded += 1

                out_f.write(json.dumps(decoded) + "\n")
                total_written += 1

    print(f"\nDone.")
    print(f"  Total samples : {total_samples}")
    print(
        f"  Decoded OK    : {total_decoded}  "
        f"({100*total_decoded/max(total_samples,1):.1f}%)"
    )
    print(f"  Written       : {total_written}")
    print(f"  Failed/invalid: {total_failed}")


if __name__ == "__main__":
    main()
