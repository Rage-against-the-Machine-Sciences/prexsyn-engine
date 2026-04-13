from pathlib import Path

import click
import numpy as np
from tqdm import tqdm

import prexsyn_engine

BATCH_SIZE = 1024


@click.command()
@click.option(
    "--cs-path",
    type=str,
    required=True,
    default="./data/enamine_rxn115_chemspace.bin",
    help="Path to chemical space cache file",
)
@click.option(
    "--num-samples",
    type=int,
    required=True,
    default=1000,
    help="Number of samples to generate",
)
@click.option(
    "--out-dir",
    type=str,
    required=True,
    default="./results/samples/",
    help="Path to store samples in batches",
)
@click.option(
    "--bb-weights",
    type=str,
    default=None,
    help="Path to bb_dist .npy weights file (optional, defaults to uniform)",
)
def main(cs_path, num_samples, out_dir, bb_weights):
    cs_path = Path(cs_path)
    if not cs_path.exists():
        print(f"Chemical space cache file not found at {cs_path}")
        return
    print(f"Loading chemical space from cache at {cs_path}")
    chemspace = prexsyn_engine.chemspace.ChemicalSpace.deserialize(cs_path)

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"Chemical space at {cs_path} contains:")
    print("- Number of building blocks:", len(chemspace.bb_lib()))
    print("- Number of reactions:", len(chemspace.rxn_lib()))
    print("- Number of intermediates:", len(chemspace.int_lib()))

    weights = []
    if bb_weights is not None:
        import numpy as np

        weights = np.load(bb_weights).tolist()
        print(
            f"Loaded BB weights from {bb_weights}, non-zero: {sum(1 for w in weights if w > 0)}"
        )

    datapipe_line = prexsyn_engine.datapipe.DataPipeline(
        chemspace,
        {
            "ecfp4": prexsyn_engine.descriptor.MorganFingerprint.ecfp4(),
            "fcfp4": prexsyn_engine.descriptor.MorganFingerprint.fcfp4(),
        },
        {"synthesis": prexsyn_engine.descriptor.SynthesisPostfixNotation.create(16)},
        bb_weights=weights,  # empty list -> uniform
        smoothing_alpha=1.0,
    )

    num_batches = int(np.ceil(num_samples / BATCH_SIZE))

    datapipe_line.start_workers(list(range(16)))
    try:
        for i in tqdm(range(num_batches)):
            batch = datapipe_line.get(BATCH_SIZE)
            np.savez_compressed(
                out_dir / f"batch_{i:06d}.npz",
                synthesis=batch["synthesis"],  # (1024, 16, 3) int64
                ecfp4=batch["ecfp4"],  # (1024, 2048) bool
            )
            del batch
    except KeyboardInterrupt:
        pass
    finally:
        datapipe_line.stop_workers()


if __name__ == "__main__":
    main()
