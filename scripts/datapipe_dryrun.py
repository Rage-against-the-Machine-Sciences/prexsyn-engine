from pathlib import Path

import click
from tqdm import tqdm

import prexsyn_engine


@click.command()
@click.option(
    "--cs-path",
    type=str,
    required=True,
    default="./data/enamine_rxn115_chemspace.bin",
    help="Path to chemical space cache file",
)
def main(cs_path):
    cs_path = Path(cs_path)
    if not cs_path.exists():
        print(f"Chemical space cache file not found at {cs_path}")
        return
    print(f"Loading chemical space from cache at {cs_path}")
    chemspace = prexsyn_engine.chemspace.ChemicalSpace.deserialize(cs_path)

    print(f"Chemical space at {cs_path} contains:")
    print("- Number of building blocks:", len(chemspace.bb_lib()))
    print("- Number of reactions:", len(chemspace.rxn_lib()))
    print("- Number of intermediates:", len(chemspace.int_lib()))

    datapipe_line = prexsyn_engine.datapipe.DataPipeline(
        chemspace,
        {
            "ecfp4": prexsyn_engine.descriptor.MorganFingerprint.ecfp4(),
            "fcfp4": prexsyn_engine.descriptor.MorganFingerprint.fcfp4(),
        },
        {"synthesis": prexsyn_engine.descriptor.SynthesisPostfixNotation.create(16)},
    )

    try:
        datapipe_line.start_workers(list(range(16)))
        for i in tqdm(range(1000)):
            batch = datapipe_line.get(1024)
            del batch
        input("Data pipeline workers started. Press Enter to stop...")
    except KeyboardInterrupt:
        print("Keyboard interrupt received. Stopping workers...")
    finally:
        datapipe_line.stop_workers()


if __name__ == "__main__":
    main()
