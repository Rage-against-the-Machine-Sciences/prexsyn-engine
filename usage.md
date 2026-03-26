# to generate training data


# install the thing

install pixi

for runpod specific setup, add `export PATH="/workspace/.pixi/bin:$PATH"` to `/root/.bashrc` and source it so install persists on volume.

Build the package:
1. `pixi shell`
2. `pixi run configure-release` (will take long)
3. `pixi run build-release`
4. `pixi run pyinstall` (will take long)


Get Enamine US Stock (Oct. 2023) version from SynFormer (can also use latest version's SDF directly after requesting from [Enamine](https://enamine.net/building-blocks/building-blocks-catalog)):
1. `git clone https://github.com/wenhao-gao/synformer.git`
2. Download the processed fingerprints of the BBs: `curl -L -o fpindex.pkl "https://huggingface.co/whgao/synformer/resolve/main/fpindex.pkl?download=true"`
3. `mv fpindex.pkl data/`
4. Generate raw SDF with IDs:
```bash
cd /workspace/synformer
python -c "
import pickle
from rdkit import Chem

with open('data/processed/comp_2048/fpindex.pkl', 'rb') as f:
    fpindex = pickle.load(f)

mols = fpindex._molecules
print(f'Loaded {len(mols)} BBs')

writer = Chem.SDWriter('/workspace/prexsyn-engine/data/Enamine_Rush-Delivery_Building_Blocks-US_223244cmpd_20231001.sdf')
skipped = 0
for i, mol in enumerate(mols):
    rdmol = Chem.MolFromSmiles(mol.smiles)
    if rdmol is None:
        skipped += 1
        continue
    # Set molecule name as the 'id' property — prexsyn reads _Name as id
    rdmol.SetProp('_Name', f'BB{i}')
    rdmol.SetProp('id', f'BB{i}')
    writer.write(rdmol)
writer.close()
print(f'Written {len(mols) - skipped}, skipped {skipped}')
"
```
5. Switch to prexsyn_engine repo and shell
6. Run chemical space builder: `python -m scripts.enamine_rxn115`
6. Verify Chemical Space (should show counts: 211220 / 115 / 317193.):
```bash
# Quick sanity check
python -c "
import prexsyn_engine
cs = prexsyn_engine.chemspace.ChemicalSpace.deserialize('data/enamine_rxn115_chemspace.bin')
print('BBs:', len(cs.bb_lib()))
print('Reactions:', len(cs.rxn_lib()))
print('Intermediates:', len(cs.int_lib()))
"
```
7. To verify, run datapipe_dryrun script to verify in inspect mode: `python -m scripts.datapipe_dryrun --inspect --cd-path ./data/enamine_rxn115_chemspace.bin`
8. Run modified data generation script for our requirements.
