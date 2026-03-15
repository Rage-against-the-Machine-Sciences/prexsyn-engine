import pickle

import pytest
import prexsyn_engine


chemistry = prexsyn_engine.chemistry
Molecule = chemistry.Molecule
MoleculeError = chemistry.MoleculeError


def test_from_smiles_success():
    mol = Molecule.from_smiles("CCO")

    assert mol.smiles() == "CCO"
    assert mol.num_heavy_atoms() == 3


def test_from_smiles_invalid_raises_molecule_error():
    with pytest.raises(MoleculeError, match="Failed to parse SMILES"):
        Molecule.from_smiles("not-a-smiles")


def test_repr_contains_smiles():
    mol = Molecule.from_smiles("CO")

    assert repr(mol) == "<Molecule CO>"


def test_pickle_roundtrip_preserves_smiles_and_properties():
    mol = Molecule.from_smiles("CCN")

    dumped = pickle.dumps(mol)
    loaded = pickle.loads(dumped)

    assert loaded.smiles() == "CCN"
    assert loaded.num_heavy_atoms() == 3


def test_largest_fragment_returns_main_component():
    # Ethane + water separated by dot notation; ethane is the largest fragment.
    mol = Molecule.from_smiles("CC.O")

    largest = mol.largest_fragment()

    assert largest.smiles() == "CC"
    assert largest.num_heavy_atoms() == 2
