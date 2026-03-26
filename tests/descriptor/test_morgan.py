from pathlib import Path

import numpy as np

from prexsyn_engine import chemistry, descriptor


Molecule = chemistry.Molecule


def _reference_csv_path() -> Path:
    return Path(__file__).resolve().parents[2] / "resources" / "test" / "fp_ref.csv"


def _parse_indices(field: str) -> list[int]:
    if not field:
        return []
    return [int(token) for token in field.split(";") if token]


def _parse_reference_line(line: str) -> tuple[str, list[int], list[int]]:
    first_comma = line.find(",")
    second_comma = line.find(",", first_comma + 1)
    third_comma = line.find(",", second_comma + 1)

    if first_comma == -1 or second_comma == -1 or third_comma == -1:
        raise ValueError("Malformed fingerprint reference CSV row")

    smiles = line[:first_comma]
    ecfp4_field = line[first_comma + 1 : second_comma]
    fcfp4_field = line[second_comma + 1 : third_comma]
    return smiles, _parse_indices(ecfp4_field), _parse_indices(fcfp4_field)


def _load_reference_rows() -> list[tuple[str, list[int], list[int]]]:
    path = _reference_csv_path()
    with path.open("r", encoding="utf-8") as f:
        # Skip header: smiles,ecfp4,fcfp4,rdkit
        next(f)
        rows = []
        for raw_line in f:
            line = raw_line.rstrip("\r\n")
            if line:
                rows.append(_parse_reference_line(line))
    return rows


def _nonzero_indices(fp_array: np.ndarray) -> list[int]:
    return np.flatnonzero(fp_array).astype(int).tolist()


def test_ecfp4_matches_reference_bits_from_csv():
    rows = _load_reference_rows()
    assert rows

    fp = descriptor.MorganFingerprint.ecfp4()
    for smiles, ecfp4_indices, _ in rows:
        result = fp(Molecule.from_smiles(smiles))
        assert _nonzero_indices(result) == ecfp4_indices


def test_fcfp4_matches_reference_bits_from_csv():
    rows = _load_reference_rows()
    assert rows

    fp = descriptor.MorganFingerprint.fcfp4()
    for smiles, _, fcfp4_indices in rows:
        result = fp(Molecule.from_smiles(smiles))
        assert _nonzero_indices(result) == fcfp4_indices
