import csv
from argparse import Namespace

from openbabel import pybel

from smiles_to_mol2 import convert, safe_filename


def test_safe_filename():
    assert safe_filename("CID 1 / aspirin") == "CID_1_aspirin"


def test_end_to_end(tmp_path):
    source = tmp_path / "input.csv"
    source.write_text("cid,Name,SMILES\n2244,Aspirin,CC(=O)OC1=CC=CC=C1C(=O)O\n")
    output = tmp_path / "results"
    args = Namespace(
        input_csv=str(source), output_dir=str(output), smiles_column=None,
        id_column=None, name_column=None, num_confs=2, max_iters=100, seed=42,
    )
    assert convert(args) == 0
    mol2_files = list((output / "mol2_structures").glob("*.mol2"))
    assert len(mol2_files) == 1
    assert len(list(pybel.readfile("mol2", str(mol2_files[0])))) == 1
    with (output / "conversion_report.csv").open() as handle:
        report = list(csv.DictReader(handle))
    assert report[0]["status"] == "success"
