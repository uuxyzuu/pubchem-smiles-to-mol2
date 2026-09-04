# PubChem SMILES to MOL2

Batch-convert a PubChem CSV export—or another CSV containing SMILES—into separate, energy-minimized 3D MOL2 files.

PubChem commonly exports SMILES in CSV or SDF form but does not provide a convenient one-click batch export of separate MOL2 files. This tool fills that workflow gap while producing a transparent row-by-row conversion report.

## What it does

- Automatically detects common `SMILES`, compound-ID, and name columns.
- Preserves defined SMILES stereochemistry and formal charges.
- Adds explicit hydrogens.
- Generates multiple 3D conformers with RDKit ETKDGv3.
- Minimizes conformers with MMFF94s, falling back to UFF when needed.
- Keeps the lowest calculated-energy conformer.
- Writes one MOL2 file per compound using Open Babel.
- Creates `conversion_report.csv` and a ZIP archive of all successful structures.
- Continues processing if an individual molecule fails.

## Installation

Python 3.10 or newer is recommended.

```bash
git clone https://github.com/uuxyzuu/pubchem-smiles-to-mol2.git
cd pubchem-smiles-to-mol2
python -m venv .venv
source .venv/bin/activate  # Windows: .venv\\Scripts\\activate
python -m pip install -r requirements.txt
```

## Usage

For a standard PubChem CSV export containing `cid`, `Name`, and `SMILES`:

```bash
python smiles_to_mol2.py PubChem_compounds.csv -o results
```

Output:

```text
results/
├── conversion_report.csv
├── PubChem_compounds_mol2.zip
└── mol2_structures/
    ├── CID_123_Compound_A.mol2
    └── CID_456_Compound_B.mol2
```

For a CSV with custom column names:

```bash
python smiles_to_mol2.py compounds.csv -o results \
  --smiles-column canonical_smiles \
  --id-column compound_id \
  --name-column compound_name
```

Adjust conformer sampling or minimization:

```bash
python smiles_to_mol2.py compounds.csv -o results \
  --num-confs 20 --max-iters 2000 --seed 42
```

The program returns exit code `0` when every row succeeds, `2` when conversion is partially successful, and `1` for an input/configuration error.

## Chemical assumptions

- Protonation, tautomeric form, salts, and disconnected fragments are retained as encoded in the input SMILES. The tool does **not** predict protonation at a target pH.
- Undefined stereocenters remain undefined; they are not enumerated.
- MMFF94s/UFF minimization gives a practical starting geometry, not a quantum-mechanical global minimum.
- MOL2 atom types and charges should be checked against the requirements of the downstream docking or molecular-dynamics program.
- Metal-containing compounds and unusual covalent species may require manual parameterization.

## Input example

```csv
cid,Name,SMILES
2244,Aspirin,CC(=O)OC1=CC=CC=C1C(=O)O
2519,Caffeine,CN1C=NC2=C1C(=O)N(C(=O)N2C)C
```

## Testing

```bash
python -m pip install pytest
pytest -q
```

## License

MIT License. See [LICENSE](LICENSE).
