#!/usr/bin/env python3
"""Batch-convert a CSV containing SMILES into minimized 3D MOL2 files."""

from __future__ import annotations

import argparse
import csv
import re
import sys
import zipfile
from collections import Counter
from pathlib import Path

from openbabel import pybel
from rdkit import Chem
from rdkit.Chem import AllChem


def safe_filename(value: str, fallback: str = "compound") -> str:
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", str(value).strip()).strip("._")
    return cleaned[:120] or fallback


def find_column(fieldnames: list[str], requested: str | None, candidates: tuple[str, ...]) -> str:
    if requested:
        if requested not in fieldnames:
            raise ValueError(f"Column {requested!r} was not found. Available: {', '.join(fieldnames)}")
        return requested
    lower_map = {name.casefold(): name for name in fieldnames}
    for candidate in candidates:
        if candidate.casefold() in lower_map:
            return lower_map[candidate.casefold()]
    raise ValueError(f"Could not identify a required column. Available: {', '.join(fieldnames)}")


def build_3d(smiles: str, title: str, num_confs: int, max_iters: int, seed: int):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("SMILES parsing failed")
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    mol = Chem.AddHs(mol)

    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    params.pruneRmsThresh = 0.5
    params.useSmallRingTorsions = True
    requested = min(num_confs, 3) if mol.GetNumHeavyAtoms() > 80 else num_confs
    conf_ids = list(AllChem.EmbedMultipleConfs(mol, numConfs=requested, params=params))
    if not conf_ids:
        params.useRandomCoords = True
        conf_id = AllChem.EmbedMolecule(mol, params)
        if conf_id < 0:
            raise ValueError("3D embedding failed")
        conf_ids = [conf_id]

    energies: list[tuple[float, int]] = []
    method = "MMFF94s"
    if AllChem.MMFFHasAllMoleculeParams(mol):
        props = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant="MMFF94s")
        for conf_id in conf_ids:
            try:
                ff = AllChem.MMFFGetMoleculeForceField(mol, props, confId=conf_id)
                ff.Minimize(maxIts=max_iters)
                energies.append((ff.CalcEnergy(), conf_id))
            except Exception:
                pass

    if not energies:
        method = "UFF"
        for conf_id in conf_ids:
            try:
                ff = AllChem.UFFGetMoleculeForceField(mol, confId=conf_id)
                ff.Minimize(maxIts=max_iters)
                energies.append((ff.CalcEnergy(), conf_id))
            except Exception:
                pass

    if energies:
        energy, best_conf = min(energies)
    else:
        method, energy, best_conf = "unoptimized", None, conf_ids[0]
    mol.SetProp("_Name", title)
    return mol, best_conf, method, energy


def write_mol2(mol, conf_id: int, title: str, destination: Path) -> None:
    sdf_block = Chem.MolToMolBlock(mol, confId=conf_id)
    converted = pybel.readstring("sdf", sdf_block)
    converted.title = title
    converted.write("mol2", str(destination), overwrite=True)
    if not destination.exists() or destination.stat().st_size == 0:
        raise ValueError("MOL2 export failed")


def convert(args: argparse.Namespace) -> int:
    input_path = Path(args.input_csv)
    output_dir = Path(args.output_dir)
    structures_dir = output_dir / "mol2_structures"
    structures_dir.mkdir(parents=True, exist_ok=True)

    with input_path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            raise ValueError("The CSV has no header")
        rows = list(reader)
        smiles_col = find_column(reader.fieldnames, args.smiles_column, ("SMILES", "CanonicalSMILES", "IsomericSMILES"))
        id_col = find_column(reader.fieldnames, args.id_column, ("cid", "CID", "id", "ID"))
        name_col = None
        try:
            name_col = find_column(reader.fieldnames, args.name_column, ("Name", "name", "Title", "title"))
        except ValueError:
            if args.name_column:
                raise

    report = []
    used_names: Counter[str] = Counter()
    for index, row in enumerate(rows, start=1):
        compound_id = (row.get(id_col) or "").strip() or f"row_{index}"
        compound_name = (row.get(name_col) or "").strip() if name_col else compound_id
        compound_name = compound_name or compound_id
        smiles = (row.get(smiles_col) or "").strip()
        stem = f"{safe_filename(id_col.upper())}_{safe_filename(compound_id)}_{safe_filename(compound_name)}"
        used_names[stem] += 1
        if used_names[stem] > 1:
            stem = f"{stem}_{used_names[stem]}"
        destination = structures_dir / f"{stem}.mol2"
        status, method, energy, message = "failed", "", None, ""
        try:
            if not smiles:
                raise ValueError("Missing SMILES")
            title = f"{id_col} {compound_id} | {compound_name}"
            mol, conf_id, method, energy = build_3d(
                smiles, title, args.num_confs, args.max_iters, args.seed
            )
            write_mol2(mol, conf_id, title, destination)
            status = "success"
        except Exception as exc:
            message = str(exc)
        report.append({
            "row": index,
            "compound_id": compound_id,
            "name": compound_name,
            "smiles": smiles,
            "status": status,
            "optimization": method,
            "energy_kcal_mol": "" if energy is None else f"{energy:.6f}",
            "output_file": destination.name if status == "success" else "",
            "message": message,
        })
        if index % 25 == 0 or index == len(rows):
            print(f"Processed {index}/{len(rows)}", flush=True)

    report_path = output_dir / "conversion_report.csv"
    with report_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=report[0].keys() if report else ["row", "status"])
        writer.writeheader()
        writer.writerows(report)

    archive_path = output_dir / f"{safe_filename(input_path.stem)}_mol2.zip"
    with zipfile.ZipFile(archive_path, "w", zipfile.ZIP_DEFLATED) as archive:
        archive.write(report_path, report_path.name)
        for path in sorted(structures_dir.glob("*.mol2")):
            archive.write(path, f"mol2_structures/{path.name}")

    success = sum(item["status"] == "success" for item in report)
    methods = Counter(item["optimization"] for item in report if item["status"] == "success")
    print(f"Completed: {success}/{len(report)} successful; methods={dict(methods)}")
    print(f"Archive: {archive_path}")
    return 0 if success == len(report) else 2


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_csv", help="CSV containing compound identifiers and SMILES")
    parser.add_argument("-o", "--output-dir", default="mol2_output")
    parser.add_argument("--smiles-column")
    parser.add_argument("--id-column")
    parser.add_argument("--name-column")
    parser.add_argument("--num-confs", type=int, default=10)
    parser.add_argument("--max-iters", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20260903)
    args = parser.parse_args(argv)
    if args.num_confs < 1 or args.max_iters < 1:
        parser.error("--num-confs and --max-iters must be positive")
    return args


def main(argv=None) -> int:
    try:
        return convert(parse_args(argv))
    except (OSError, ValueError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
