#!/usr/bin/env python3
"""Convert molecules in an SDF file to SMILES strings."""

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator, TextIO

try:
    from rdkit import Chem
except ImportError:  # pragma: no cover - depends on local environment
    Chem = None


@dataclass(frozen=True)
class SmilesRecord:
    index: int
    name: str
    smiles: str


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Read an SDF file and emit SMILES strings."
    )
    parser.add_argument("input", type=Path, help="Path to the input .sdf file")
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        help="Optional output file. Defaults to stdout.",
    )
    parser.add_argument(
        "--delimiter",
        default="\t",
        help="Delimiter between SMILES and molecule name. Default: tab.",
    )
    parser.add_argument(
        "--implicit-h",
        action="store_true",
        help="Remove explicit hydrogens before export.",
    )
    parser.add_argument(
        "--no-canonical",
        action="store_true",
        help="Do not canonicalize SMILES output.",
    )
    parser.add_argument(
        "--include-header",
        action="store_true",
        help="Write a header line: smiles<delimiter>name",
    )
    return parser


def ensure_rdkit() -> None:
    if Chem is not None:
        return

    raise SystemExit(
        "RDKit is required for SDF to SMILES conversion. "
        "Install it first, for example: conda install -c conda-forge rdkit"
    )


def iter_smiles_records(
    sdf_path: Path,
    *,
    implicit_hydrogens: bool,
    canonical: bool,
) -> Iterator[SmilesRecord]:
    ensure_rdkit()

    supplier = Chem.SDMolSupplier(str(sdf_path), removeHs=implicit_hydrogens)
    for index, mol in enumerate(supplier, start=1):
        if mol is None:
            print(
                f"warning: skipped invalid molecule at record {index}",
                file=sys.stderr,
            )
            continue

        name = mol.GetProp("_Name").strip() if mol.HasProp("_Name") else ""
        smiles = Chem.MolToSmiles(mol, canonical=canonical)
        yield SmilesRecord(index=index, name=name, smiles=smiles)


def format_records(
    records: Iterable[SmilesRecord],
    *,
    delimiter: str,
    include_header: bool,
) -> Iterator[str]:
    if include_header:
        yield f"smiles{delimiter}name"

    for record in records:
        yield f"{record.smiles}{delimiter}{record.name}".rstrip()


def write_lines(lines: Iterable[str], destination: TextIO) -> None:
    for line in lines:
        destination.write(line)
        destination.write("\n")


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()

    if not args.input.is_file():
        parser.error(f"input file does not exist: {args.input}")

    records = iter_smiles_records(
        args.input,
        implicit_hydrogens=args.implicit_h,
        canonical=not args.no_canonical,
    )
    lines = format_records(
        records,
        delimiter=args.delimiter,
        include_header=args.include_header,
    )

    if args.output is None:
        write_lines(lines, sys.stdout)
        return 0

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="\n") as handle:
        write_lines(lines, handle)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
