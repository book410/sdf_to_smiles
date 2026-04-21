#!/usr/bin/env python3
"""Validate sdf_to_smiles.py against real PubChem compound records."""

from __future__ import annotations

import argparse
import csv
import json
import re
import time
import subprocess
import sys
import urllib.error
import urllib.request
from dataclasses import dataclass
from pathlib import Path

from rdkit import Chem

PROJECT_ROOT = Path(__file__).resolve().parent

PUBCHEM_SDF_URL = (
    "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/SDF?record_type=2d"
)
PUBCHEM_PROPERTY_URL = (
    "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/"
    "{cid}/property/Title,SMILES/CSV"
)
DEFAULT_CIDS = [
    702,    # Ethanol
    176,    # Acetic Acid
    241,    # Benzene
    1140,   # Toluene
    8078,   # Cyclohexane
    996,    # Phenol
    240,    # Benzaldehyde
    243,    # Benzoic Acid
    750,    # Glycine
    5950,   # L-Alanine
    612,    # Lactic Acid
    5793,   # D-Glucose
    2244,   # Aspirin
    1983,   # Acetaminophen
    2519,   # Caffeine
    3672,   # Ibuprofen
    5329,   # Sulfamethoxazole
    5997,   # Cholesterol
    33613,  # Amoxicillin
    60823,  # Atorvastatin
]
DEFAULT_TIMEOUT = 90
DEFAULT_RETRIES = 4
DEFAULT_RETRY_BACKOFF = 3.0


@dataclass(frozen=True)
class TestCaseResult:
    cid: int
    title: str
    sdf_path: Path
    pubchem_smiles: str
    converted_smiles: str
    raw_match: bool
    normalized_pubchem: str
    normalized_converted: str
    normalized_match: bool


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Download real compound SDF files from PubChem, convert them with "
            "sdf_to_smiles.py, and compare against PubChem SMILES. "
            "By default this runs a 20-compound simple-to-complex test set."
        )
    )
    parser.add_argument(
        "--cid",
        action="append",
        type=int,
        dest="cids",
        help="PubChem CID to test. Repeat to add more than one CID.",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=PROJECT_ROOT / "pubchem_test_data",
        help="Directory for downloaded SDF files and the summary TSV.",
    )
    parser.add_argument(
        "--python",
        default=sys.executable,
        help="Python executable used to run sdf_to_smiles.py.",
    )
    parser.add_argument(
        "--converter",
        type=Path,
        default=PROJECT_ROOT / "sdf_to_smiles.py",
        help="Path to the sdf_to_smiles.py script under test.",
    )
    parser.add_argument(
        "--timeout",
        type=int,
        default=DEFAULT_TIMEOUT,
        help="Per-request timeout in seconds. Default: %(default)s.",
    )
    parser.add_argument(
        "--retries",
        type=int,
        default=DEFAULT_RETRIES,
        help="Number of attempts for each PubChem request. Default: %(default)s.",
    )
    parser.add_argument(
        "--retry-backoff",
        type=float,
        default=DEFAULT_RETRY_BACKOFF,
        help=(
            "Base backoff in seconds between retries. Sleep time grows linearly "
            "with the attempt number. Default: %(default)s."
        ),
    )
    return parser


def fetch_text(
    url: str,
    *,
    timeout: int,
    retries: int,
    retry_backoff: float,
) -> str:
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "smiles-test/1.0"},
    )
    last_error: OSError | None = None
    for attempt in range(1, retries + 1):
        try:
            with urllib.request.urlopen(request, timeout=timeout) as response:
                return response.read().decode("utf-8")
        except (TimeoutError, urllib.error.URLError, OSError) as exc:
            last_error = exc
            if attempt == retries:
                break
            time.sleep(retry_backoff * attempt)
    assert last_error is not None
    raise last_error


def fetch_bytes(
    url: str,
    *,
    timeout: int,
    retries: int,
    retry_backoff: float,
) -> bytes:
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "smiles-test/1.0"},
    )
    last_error: OSError | None = None
    for attempt in range(1, retries + 1):
        try:
            with urllib.request.urlopen(request, timeout=timeout) as response:
                return response.read()
        except (TimeoutError, urllib.error.URLError, OSError) as exc:
            last_error = exc
            if attempt == retries:
                break
            time.sleep(retry_backoff * attempt)
    assert last_error is not None
    raise last_error


def slugify(text: str) -> str:
    slug = re.sub(r"[^A-Za-z0-9]+", "_", text.strip()).strip("_").lower()
    return slug or "compound"


def normalize_smiles(smiles: str) -> str:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"RDKit could not parse SMILES: {smiles}")
    mol = Chem.RemoveHs(mol)
    return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)


def metadata_cache_path(cid: int, out_dir: Path) -> Path:
    return out_dir / f"{cid}_metadata.json"


def fetch_pubchem_metadata(
    cid: int,
    *,
    out_dir: Path,
    timeout: int,
    retries: int,
    retry_backoff: float,
) -> tuple[str, str]:
    cache_path = metadata_cache_path(cid, out_dir)
    if cache_path.is_file():
        payload = json.loads(cache_path.read_text(encoding="utf-8"))
        return payload["Title"], payload["SMILES"]

    csv_text = fetch_text(
        PUBCHEM_PROPERTY_URL.format(cid=cid),
        timeout=timeout,
        retries=retries,
        retry_backoff=retry_backoff,
    )
    rows = list(csv.DictReader(csv_text.splitlines()))
    if len(rows) != 1:
        raise ValueError(f"expected one property row for CID {cid}, got {len(rows)}")
    row = rows[0]
    cache_path.write_text(
        json.dumps({"Title": row["Title"], "SMILES": row["SMILES"]}, ensure_ascii=True),
        encoding="utf-8",
    )
    return row["Title"], row["SMILES"]


def download_sdf(
    cid: int,
    title: str,
    out_dir: Path,
    *,
    timeout: int,
    retries: int,
    retry_backoff: float,
) -> Path:
    path = out_dir / f"{cid}_{slugify(title)}.sdf"
    if path.is_file() and path.stat().st_size > 0:
        return path
    sdf_bytes = fetch_bytes(
        PUBCHEM_SDF_URL.format(cid=cid),
        timeout=timeout,
        retries=retries,
        retry_backoff=retry_backoff,
    )
    path.write_bytes(sdf_bytes)
    return path


def run_converter(python_bin: str, converter: Path, sdf_path: Path) -> str:
    completed = subprocess.run(
        [
            python_bin,
            str(converter),
            str(sdf_path),
            "--implicit-h",
        ],
        text=True,
        capture_output=True,
        check=False,
    )
    if completed.returncode != 0:
        raise RuntimeError(
            f"converter failed for {sdf_path} with stderr:\n{completed.stderr}"
        )

    lines = [line for line in completed.stdout.splitlines() if line.strip()]
    if len(lines) != 1:
        raise ValueError(
            f"expected one converted record for {sdf_path}, got {len(lines)}"
        )

    converted_smiles, *_rest = lines[0].split("\t", 1)
    return converted_smiles


def test_one_case(
    cid: int,
    *,
    python_bin: str,
    converter: Path,
    out_dir: Path,
    timeout: int,
    retries: int,
    retry_backoff: float,
) -> TestCaseResult:
    title, pubchem_smiles = fetch_pubchem_metadata(
        cid,
        out_dir=out_dir,
        timeout=timeout,
        retries=retries,
        retry_backoff=retry_backoff,
    )
    sdf_path = download_sdf(
        cid,
        title,
        out_dir,
        timeout=timeout,
        retries=retries,
        retry_backoff=retry_backoff,
    )
    converted_smiles = run_converter(python_bin, converter, sdf_path)
    normalized_pubchem = normalize_smiles(pubchem_smiles)
    normalized_converted = normalize_smiles(converted_smiles)
    return TestCaseResult(
        cid=cid,
        title=title,
        sdf_path=sdf_path,
        pubchem_smiles=pubchem_smiles,
        converted_smiles=converted_smiles,
        raw_match=(pubchem_smiles == converted_smiles),
        normalized_pubchem=normalized_pubchem,
        normalized_converted=normalized_converted,
        normalized_match=(normalized_pubchem == normalized_converted),
    )


def write_summary(results: list[TestCaseResult], out_dir: Path) -> Path:
    summary_path = out_dir / "summary.tsv"
    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "cid",
                "title",
                "sdf_path",
                "pubchem_smiles",
                "converted_smiles",
                "raw_match",
                "normalized_pubchem",
                "normalized_converted",
                "normalized_match",
            ]
        )
        for result in results:
            writer.writerow(
                [
                    result.cid,
                    result.title,
                    str(result.sdf_path),
                    result.pubchem_smiles,
                    result.converted_smiles,
                    result.raw_match,
                    result.normalized_pubchem,
                    result.normalized_converted,
                    result.normalized_match,
                ]
            )
    return summary_path


def main() -> int:
    args = build_parser().parse_args()
    cids = args.cids or DEFAULT_CIDS
    args.out_dir.mkdir(parents=True, exist_ok=True)

    results: list[TestCaseResult] = []
    failures = 0

    for cid in cids:
        try:
            result = test_one_case(
                cid,
                python_bin=args.python,
                converter=args.converter,
                out_dir=args.out_dir,
                timeout=args.timeout,
                retries=args.retries,
                retry_backoff=args.retry_backoff,
            )
        except (OSError, RuntimeError, ValueError, urllib.error.URLError) as exc:
            failures += 1
            print(f"[FAIL] CID {cid}: {exc}", file=sys.stderr)
            continue

        results.append(result)
        status = "PASS" if result.normalized_match else "FAIL"
        print(
            f"[{status}] CID {result.cid} {result.title} | "
            f"raw_match={result.raw_match} normalized_match={result.normalized_match}"
        )

    summary_path = write_summary(results, args.out_dir)
    print(f"summary_tsv={summary_path}")

    if failures:
        print(f"failed_cases={failures}", file=sys.stderr)

    return 0 if failures == 0 and all(r.normalized_match for r in results) else 1


if __name__ == "__main__":
    raise SystemExit(main())
