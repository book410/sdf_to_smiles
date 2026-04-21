# Project Memory

## Overview

- Project name: `sdf_to_smiles`
- Goal: convert SDF records to SMILES and validate the conversion against real
  PubChem compounds.

## Core Files

- `sdf_to_smiles.py`
  Convert local `.sdf` files to SMILES.
- `pubchem_smiles_roundtrip_test.py`
  Download real PubChem records, convert their SDF with the local converter,
  normalize both SMILES strings with RDKit, and compare the normalized result.
- `pubchem_test_data/summary.tsv`
  Round-trip validation summary for the current cached PubChem test set.

## Current Behavior

- `sdf_to_smiles.py` keeps explicit hydrogens by default.
- Pass `--implicit-h` to remove explicit hydrogens before SMILES export.
- SMILES output is canonical by default.
- Pass `--no-canonical` to disable canonicalization.
- Pass `--include-header` to emit `smiles<delimiter>name` as the first row.
- The default output delimiter is a tab.

## Validation Rules

- Compare PubChem and converted SMILES after RDKit normalization.
- Treat `normalized_match=True` as the primary correctness signal.
- Do not treat `raw_match=False` as a failure by itself.
- Use the PubChem test script with retry and timeout controls when the network
  is slow.

## PubChem Test Set

- The default PubChem validation set contains 20 real compounds.
- The test script caches downloaded SDF files and metadata under
  `pubchem_test_data/`.
- Recommended network settings:
  `--timeout 90 --retries 4 --retry-backoff 3`

## Verified Results

- The 20-compound PubChem validation set completed successfully.
- Final validation status: `20/20 normalized_match=True`
- The summary file is stored at `pubchem_test_data/summary.tsv`.

## Common Commands

```bash
python3 sdf_to_smiles.py input.sdf
python3 sdf_to_smiles.py input.sdf --implicit-h
python3 sdf_to_smiles.py input.sdf -o output.txt --include-header
python3 pubchem_smiles_roundtrip_test.py --timeout 90 --retries 4 --retry-backoff 3
```

## Repository Notes

- The project directory has been cleaned to keep only the main scripts, README,
  and cached PubChem validation data.
- `pubchem_smiles_roundtrip_test.py` no longer uses hard-coded absolute project
  paths; it derives defaults from the script directory.

## Next Steps

- Add grouped reporting by test level or compound class if needed.
- Add a `tests/` directory for unit and regression tests.
- Extend validation with charged, salt-form, or multi-component compounds if
  broader coverage is needed.
