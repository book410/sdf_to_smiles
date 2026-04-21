# smiles

This directory contains one conversion script and one validation script.

## Files

- `sdf_to_smiles.py`: convert SDF records to SMILES.
- `pubchem_smiles_roundtrip_test.py`: download real PubChem records and compare
  converted SMILES against PubChem SMILES after RDKit normalization.
- `pubchem_test_data/`: cached PubChem SDF files, metadata, and test summaries.

## Usage

Convert a local SDF file:

```bash
python3 sdf_to_smiles.py input.sdf
python3 sdf_to_smiles.py input.sdf --implicit-h
python3 sdf_to_smiles.py input.sdf -o output.txt --include-header
```

Run the default 20-compound PubChem validation set:

```bash
python3 pubchem_smiles_roundtrip_test.py --timeout 90 --retries 4 --retry-backoff 3
```
