# Reproducibility Runbook

This runbook is a compact public-release guide for the LUAD deployability audit package.

## Before Running Scripts

1. Install Python dependencies from `requirements.txt` or `environment.yml`.
2. Download required raw public datasets from GEO and TCGA/GDC.
3. Place local raw and processed inputs under the expected directories:
   - `data/raw/`
   - `data/processed/`
   - `work/`
   - `outputs/`

These local working directories are intentionally ignored by `.gitignore`.

## Script Families

- `extract_*` and `prepare_*`: extraction and preprocessing helpers.
- `run_*`: validation, scoring and sensitivity analyses.
- `build_*`: table, figure, workbook and package builders.
- `audit_*`, `review_*`, `compare_*`, `qc_*`: audit and quality-control helpers.

## Public Release Checks

Before uploading or after cloning, run:

```powershell
python -m py_compile (Get-ChildItem scripts -Filter *.py).FullName
foreach ($f in Get-ChildItem scripts -Filter *.mjs) { node --check $f.FullName }
```

The upload-ready package was checked this way during final cleanup.

## Outputs To Compare

- Main figures: `figures/`
- Figure/table source data: `source_data/`
- Supplementary workbooks and figures: `supplementary/`

Raw third-party matrices are not redistributed.

## Public Archive

- Repository: https://github.com/zhang874681817/luad-deployability-audit
- Release tag: `v1.0.0`
- Zenodo concept DOI: https://doi.org/10.5281/zenodo.20507011
- Previous `v0.1.0` version DOI: https://doi.org/10.5281/zenodo.20507012
