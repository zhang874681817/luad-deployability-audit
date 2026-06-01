# Reproducibility

## Environment

Create a Python environment with either:

```powershell
conda env create -f environment.yml
conda activate luad-deployability-audit
```

or:

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
pip install -r requirements.txt
```

The final packaging pass checked Python script syntax with `py_compile` and JavaScript module syntax with `node --check`. Full end-to-end regeneration requires raw public GEO/TCGA inputs that are not redistributed in this repository.

## Data

Raw third-party expression matrices and clinical data are not included. Retrieve public data from GEO and TCGA/GDC using the accessions documented in `supplementary/` and `source_data/`.

Derived data needed to inspect the manuscript-level results are included:

- `source_data/`: figure and table source data.
- `figures/`: main figure PDFs/PNGs.
- `supplementary/`: supplementary workbooks and supplementary figures.

## Suggested Reproduction Path

1. Inspect `source_data/Table1_primary_candidate_characteristics.csv` and `source_data/Table2_external_performance_primary_candidates.csv` for the main tables.
2. Inspect `source_data/Fig*.csv` files for the main figure inputs.
3. Use `scripts/` to rerun selected analysis stages after placing raw/processed inputs in the expected local directories:
   - `data/raw/`
   - `data/processed/`
   - `work/`
   - `outputs/`
4. Compare regenerated outputs against the supplied `source_data/`, `figures/` and `supplementary/` files.

## Important Interpretation Constraints

- TCGA-LUAD GDC STAR TPM rescoring is a formula-direction sanity check, not a clean independent external validation.
- Median fallback grouping supports score portability only, not original cutoff validation.
- The 282 pre-audit full-text-target exclusions are included as accurate author-performed manual table-screening aggregate PRISMA-style categories.
- The audit supports further deployability evaluation of selected candidates, not immediate clinical deployment.

## Known Limitations

Some historical extraction artifacts and raw data files were not redistributed. This affects full raw-data regeneration but not inspection of the final derived source data supplied in this release.
