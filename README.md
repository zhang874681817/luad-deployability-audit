# LUAD transcriptomic prognostic signature deployability audit

This repository contains code, figure source data, derived tables and supplementary materials for a systematic computational deployability audit of published lung adenocarcinoma (LUAD) transcriptomic prognostic signatures.

## Scope

This is not a new LUAD prognostic-signature discovery study. The audit evaluates published signatures for formula reproducibility, external scoreability, score portability, reported-cutoff transportability and clinical deployability readiness.

## What Is Included

- `scripts/`: analysis, audit, scoring and figure-generation scripts used during the project.
- `source_data/`: source data for the main figures and tables.
- `figures/`: publication figure files supplied with the submission package.
- `supplementary/`: supplementary workbooks and supplementary figures.
- `REPRODUCIBILITY.md`: practical notes for environment setup, data retrieval and rerunning analyses.
- `FILE_MANIFEST_SHA256.csv`: file inventory with SHA256 hashes for the upload-ready package.

## What Is Not Included

Raw third-party expression matrices are not redistributed. Users should retrieve GEO and TCGA/GDC data from the original public repositories under their applicable terms.

Some scripts expect raw or processed input directories such as `data/`, `work/` and `outputs/` that were used during the local analysis build. The final derived source data and supplementary files needed to inspect the reported manuscript results are included here.

## Main Reproducibility Entry Points

Use `REPRODUCIBILITY.md` first. In brief:

1. Create a Python environment from `environment.yml` or install packages from `requirements.txt`.
2. Retrieve raw public data from GEO and TCGA/GDC as needed.
3. Use the scripts in `scripts/` to regenerate processed tables and figures.
4. Compare regenerated outputs with `source_data/`, `figures/` and `supplementary/`.

Several `.mjs` workbook-builder scripts were created in the Codex artifact-tool environment and are retained for provenance. The corresponding workbook outputs are already included in `supplementary/`.

## Citation

Please cite the associated manuscript once available:

Deployability audit of published lung adenocarcinoma transcriptomic prognostic signatures.

## License

Code is released under the MIT License. Derived documentation, figures and tables are released under CC BY 4.0. Raw third-party datasets remain governed by GEO, TCGA/GDC and original-source terms and are not redistributed here.
