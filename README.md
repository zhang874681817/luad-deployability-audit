# LUAD transcriptomic prognostic signature deployability audit

This repository contains code, figure source data, derived tables and supplementary materials for a systematic computational deployability audit of published lung adenocarcinoma (LUAD) transcriptomic prognostic signatures.

## Scope

This is not a new LUAD prognostic-signature discovery study. The audit evaluates published signatures for formula reproducibility, external scoreability, score portability, reported-cutoff transportability and clinical deployability readiness.

## Repository and Archive

- Public repository: https://github.com/zhang874681817/luad-deployability-audit
- Current release: `v1.0.1`
- Concept DOI (resolves to the latest Zenodo version): https://doi.org/10.5281/zenodo.20507011
- Previous `v0.1.0` version DOI: https://doi.org/10.5281/zenodo.20507012

## What Is Included

- `scripts/`: analysis, audit, scoring, search-merge and figure-generation scripts used during the project.
- `source_data/`: source data for the main figures and tables, including the final PubMed/Embase/Web of Science merge outputs and the authoritative Appendix 3 snapshot.
- `figures/`: publication figure files supplied with the submission package.
- `supplementary/`: final appendices, supplementary workbooks and supplementary figures.
- `REPRODUCIBILITY.md`: practical notes for environment setup, data retrieval and rerunning analyses.
- `FILE_MANIFEST_SHA256.csv`: file inventory with SHA256 hashes for the upload-ready package.

## What Is Not Included

Raw third-party expression matrices are not redistributed. Users should retrieve GEO and TCGA/GDC data from the original public repositories under their applicable terms.

Some scripts expect raw or processed input directories such as `data/`, `work` and `outputs` that were used during the local analysis build. The final derived source data and supplementary files needed to inspect the reported manuscript results are included here.

## Main Reproducibility Entry Points

Use `REPRODUCIBILITY.md` first. In brief:

1. Create a Python environment from `environment.yml` or install packages from `requirements.txt`.
2. Retrieve raw public data from GEO and TCGA/GDC as needed.
3. Use the scripts in `scripts/` to regenerate processed tables and figures.
4. Compare regenerated outputs with `source_data/`, `figures` and `supplementary`.

The release also includes the scripts and derived outputs for the deployability scorecard, sensitivity analyses, exploratory neural-network benchmark and post hoc stability-guided feature analyses. These analyses do not constitute independent validation of a new clinical model.

## Citation

Please cite the associated manuscript and archived repository:

Deployability of published lung adenocarcinoma transcriptomic prognostic signatures: a systematic audit of computational reproducibility and external transportability.

Archived concept DOI: https://doi.org/10.5281/zenodo.20507011

Release tag: `v1.0.1`

## License

Code is released under the MIT License. Derived documentation, figures and tables are released under CC BY 4.0. Raw third-party datasets remain governed by GEO, TCGA/GDC and original-source terms and are not redistributed here.
