# Public Release Notes

This upload-ready directory was prepared from the manuscript repository package.

Final cleanup actions:

- Removed author-completion placeholder files from the upload-ready package.
- Removed Python bytecode caches.
- Removed internal publication-readiness CSV files that tracked missing repository URL/DOI status.
- Replaced provisional environment pins with a reproducible minimum-version environment.
- Added `.gitignore`, GitHub upload instructions and a PowerShell upload helper.
- Preserved scripts, main figures, figure/table source data and supplementary data.

Validation performed:

- Python scripts: syntax checked with `py_compile`.
- JavaScript modules: syntax checked with `node --check`.
- No files over the GitHub 100 MB single-file limit were found.

Items to complete after upload:

- Record the GitHub URL.
- Record the frozen commit hash printed by `upload_to_github.ps1`.
- Archive the GitHub release through Zenodo or OSF and record the DOI.
