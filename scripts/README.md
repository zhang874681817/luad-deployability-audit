# Scripts

The Python scripts in this directory passed syntax checks with `py_compile` during final packaging.

The `.mjs` scripts passed `node --check`. Several `.mjs` files depend on the Codex artifact-tool environment used to generate workbook artifacts; the generated workbooks are included in `../supplementary/`.

Most analysis scripts assume the repository root contains local working directories such as `data/`, `work/` and `outputs/`. Raw third-party expression matrices are not redistributed and must be retrieved from GEO, TCGA/GDC or original repositories when full regeneration is needed.
