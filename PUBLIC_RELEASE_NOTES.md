# Public Release Notes

## v1.0.1 (2026-10-12)

This visual-only correction removes the former non-data Figure 5D panel and
reflows Figure 5 as three data-driven panels. The clinical increment results,
likelihood-ratio test proportions, deployability scorecard values, source data,
study conclusions and interpretation boundaries are unchanged.

## v1.0.0 (2026-10-12)

This release synchronizes the public repository with the final manuscript audit package.

### Added or updated

- Final three-database reconciliation: 10,061 records before cross-database deduplication, 2,799 cross-database duplicates removed and 7,262 unique verification records.
- Reconciled 178-article audit corpus and final PRISMA source data.
- Revised Figures 1-6 and their source data.
- Ten-domain deployability scorecard and sensitivity analyses.
- Exploratory neural-network benchmark, optimization summaries and stability-guided 5-gene/40-gene analyses.
- Final Appendix 3 data snapshot, clinical-model table and figure-source package.
- Updated author, citation, licensing and reproducibility metadata.

### Interpretation boundary

This is a systematic computational deployability, reproducibility and external transportability audit of published LUAD transcriptomic prognostic signatures. It is not a new prognostic-model discovery study and does not constitute clinical deployment validation. The AI and core-score analyses are exploratory and post hoc.

### Validation

- Python scripts were syntax checked.
- The file manifest was regenerated with SHA256 checksums.
- No release file exceeds the GitHub 100 MB single-file limit.
- Main-figure files were checked against the final manuscript package.
