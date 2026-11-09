# DMG Spatial Study — Analysis Scripts

Numbered in run order. Each script is self-documenting (see module docstring for
inputs/outputs). Raw data sources are listed in `../tables/table01_dataset_overview.csv`
(download URLs verified); see `../../protocol_dmg_spatial.md` for the full study protocol.

| # | Script | Purpose |
|---|--------|---------|
| 01 | `01_visium_loader.py` | Loaders for Ren / Kordowski / Huchede Visium and Xenium datasets |
| 02 | `02_collot_loader.py` | Loader for held-out Collot 2026 Visium sections (Zenodo 10.5281/zenodo.17922773) |
| 03 | `03_c2l_run.py` | cell2location reference training + per-section mapping (resumable; 24 discovery + 3 Collot sections) |
| 04 | `04_qc_domains.py` | Per-section Visium QC and spatial domains |
| 05 | `05_niches.py` | Niche calling (neighbor-smoothed composition, Leiden) + colocalization statistics |
| 06 | `06_lr_analysis.py` | LIANA ligand–receptor inference per section on dominant-population labels |
| 07 | `07_fix_spatial_z.py` | Spatial support (neighborhood enrichment z) for LR pairs; re-filter |
| 08 | `08_aggregate_validation.py` | Cross-cohort validation: colocalization concordance, rostro-caudal profiles |
| 09 | `09_aggregate_lr.py` | Cross-section/cross-cohort LR recurrence; overlap with Kordowski published LR set |
| 10 | `10_finalize_lr.py` | Final LR tables: 250 same-population headline set + 79 interface crosstalk pairs |
| 11 | `11_xenium_analysis.py` | Xenium marker-score state assignment + neighborhood enrichment |
| 12 | `12_figures_publication.py` | All 9 publication figures (PNG + editable SVG); run `python3 12_figures_publication.py all` or pass figure keys (e.g. `fig04 fig05`) |
| 13 | `13_figures_supplementary.py` | Supplementary figures S1 (DMG vs GBM abundance) and S2 (state verification) |

## Notes

- Scripts 03–11 write intermediates to a shared working volume (`c2l_out/`, `niches/`,
  `lr_out/`, `validation/`, `audit/`); script 12 reads those intermediates and writes
  final figures to `../figures/`.
- The reference atlas (Liu + Jessa harmonized scRNA, 13 labels) is built upstream of
  script 03; its construction is documented in the report Methods and the execution-trace
  notebooks (`../../execution_trace/`).
- Figure style (script 12): Liberation Sans 7 pt base, no panel titles (captions carry
  descriptions), lowercase bold panel tags (a, b, c), 300 dpi PNG + editable-text SVG.
