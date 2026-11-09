# External validation dataset candidates — DMG/DIPG (H3K27-altered) spatial & single-cell audit

**Search date:** 2026-09-26
**Already in use (excluded, overlap flagged only):** GSE194329 (Ren, Visium), GSE280990 (Kordowski, Visium), GSE268577 (Huchedé, Visium), GSE271936 (Xenium pHGG), GSE184357 (Liu 2022, Smart-seq2 sc/snRNA + scATAC), GSE210567 (Jessa 2022, 10x scRNA; superseries = GSE210568).

**Verification method:** every GEO accession below was retrieved and inspected directly on GEO (SOFT record via `acc.cgi` / E-utilities `esummary`), not inferred from papers.

---

## 1. New spatial transcriptomics datasets (human DMG)

| Accession / DOI | Study | Platform | Patients / samples | H3K27 status | Species | Raw data public? | Overlap with the six in-use datasets |
|---|---|---|---|---|---|---|---|
| **Zenodo: 10.5281/zenodo.17922773** (no GEO) | Collot & Ruiz-Moreno et al., *Cancer Cell* 2026 (IGSF11–VISTA; Princess Máxima Ctr, Rios lab). PMID 41576930 | **Visium** (16 samples, 51k spots) + **CosMx** (10 DMG cases / 23 patient samples, 26.5k cells) + snRNA-seq; also mouse IUE Visium + scRNA | 44 DMG patients total across modalities (new + re-analyzed public) | H3K27-altered DMG (autopsy program cohort) | Human + mouse model | **Yes** — open Zenodo, 2.7 GB: QC'd sc/snRNA matrices + Space Ranger outputs (human & mouse). CosMx included per paper's deposition table | **Partial** — they *re-analyzed* GSE194329, GSE280990, GSE210568 (=Jessa superseries), E-MTAB-15147. The 16 Visium + 23 CosMx samples are **newly generated** and not in our list |
| **GSE242352** | Webb et al. (USC, Craig lab) 2025, PMID 40514380 — spatial subclonal LOH/ecDNA in glioma | 10x Visium (fresh) + Visium CytAssist FFPE | 17 gliomas total; **exactly 1 DMG** ("Diffuse Midline Glioma A1", GSM7759704) | DMG sample present; H3 status not annotated on GEO | Human | Processed only (spaceranger outputs + integrated .rds); authors state raw withheld for privacy | None |
| **Damodharan et al. 2024**, *Sci Rep*, DOI 10.1038/s41598-024-73199-w | Dey lab (UW–Madison) | NanoString **GeoMx DSP** (transcriptome + proteome, ROI-based) | 8 DMG patients (4 pediatric, 4 adult; pons n=3, thalamus n=5), 95 ROIs | All H3 K27-altered (H3K27M IHC-guided ROIs) | Human | **No accession found.** Paper says "will be deposited into GEO… upon publication" (Sept 2024); targeted GEO searches (GeoMx/DSP + glioma; author; title) return nothing as of 2026-09-26. Effectively "on request" | None |

## 2. New spatial datasets (mouse DMG models)

| Accession | Study | Platform | Samples | Model | Raw public? | Overlap |
|---|---|---|---|---|---|---|
| **GSE299434** | Holliday et al., Children's Cancer Institute Australia 2025, PMID 41533777 (FACT+BET inhibitor therapy) | **Xenium Prime** (5,000-gene panel) | n=4 (treatment arms: CBL0137/JQ1 combo in immunocompetent setting) | IUE 24B7 H3.3K27M brainstem tumors, C57BL/6 | Yes (RAW.tar + integrated Seurat .rds) | None. Sibling series GSE299509/GSE272131 (human DMG cell-line RNA-seq) & GSE272129/GSE272130 (ChIP) are not spatial |
| **GSE246584** | Reitman lab (Duke) 2023, PMID 37904990 (ATM disruption + radiation) | **Xenium** (GPL33896, mouse) | n=4 (Atm intact/lost ± focal irradiation) | RCAS/tv-a H3.3K27M + p53 loss primary murine DMG | Yes (RAW.tar) | None |

## 3. Large pHGG cohorts with spatial component that include H3K27M

| Accession | Study | Spatial modality | Cohort | H3K27M content | Raw public? | Overlap |
|---|---|---|---|---|---|---|
| **GSE292623** | Sussman, De Raedt, Tan (CHOP/UPenn), *Cell Rep Med* 2026, DOI 10.1016/j.xcrm.2026.102766, PMID 42030938 | **CODEX spatial proteomics** (not transcriptomics); plus snRNA-seq + snATAC + WGS | 16 pHGG patients, longitudinal (dx + post-therapy) | Includes H3K27M-altered pHGG (molecularly diverse cohort) | GSE292623 hosts only the 18 cell-line irradiation RNA-seq samples; methylation on Zenodo 10.5281/zenodo.15092768; patient snRNA/snATAC availability to be confirmed from paper (likely mixed open/controlled) | None |
| GSE271936 | — | Xenium | n=4 pHGG | — | — | **Already in use** |

## 4. New sc/snRNA cohorts (reference augmentation candidates)

| Accession / DOI | Study | Platform | DMG-relevant content | Species | Raw public? | Overlap |
|---|---|---|---|---|---|---|
| **Okada et al. 2026**, *Neuro-Oncology* noag119, PMID 42178382 | Diaz/Nazarian/Wang labs | snRNA-seq + snATAC-seq | **N=22 postmortem DMG, H3K27-altered** (imipridone-treated and -naive); +IF myeloid n=46 | Human | **No GEO deposition found as of 2026-09-26** (paper published 2026-05-22; likely pending). Strong candidate once released | None |
| **GSE231859** (superseries GSE231860, + bulk GSE231858) | DeSisto/Green, U. Colorado, *Neuro-Oncology* 2024, PMID 37934854 | 10x v3 scRNA | 19 pHGG, but per the authors' own abstract the cohort is **cortical pHGG–dominant with ~1 DMG** | Human | Yes | None — but DMG content too small to matter |
| **GSE245674** | Zang et al., *BMC Neurol* 2024, PMID 38383423 | scRNA | **1 H3K27M DMG** + 1 anaplastic ependymoma (n=2 total) | Human | Yes | None |
| **GSE202174** | (anti-CD47/RT macrophage study), public 2025-05 | 10x scRNA of sorted phagocytosing vs non-phagocytosing macrophages | DMG + GBM tumor-derived macrophages, n=9; myeloid-only, likely PDX-derived | Human | Yes | None |
| **GSE222398** | Patch-seq of patient-derived DMG xenografts, 2023 | Patch-seq | n=134 cells, DMG PDX (xenograft context) | Human (PDX) | Yes | None |
| **GSE262775** | "Coexisting DMG subpopulation dependencies", public 2025-12, PMID 42020604 | scRNA | SU-DIPG-XVII PDX, 9-drug panel, n=24 samples — single-line drug screen | Human (PDX) | Yes | None |
| **PEDIA-BRAIN**, bioRxiv 2025.09.21.677597 | Wedemeyer/Mardis, Nationwide Children's | snRNA + snATAC multiome | 149,771 normal human pons nuclei (fetal→adult); DMG snRNA mapped onto reference as exemplar | Human | Browsable/downloadable at pediabrain.nchgenomics.org | DMG exemplar data appear to be existing cohorts; value here is the **normal pons reference** |
| **Project HOPE** (Danielli, Neyazi, Hack, Jiang… Filbin lab), *Neuro-Oncol Pediatrics* 2025 abstract, DOI 10.1093/neuped/wuaf001.148 | Filbin lab | "spatiotemporal single-cell" | 76 paired pHGG samples incl. DMG H3K27M | Human | **Not public yet** (conference abstract only; no GEO found) | Unknown — monitor |
| E-MTAB-15147 | Brainstem-organoid DMG (DMGO) + GD2 CAR-T, *Nat Cancer* 2026 | scRNA | H3.3K27M engineered organoid tumors | Human organoid model | Yes (ArrayExpress) | Re-analyzed by Collot et al.; model, not patient |

## 5. Overlap / same-study flags for the six in-use datasets

- **GSE184357 (Liu 2022)**: the series also underlies the paper's *spatial* component, but that component is **HybISS targeted in situ sequencing** (116-gene panel; 16 sections / 14 tumors), **not** transcriptome-wide Visium. No separate Visium accession exists for Liu 2022. Same-study overlap; HybISS is too targeted for spatial validation of transcriptome-wide signatures.
- **GSE210568** = Jessa 2022 **superseries** (scRNA GSE210567 + scATAC). Collot et al. re-analyzed GSE210568 — same data you already use.
- **Collot Zenodo** deposition bundles re-analyzed GSE194329 + GSE280990 + GSE210568 + E-MTAB-15147 alongside new data — use only the newly generated samples to avoid double-counting.

---

## Compatibility assessment

**Suitable as external *spatial* validation (human DMG):**
1. **Collot et al. Zenodo (Visium + CosMx)** — the only substantial *new* human DMG spatial resource. 16 Visium + 23 CosMx patient samples, autopsy/programmatic cohort, open processed data (Space Ranger outputs; note: filtered/processed matrices, not raw FASTQ). Compatible with Visium-based findings; CosMx adds single-cell-resolution cross-platform validation. **Action: download, exclude any spots/samples traceable to GSE194329/GSE280990 re-analysis.**
2. **GSE242352** — only 1 DMG Visium section and processed-only data; usable at most as anecdotal confirmation, not a validation cohort.
3. **Damodharan GeoMx** — ROI-targeted (not unbiased mapping) and **not actually public**; do not plan around it unless authors respond.

**Suitable as *mouse-model* spatial validation:** GSE299434 and GSE246584 (both Xenium, H3.3K27M murine DMG, treated arms) — appropriate for cross-species corroboration of niche/ligand–receptor findings, not for human validation.

**Suitable for *scRNA reference augmentation*:** Okada 2026 snRNA (n=22 postmortem DMG) is the strongest upcoming addition — **deposition not yet visible; re-check GEO or email Diaz/Nazarian labs**. GSE245674 (1 DMG), GSE202174 (myeloid-only), GSE222398/GSE262775 (PDX) are marginal. PEDIA-BRAIN is valuable as a **normal pons snRNA reference** for distinguishing tumor-specific vs regional-normal signals.

**Pooling justified?** Pooling the new Collot Visium samples with the three in-use Visium cohorts is technically feasible (same platform) but **not justified without de-duplication**, because Collot re-analyzed GSE194329/GSE280990 — risk of sample leakage between discovery and validation. Keep Collot strictly as an held-out external set. DeSisto (GSE231859) pooling is not justified for DMG (≈1 DMG among 19 cortical pHGG). Mouse Xenium datasets should never be pooled with human Visium; use for orthogonal confirmation only.

**Bottom line:** beyond the six datasets already in use, there is exactly **one** new human DMG spatial cohort of validation grade (Collot/*Cancer Cell* 2026, Zenodo-hosted Visium+CosMx), two mouse Xenium DMG datasets, and one large upcoming postmortem DMG snRNA cohort (Okada 2026, deposition pending). No Stereo-seq/MERFISH/Slide-seq DMG datasets exist in GEO.

## Search terms used (evidence of diligence)

- GEO E-utilities `db=gds`: `(diffuse midline glioma OR DIPG OR H3K27M OR H3K27-altered) AND (spatial OR Visium OR Xenium OR Slide-seq OR MERFISH OR Stereo-seq OR CosMx OR GeoMx)` → 27 records, all inspected.
- `(Stereo-seq OR MERFISH OR Slide-seq OR CosMx) AND (diffuse midline OR DIPG OR H3K27M)` → **0 records**.
- `(diffuse midline glioma) AND (single cell OR single nucleus) AND 2023:2026[pdat]` → 179 records; all series-level hits screened by title/summary.
- `(GeoMx OR digital spatial profil* OR DSP) AND (glioma OR DMG OR DIPG)` → 241 records; 12 series screened (none DMG-specific/public).
- `(imipridone OR ONC201) AND (DMG OR H3K27) AND (single OR nucleus)`; `postmortem AND diffuse midline glioma`; `(Project HOPE OR Danielli) AND (glioma OR DMG)`; `(DMG OR H3K27) AND (snATAC OR single nucleus OR postmortem) AND 2025:2026[pdat]` → screened; no public Okada/Project HOPE data yet.
- Literature: Consensus/Exa searches "diffuse midline glioma DIPG H3K27M spatial transcriptomics" (2022–2026) and "…single-cell single-nucleus RNA sequencing cohort patient" (2024–2026); targeted web searches for Damodharan GeoMx deposition, IGSF11–VISTA data, CHOP pHGG atlas, DeSisto cohort composition, PEDIA-BRAIN availability.
