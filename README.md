# TNBC Chemoresistance Biomarkers — MAP3K8 and IL-1β

A reproducible R-based bioinformatics pipeline that evaluates **MAP3K8** and **IL-1β** as candidate chemoresistance / distant-relapse biomarkers in chemo-treated triple-negative breast cancer (TNBC), using public GEO microarray cohorts.

Each analysis is self-contained in its own subfolder: an R script, the output figures (PDF / PNG), CSV statistics, and a Word `Methods` document with embedded figures, results tables, interpretation, and limitations.

---

## What's in here

| Subfolder | Question | Key outputs |
|---|---|---|
| `01_pCR_RCB_Analysis/` | Do MAP3K8 / IL-1β predict pathological complete response (pCR) or RCB class? | Wilcoxon, Kruskal-Wallis, logistic regression; boxplots |
| `02_FourGroup_Interaction_KM/` | Do MAP3K8 and IL-1β interact, or carry independent prognostic signal? | 4-group Kaplan-Meier (DRFS), pairwise BH log-rank, Cox interaction LRT, group-vs-reference forest |
| `03_Continuous_Cox/` | Does the signal survive without dichotomisation, and after adjustment for age + stage + RCB? | Univariable / bivariable / multivariable continuous Cox (per-SD HR), restricted-cubic-spline HR curves, Schoenfeld PH test |
| `04_GSE41998_Validation/` | Does an independent neoadjuvant-chemo TNBC cohort reproduce the findings? | Pulls GSE41998 via `GEOquery`, restricts to TNBC, repeats the pCR analysis |
| `Summary.docx` | Project-level integrated narrative across all four analyses | — |
| `Build_Summary_Doc.R` | Rebuilds `Summary.docx` from the per-analysis CSV outputs (no model re-fitting) | — |

## Cohorts

- **Discovery:** GSE25066 (Hatzis & Pusztai, *JAMA* 2011; Booser series) — 508 women, neoadjuvant anthracycline-taxane, GPL96 / U133A. TNBC subset *n = 178* (ER- / PR- / HER2-).
- **Validation:** GSE41998 (Horak et al., SWOG) — 279-patient neoadjuvant ixabepilone-vs-paclitaxel trial, GPL96. TNBC subset *n = 140*.

## Data dependency

Analyses 01–03 read a workspace file `GSE25066_TNBC_MAP3K8_workspace.RData` containing the parsed `clinical_tnbc` (178 × 82) and `expr_tnbc` (22283 × 178) objects derived from GEO accession GSE25066. That workspace is **not** committed (218 MB; the underlying GSE25066 raw data is freely available on GEO). Re-running 01–03 requires the workspace; analysis 04 only needs internet access (it pulls live from GEO via `GEOquery`).

By default, scripts look for the workspace at `../../GSE25066_TNBC_MAP3K8_workspace.RData` (i.e., one directory above the repo root). You can override this with an environment variable:

```bash
export TNBC_WORKSPACE=/absolute/path/to/GSE25066_TNBC_MAP3K8_workspace.RData
```

## How to run

```bash
# from a terminal, repo root = the Analyses/ folder:
Rscript 01_pCR_RCB_Analysis/pCR_RCB_Analysis.R
Rscript 02_FourGroup_Interaction_KM/FourGroup_Interaction_KM.R
Rscript 03_Continuous_Cox/Continuous_Cox.R
Rscript 04_GSE41998_Validation/GSE41998_Validation.R
Rscript Build_Summary_Doc.R
```

Each script auto-detects its own location and writes its outputs to the same subfolder. Analysis 04 caches the GSE41998 ExpressionSet in `04_GSE41998_Validation/GSE41998_eset.rds` after the first run (cache file is gitignored).

## R packages used

`survival`, `survminer`, `rms`, `ggplot2`, `ggpubr`, `officer`, `flextable`, `dplyr`, `tidyr`, `GEOquery`, `Biobase`, `httr`, `jsonlite`.

## Key findings (one-line summary)

- MAP3K8 and IL-1β both behave as **risk markers** for distant relapse-free survival in chemo-treated TNBC (high expression → worse DRFS).
- Effect is on **late relapse**, not on initial chemoresponse — neither gene predicts pCR in either GSE25066 or GSE41998.
- Joint effect is **additive, not synergistic** (interaction LRT P = 0.86).
- Both effects survive adjustment for age, stage, and RCB; IL-1β shows a time-varying hazard (PH P = 0.022) consistent with a late-relapse signal.

See `Summary.docx` for the full integrated narrative.
