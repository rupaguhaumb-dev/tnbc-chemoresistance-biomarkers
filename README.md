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

Analyses 01–03 read a parent-folder workspace `GSE25066_TNBC_MAP3K8_workspace.RData` containing the parsed `clinical_tnbc` (178 × 82) and `expr_tnbc` (22283 × 178) objects. That workspace is **not** committed to this repo (the file is 218 MB and not redistributable in raw form). Re-running 01–03 requires that workspace; re-running 04 only requires internet access (it pulls live from GEO).

## How to run

```r
# inside R, with this folder as the working directory
source("01_pCR_RCB_Analysis/pCR_RCB_Analysis.R")
source("02_FourGroup_Interaction_KM/FourGroup_Interaction_KM.R")
source("03_Continuous_Cox/Continuous_Cox.R")
source("04_GSE41998_Validation/GSE41998_Validation.R")
source("Build_Summary_Doc.R")
```

## R packages used

`survival`, `survminer`, `rms`, `ggplot2`, `ggpubr`, `officer`, `flextable`, `dplyr`, `tidyr`, `GEOquery`, `Biobase`, `httr`, `jsonlite`.

## Key findings (one-line summary)

- MAP3K8 and IL-1β both behave as **risk markers** for distant relapse-free survival in chemo-treated TNBC (high expression → worse DRFS).
- Effect is on **late relapse**, not on initial chemoresponse — neither gene predicts pCR in either GSE25066 or GSE41998.
- Joint effect is **additive, not synergistic** (interaction LRT P = 0.86).
- Both effects survive adjustment for age, stage, and RCB; IL-1β shows a time-varying hazard (PH P = 0.022) consistent with a late-relapse signal.

See `Summary.docx` for the full integrated narrative.
