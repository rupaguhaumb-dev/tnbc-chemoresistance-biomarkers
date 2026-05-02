# 05 — Single-gene Kaplan-Meier curves

Standalone median-split DRFS Kaplan-Meier curves for each gene
individually, complementing the joint four-group analysis in folder
`02_FourGroup_Interaction_KM/`.

| Script | Output |
|---|---|
| `MAP3K8_TNBC_KM.R` | `KMplot_MAP3K8.{pdf,png,tiff}` + `MAP3K8_KM_stats.csv` |
| `IL1B_TNBC_KM.R`   | `KMplot_IL1B.{pdf,png,tiff}` + `IL1B_KM_stats.csv` |

Both scripts use the median expression in the TNBC cohort as the cutoff
and report HR for **high vs low** (low = reference). DRFS time and
event are pulled from the same clinical columns as analyses 01–03.

## Run

```bash
Rscript MAP3K8_TNBC_KM.R
Rscript IL1B_TNBC_KM.R
```

These need the parent `GSE25066_TNBC_MAP3K8_workspace.RData` workspace
(see top-level `README.md` for the `TNBC_WORKSPACE` env-var contract).
