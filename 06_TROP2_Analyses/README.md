# 06 — TROP2 / TACSTD2 analyses

DRFS Kaplan-Meier curves and a four-group interaction analysis with
MAP3K8, in the same chemo-treated TNBC cohort (GSE25066, n = 178)
used in folders 01–05.

| Script | Outputs |
|---|---|
| `TROP2_TNBC_KM.R` | Single-gene median-split DRFS KM for TROP2 → `KMplot_TROP2.{pdf,png,tiff}` + `TROP2_KM_stats.csv` |
| `TROP2_MAP3K8_FourGroup_KM.R` | TROP2 × MAP3K8 four-group KM with pairwise BH log-rank, Cox interaction LRT, and group-vs-reference forest → `KMplot_TROP2_MAP3K8_4groups.{pdf,png,tiff}` + `forest_TROP2_MAP3K8_4groups.{pdf,png}` + four CSVs |

## Probe selection

TROP2 is annotated as **TACSTD2** on the Affymetrix HG-U133A platform
(GPL96), with three available probes (202285_s_at, 202286_s_at,
202287_s_at). Following the convention used in folder 05 for IL-1β,
the probe with the highest mean log2 expression in the TNBC cohort
is selected at run time.

## Run

```bash
Rscript TROP2_TNBC_KM.R
Rscript TROP2_MAP3K8_FourGroup_KM.R
```

Both scripts need the parent `GSE25066_TNBC_MAP3K8_workspace.RData`
(see top-level `README.md` for the `TNBC_WORKSPACE` env-var contract).
