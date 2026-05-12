# 07 — MAP3K8 as a chemoresistance / chemosensitivity predictor

Evaluates MAP3K8 expression in the GSE25066 TNBC cohort (n = 178)
against five different chemoresponse readouts:

| Endpoint | Column | Outcome |
|---|---|---|
| Actual pCR vs RD | `pathologic_response_pcr_rd:ch1` | pCR (57), RD (113) |
| Actual RCB class | `pathologic_response_rcb_class:ch1` | RCB-0/I, RCB-II, RCB-III |
| Hatzis genomic signature | `chemosensitivity_prediction:ch1` | Rx Sensitive vs Rx Insensitive |
| DLDA-30 prediction | `dlda30_prediction:ch1` | pCR vs RD |
| RCB-0/I prediction | `rcb_0_i_prediction:ch1` | RCB-0/I vs RCB-II/III |

## What it computes

1. **Distribution tests** — Wilcoxon (2-class endpoints) and Kruskal-Wallis (RCB class).
2. **Univariable logistic regression** — pCR ~ MAP3K8 (continuous log2).
3. **ROC analysis** — AUC for MAP3K8 alone, the Hatzis signature alone, and a combined logistic model; DeLong test for incremental value of MAP3K8.
4. **Cutoff performance** — sensitivity, specificity, PPV, NPV at the median cut and at the Youden-optimal cut.
5. **Concordance** — 2×2 cross-tab + Cohen's kappa between MAP3K8 dichot and the Hatzis Rx-Sensitive / Insensitive call.

## Outputs

```
boxplots_MAP3K8_by_class.{pdf,png}        # 5-panel composite (A-E)
ROC_MAP3K8_vs_signature.{pdf,png}         # MAP3K8 vs Hatzis vs combined
MAP3K8_chemoresponse_stats_2class.csv      # Wilcoxon per binary endpoint
MAP3K8_chemoresponse_stats_RCB.csv         # Kruskal-Wallis across RCB classes
MAP3K8_logistic_pCR.csv                    # Logistic regression coefficients
MAP3K8_classifier_metrics.csv              # Performance at median + Youden cutoffs
MAP3K8_vs_Hatzis_concordance.csv           # 2x2 table + kappa
MAP3K8_chemoresponse_SUMMARY.csv           # One-line summary metrics
```

## Run

```bash
Rscript MAP3K8_Chemoresponse.R
```

Needs the parent `GSE25066_TNBC_MAP3K8_workspace.RData` (see top-level `README.md` for the `TNBC_WORKSPACE` env-var contract).
