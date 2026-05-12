# 08 — MAP3K8 as a chemoresistance prognostic marker in chemo-treated TNBC

## Hypothesis

Among triple-negative breast cancer patients treated with neoadjuvant
anthracycline-taxane chemotherapy, MAP3K8-high tumours carry a
**chemoresistance phenotype**: they relapse distantly after chemotherapy
at a higher rate than MAP3K8-low tumours, regardless of their initial
pathological response.

This complements Analysis 07 (which showed MAP3K8 does **not** predict
initial chemoresponse / pCR). The two analyses together establish
MAP3K8 as a **prognostic marker of post-chemo relapse**, not a
**predictive marker of primary chemoresponse**.

## Three converging lines of evidence

| | Test | What it demonstrates |
|---|---|---|
| **A** | DRFS in the full cohort, MAP3K8-high vs -low | The base prognostic signal |
| **B** | Multivariable Cox: DRFS ~ MAP3K8 + age + AJCC stage + RCB | The MAP3K8 effect is **independent of initial chemoresponse** (RCB) |
| **C** | Pre-specified subgroup KM in **RCB-0/I (chemo responders)** | Even among patients who responded well, MAP3K8-high still marks relapse — the strongest direct test of the chemoresistance hypothesis |
| **D** | Subgroup KM in RCB-II/III (residual disease) | Whether MAP3K8 further stratifies already-resistant tumours |

## Outputs

```
KM_MAP3K8_overall.{pdf,png,tiff}     # Full-cohort Kaplan-Meier
KM_MAP3K8_RCB0I.{pdf,png}            # RCB-0/I subgroup KM
KM_MAP3K8_RCBII_III.{pdf,png}        # RCB-II/III subgroup KM
forest_subgroup_HRs.{pdf,png}        # Forest plot: HR across subgroups
cumulative_relapse_rates.csv         # 2-/5-/7-yr relapse rates by group
multivariable_Cox.csv                # Full multivariable HR table
subgroup_HR_table.csv                # Per-subgroup HR / CI / P
SUMMARY.csv                          # One-line per metric
```

## Endpoint

Distant relapse-free survival (DRFS):
- time = `drfs_even_time_years:ch1` × 12 months
- event = `drfs_1_event_0_censored:ch1`

## Run

```bash
Rscript MAP3K8_Chemoresistance.R
```

Needs the parent `GSE25066_TNBC_MAP3K8_workspace.RData` (see top-level `README.md` for the `TNBC_WORKSPACE` env-var contract).
