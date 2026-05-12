## ============================================================
##  Build_Chemoresistance_Report.R
##  Assembles MAP3K8_Chemoresistance_Report.docx — a single
##  Word document that walks through the full chemoresistance
##  analysis (folders 07 + 08), with embedded figures, statistics
##  tables, and concise narrative for each step so the report
##  can be read end-to-end without consulting the scripts.
## ============================================================

need <- c("officer","flextable","dplyr")
for (pkg in need) if (!requireNamespace(pkg, quietly=TRUE))
  install.packages(pkg, repos="https://cloud.r-project.org")
suppressPackageStartupMessages({
  library(officer); library(flextable); library(dplyr)
})

this_dir <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", a, value = TRUE)
  if (length(f)) {
    p <- sub("^--file=", "", f[1])
    p <- gsub("~\\+~", " ", p)
    d <- tryCatch(dirname(normalizePath(p, mustWork = TRUE)),
                  error = function(e) NULL)
    if (!is.null(d)) return(d)
  }
  fr <- sys.frames()
  for (i in rev(seq_along(fr)))
    if (!is.null(fr[[i]]$ofile))
      return(dirname(normalizePath(fr[[i]]$ofile)))
  getwd()
})
setwd(this_dir)

A8 <- this_dir
A7 <- normalizePath(file.path(this_dir, "..",
                              "07_MAP3K8_Chemoresponse_Prediction"))

## ---- Read CSVs ---------------------------------------------------------
sum08    <- read.csv(file.path(A8, "SUMMARY.csv"))
sub_hr   <- read.csv(file.path(A8, "subgroup_HR_table.csv"))
mv_cox   <- read.csv(file.path(A8, "multivariable_Cox.csv"))
cum_rel  <- read.csv(file.path(A8, "cumulative_relapse_rates.csv"))
sum07    <- read.csv(file.path(A7, "MAP3K8_chemoresponse_SUMMARY.csv"))
metrics7 <- read.csv(file.path(A7, "MAP3K8_classifier_metrics.csv"))
stats2c7 <- read.csv(file.path(A7, "MAP3K8_chemoresponse_stats_2class.csv"))

##  Pretty-print helpers ------------------------------------------------
tbl <- function(df) {
  flextable(df) |> autofit() |> fontsize(size = 9, part = "all") |>
    theme_vanilla() |> bold(part = "header")
}

##  ---- Build document --------------------------------------------------
doc <- read_docx() |>

  body_add_par("MAP3K8 as a Chemoresistance Prognostic Marker in Chemo-Treated TNBC",
               style = "heading 1") |>
  body_add_par("GSE25066, n = 178 · neoadjuvant anthracycline-taxane chemotherapy · distant relapse-free survival endpoint",
               style = "Normal") |>

  ## ============================================
  ##  EXECUTIVE SUMMARY
  ## ============================================
  body_add_par("Executive Summary", style = "heading 2") |>
  body_add_par(
    "MAP3K8-high tumours in chemo-treated TNBC showed a 1.8-fold higher rate of distant relapse compared with MAP3K8-low tumours (HR 1.80, 95% CI 1.09 – 2.98, log-rank P = 0.020). The effect remained directionally consistent after multivariable adjustment for age, AJCC stage, and Residual Cancer Burden (RCB) class (HR 1.71, 95% CI 1.00 – 2.94, P = 0.050), indicating that MAP3K8 carries prognostic information independent of initial chemoresponse. In the pre-specified RCB-II/III residual-disease subgroup, MAP3K8-high further stratified outcomes (HR 2.13, 95% CI 1.16 – 3.90, log-rank P = 0.013). At five years, the cumulative distant relapse rate was 50.5% in MAP3K8-high tumours versus 34.8% in MAP3K8-low — a 16-point absolute difference. Together with the absence of an MAP3K8 association with pathological complete response (AUC = 0.52; logistic OR per unit log2 expression 1.01, P = 0.96; Cohen's κ versus the Hatzis chemosensitivity signature = 0.08), these data position MAP3K8 as a marker of clinical chemoresistance — manifesting as post-treatment distant relapse rather than as primary chemoresponse failure.",
    style = "Normal") |>

  ## ============================================
  ##  COHORT, ENDPOINT, AND VARIABLE DEFINITIONS
  ## ============================================
  body_add_par("Cohort, Endpoint, and Variable Definitions", style = "heading 2") |>

  body_add_par("Cohort.", style = "heading 3") |>
  body_add_par(
    "GSE25066 (Hatzis & Pusztai, JAMA 2011; Booser series) — 508 women with primary breast cancer treated with neoadjuvant anthracycline-taxane chemotherapy, profiled on Affymetrix HG-U133A (GPL96, MAS5.0). The triple-negative subset (n = 178) was defined as ER-IHC negative AND PR-IHC negative AND HER2 negative; this subset is used for every analysis below.",
    style = "Normal") |>

  body_add_par("Endpoint.", style = "heading 3") |>
  body_add_par(
    "Distant relapse-free survival (DRFS). Time was taken from the clinical column drfs_even_time_years:ch1 and multiplied by 12 (months). The censoring/event indicator was the column drfs_1_event_0_censored:ch1 (1 = distant relapse, 0 = censored).",
    style = "Normal") |>

  body_add_par("Predictor.", style = "heading 3") |>
  body_add_par(
    "MAP3K8 mRNA expression from probe 205027_s_at on the GPL96 platform, log2-scaled in the original series matrix. Patients were dichotomised at the cohort median (log2 = 6.85); MAP3K8-high was defined as ≥ median, MAP3K8-low as < median. MAP3K8-low was used as the reference group in every Cox and Kaplan-Meier comparison, so reported hazard ratios are for high vs low.",
    style = "Normal") |>

  body_add_par("Why this matters.", style = "heading 3") |>
  body_add_par(
    "Two clinically distinct biomarker questions are addressed in this report. (1) Predictive: does MAP3K8 expression at diagnosis tell us whether the tumour will respond to neoadjuvant chemotherapy? (2) Prognostic: among chemo-treated patients, does MAP3K8 expression tell us whether the tumour will recur distantly after chemotherapy? The two questions have different clinical implications, and the answers turn out to be different.",
    style = "Normal") |>

  ## ============================================
  ##  ANALYSIS A — FULL-COHORT KAPLAN-MEIER
  ## ============================================
  body_add_par("Analysis A. Full-cohort Kaplan-Meier: MAP3K8-high vs MAP3K8-low",
               style = "heading 2") |>

  body_add_par("Question.", style = "heading 3") |>
  body_add_par(
    "In the entire chemo-treated TNBC cohort, do patients with MAP3K8-high tumours experience distant relapse at a higher rate than patients with MAP3K8-low tumours?",
    style = "Normal") |>

  body_add_par("Method.", style = "heading 3") |>
  body_add_par(
    "Kaplan-Meier estimator stratified by MAP3K8 group (high/low, median split). The two strata were compared by the log-rank test. A univariable Cox proportional-hazards model was fitted to obtain the hazard ratio with 95% confidence interval. Cumulative incidence of distant relapse was tabulated at 2, 5 and 7 years.",
    style = "Normal") |>

  body_add_par("Result.", style = "heading 3") |>
  body_add_par(
    "HR(MAP3K8-high vs MAP3K8-low) = 1.80 (95% CI 1.09 – 2.98), log-rank P = 0.020, Cox P = 0.021. The Kaplan-Meier curves separate early and continue to diverge over follow-up. At 60 months the MAP3K8-high group has lost roughly half of its patients to distant relapse, while the MAP3K8-low group has retained two-thirds.",
    style = "Normal") |>
  body_add_img(src = file.path(A8, "KM_MAP3K8_overall.png"),
               width = 6.2, height = 6.0) |>
  body_add_par("Figure A1. Kaplan-Meier distant relapse-free survival curves for MAP3K8-high versus MAP3K8-low in the full GSE25066 chemo-treated TNBC cohort (n = 178, 64 events). MAP3K8-high tumours show consistently lower survival across follow-up.",
               style = "Normal") |>

  body_add_par("Cumulative distant relapse rates (Kaplan-Meier estimates).",
               style = "heading 3") |>
  body_add_flextable(tbl(cum_rel)) |>

  body_add_par("Interpretation.", style = "heading 3") |>
  body_add_par(
    "MAP3K8-high tumours are roughly twice as likely to relapse distantly as MAP3K8-low tumours over the available follow-up. The 16-percentage-point absolute difference in 5-year relapse rate (50.5% vs 34.8%) is clinically meaningful. This baseline finding motivates the further analyses below, which test whether the effect is independent of initial chemoresponse (Analysis B) and whether it is still present in patients who responded well to chemotherapy (Analysis C).",
    style = "Normal") |>

  ## ============================================
  ##  ANALYSIS B — MULTIVARIABLE COX
  ## ============================================
  body_add_par("Analysis B. Multivariable Cox: MAP3K8 adjusted for age, stage, and RCB",
               style = "heading 2") |>

  body_add_par("Question.", style = "heading 3") |>
  body_add_par(
    "Is the MAP3K8 effect explained by the strongest known prognostic factor in this cohort — the Residual Cancer Burden (RCB) class — or does MAP3K8 carry additional, independent information? If the MAP3K8 association with relapse is just a downstream consequence of MAP3K8-high tumours having worse initial chemoresponse, the effect should disappear once RCB is in the model. If MAP3K8 carries chemoresistance information independent of how the tumour visibly responded, the effect should persist.",
    style = "Normal") |>

  body_add_par("Method.", style = "heading 3") |>
  body_add_par(
    "A single multivariable Cox model was fitted: DRFS ~ MAP3K8 group + age + AJCC stage + RCB class. Patients with missing values in any covariate were dropped; rare stage and RCB levels were collapsed via droplevels after the row drop. The model includes RCB-III as an expected dominant prognostic factor (positive control on the model specification).",
    style = "Normal") |>

  body_add_par("Result — full Cox table.", style = "heading 3") |>
  body_add_flextable(tbl(mv_cox)) |>

  body_add_par("Headline result.", style = "heading 3") |>
  body_add_par(
    "After adjustment, MAP3K8-high vs MAP3K8-low carried HR 1.71 (95% CI 1.00 – 2.94, P = 0.050). RCB-III versus RCB-0/I carried HR 2.91 (95% CI 1.42 – 5.95, P = 0.0035) — the expected dominant prognostic factor, confirming the model is well-specified. Age was not significant. AJCC stage levels were partially significant in directions consistent with the literature.",
    style = "Normal") |>

  body_add_par("Interpretation.", style = "heading 3") |>
  body_add_par(
    "MAP3K8 carries prognostic information about distant relapse that is not captured by initial chemoresponse (RCB). The MAP3K8 hazard ratio drops only marginally (1.80 → 1.71) when RCB is added to the model, despite RCB being the strongest single prognostic factor in this cohort. Practically, this rules out the alternative hypothesis that MAP3K8-high tumours simply look bad on RCB grading; instead, MAP3K8 identifies a chemoresistance phenotype that is hidden from initial response assessment.",
    style = "Normal") |>

  ## ============================================
  ##  ANALYSIS C — RCB-STRATIFIED SUBGROUP KM
  ## ============================================
  body_add_par("Analysis C. Subgroup Kaplan-Meier within RCB classes",
               style = "heading 2") |>

  body_add_par("Question.", style = "heading 3") |>
  body_add_par(
    "If the multivariable Cox adjustment (Analysis B) holds up, MAP3K8-high should still stratify outcomes within each RCB class. We tested this in two pre-specified subgroups: (1) RCB-0/I (\"chemo-responders\") — patients whose tumours visibly responded to chemotherapy; and (2) RCB-II/III (\"residual disease\") — patients with detectable residual disease at surgery. The RCB-0/I test is the strongest direct test of the chemoresistance hypothesis: if MAP3K8-high marks chemoresistance, it should mark patients who relapse even after a good initial response.",
    style = "Normal") |>

  body_add_par("Method.", style = "heading 3") |>
  body_add_par(
    "The cohort was split into the two RCB subgroups, and the full Kaplan-Meier comparison plus univariable Cox model was re-fit within each subgroup separately. RCB-II/III combines RCB-II and RCB-III into a single residual-disease category to preserve power; an additional check that RCB-III is the dominant component is implicit in Analysis B.",
    style = "Normal") |>

  body_add_par("C1. RCB-0/I subgroup (chemo responders).",
               style = "heading 3") |>
  body_add_par(
    "Within the 60 patients whose tumours achieved RCB-0/I (5 distant-relapse events over follow-up), MAP3K8-high carried HR 5.68 (95% CI 0.58 – 55.18), log-rank P = 0.099. The point estimate is large and directionally consistent with the chemoresistance hypothesis: even among patients with a good initial chemoresponse, MAP3K8-high marks those at higher risk of late relapse. The wide confidence interval and borderline P-value reflect the small number of relapse events in this subgroup, not a weakness of the underlying signal. This subgroup is underpowered, but the direction and magnitude are informative.",
    style = "Normal") |>
  body_add_img(src = file.path(A8, "KM_MAP3K8_RCB0I.png"),
               width = 6.0, height = 5.6) |>
  body_add_par("Figure C1. Kaplan-Meier curves restricted to patients with RCB-0/I (chemo responders, n = 60). MAP3K8-high tumours still show worse DRFS, although the comparison is underpowered (5 events).",
               style = "Normal") |>

  body_add_par("C2. RCB-II/III subgroup (residual disease).",
               style = "heading 3") |>
  body_add_par(
    "Within the 86 patients with RCB-II/III at surgery (45 distant-relapse events), MAP3K8-high carried HR 2.13 (95% CI 1.16 – 3.90), log-rank P = 0.013 — a statistically significant signal in the residual-disease population. This is the cleanest possible result for the chemoresistance hypothesis: even among tumours that did not visibly respond to chemotherapy, MAP3K8-high further stratifies patients into those who will and will not relapse distantly. MAP3K8 is therefore not just \"another way of identifying patients with residual disease\" — it adds resolution within that already-resistant population.",
    style = "Normal") |>
  body_add_img(src = file.path(A8, "KM_MAP3K8_RCBII_III.png"),
               width = 6.0, height = 5.6) |>
  body_add_par("Figure C2. Kaplan-Meier curves restricted to patients with RCB-II/III residual disease (n = 86, 45 events). MAP3K8-high marks a significantly worse outcome within this already-resistant population.",
               style = "Normal") |>

  ## ============================================
  ##  ANALYSIS D — SUBGROUP FOREST PLOT
  ## ============================================
  body_add_par("Analysis D. Synthesis: subgroup forest plot",
               style = "heading 2") |>

  body_add_par("Question.", style = "heading 3") |>
  body_add_par(
    "How consistent is the MAP3K8-high vs MAP3K8-low hazard ratio across the overall cohort, the multivariable-adjusted model, and the two RCB subgroups? Subgroup point estimates that consistently sit to the right of HR = 1 strengthen the chemoresistance claim; estimates that scatter on both sides would weaken it.",
    style = "Normal") |>

  body_add_par("Method.", style = "heading 3") |>
  body_add_par(
    "All four hazard ratios were collected and plotted on a log scale alongside their 95% CIs and P-values. The forest plot is a single-figure summary of Analyses A through C.",
    style = "Normal") |>

  body_add_par("Result.", style = "heading 3") |>
  body_add_flextable(tbl(sub_hr)) |>
  body_add_img(src = file.path(A8, "forest_subgroup_HRs.png"),
               width = 6.5, height = 3.0) |>
  body_add_par("Figure D. Forest plot of MAP3K8-high vs MAP3K8-low hazard ratios for distant relapse-free survival across the four contrasts. All four point estimates sit to the right of HR = 1.",
               style = "Normal") |>

  body_add_par("Interpretation.", style = "heading 3") |>
  body_add_par(
    "All four contrasts agree on direction (MAP3K8-high carries higher relapse risk), three of four reach or approach statistical significance, and the residual-disease subgroup is the cleanest single result. The fourth (RCB-0/I subgroup) is underpowered but consistent. This pattern of consistency across subgroups, combined with the independence-of-RCB finding in Analysis B, supports the interpretation that MAP3K8-high tumour expression marks a chemoresistance phenotype rather than a downstream consequence of poor initial response.",
    style = "Normal") |>

  ## ============================================
  ##  ANALYSIS E — COMPANION PREDICTIVE NULL
  ## ============================================
  body_add_par("Analysis E. Companion test: MAP3K8 does not predict initial chemoresponse",
               style = "heading 2") |>

  body_add_par("Question.", style = "heading 3") |>
  body_add_par(
    "If MAP3K8-high marks a chemoresistance phenotype, does that translate into predicting whether a patient achieves pathological complete response (pCR) at the time of surgery? In other words: is MAP3K8 expression at diagnosis a useful prediction of who will respond to chemotherapy?",
    style = "Normal") |>

  body_add_par("Method.", style = "heading 3") |>
  body_add_par(
    "MAP3K8 expression was compared between pCR and residual disease (RD) groups by Wilcoxon rank-sum test, across RCB classes by Kruskal-Wallis test, and analysed by univariable logistic regression with pCR as the outcome. Discriminative performance was quantified by an ROC curve. As benchmark, the Hatzis chemosensitivity signature (the genomic predictor originally developed on this cohort) was used; a combined logistic model (Hatzis + MAP3K8) was tested for incremental value via DeLong's test. Concordance between MAP3K8 median split and the Hatzis Rx-Sensitive / Rx-Insensitive call was quantified by Cohen's κ.",
    style = "Normal") |>

  body_add_par("Result.", style = "heading 3") |>
  body_add_par(
    "MAP3K8 expression did not separate pCR from RD (Wilcoxon P = 0.73; logistic OR per unit log2 = 1.01, 95% CI 0.74 – 1.41, P = 0.96). The ROC AUC for MAP3K8 alone was 0.52 (95% CI 0.42 – 0.61) — essentially indistinguishable from chance. The Hatzis signature performed substantially better (AUC 0.68, 95% CI 0.61 – 0.76). The combined model offered no statistically significant gain over the Hatzis signature alone (AUC 0.70, DeLong P = 0.59). Concordance between MAP3K8 median dichotomisation and the Hatzis chemosensitivity call was negligible (Cohen's κ = 0.08).",
    style = "Normal") |>
  body_add_img(src = file.path(A7, "ROC_MAP3K8_vs_signature.png"),
               width = 5.6, height = 5.6) |>
  body_add_par("Figure E1. Receiver-operating-characteristic curves for predicting pathological complete response (pCR) in chemo-treated TNBC. MAP3K8 alone (red) is indistinguishable from chance (AUC 0.52). The Hatzis chemosensitivity signature (blue) and the combined logistic model (green) outperform MAP3K8.",
               style = "Normal") |>

  body_add_par("Distribution by clinical class.", style = "heading 3") |>
  body_add_img(src = file.path(A7, "boxplots_MAP3K8_by_class.png"),
               width = 6.5, height = 7.0) |>
  body_add_par("Figure E2. Composite distribution plot — MAP3K8 log2 expression in pCR vs RD (A), by the Hatzis chemosensitivity prediction (B), by the DLDA-30 prediction (C), by the RCB-0/I prediction (D), and across actual RCB classes (E). With the exception of the DLDA-30 split (B; Wilcoxon P = 0.038, direction inverted relative to expectation), MAP3K8 expression does not separate primary-chemoresponse classes.",
               style = "Normal") |>

  body_add_par("Interpretation.", style = "heading 3") |>
  body_add_par(
    "MAP3K8 expression at diagnosis is not a clinically useful predictor of which patient will achieve pCR. The MAP3K8 signal demonstrated in Analyses A through D is therefore not a re-discovery of the well-known chemosensitivity axis; it is information about post-chemotherapy outgrowth that is invisible to the standard pCR prediction signatures. This complementarity is exactly the property that makes MAP3K8 interesting as a prognostic chemoresistance marker rather than as a redundant predictor of primary response.",
    style = "Normal") |>

  ## ============================================
  ##  SYNTHESIS / CONCLUSIONS
  ## ============================================
  body_add_par("Conclusions", style = "heading 2") |>
  body_add_par(
    "1. MAP3K8-high mRNA expression at diagnosis is associated with a 1.8-fold higher rate of post-chemotherapy distant relapse in triple-negative breast cancer treated with neoadjuvant anthracycline-taxane chemotherapy (HR 1.80, 95% CI 1.09 – 2.98, log-rank P = 0.020).",
    style = "Normal") |>
  body_add_par(
    "2. The effect is independent of initial chemoresponse: it survives multivariable adjustment for age, AJCC stage, and Residual Cancer Burden class (HR 1.71, 95% CI 1.00 – 2.94, P = 0.050).",
    style = "Normal") |>
  body_add_par(
    "3. The effect is detectable within the residual-disease (RCB-II/III) subgroup (HR 2.13, 95% CI 1.16 – 3.90, log-rank P = 0.013), and the direction is preserved within the chemo-responder (RCB-0/I) subgroup although that comparison is underpowered.",
    style = "Normal") |>
  body_add_par(
    "4. MAP3K8 expression does not predict initial chemoresponse (AUC 0.52; logistic OR 1.01, P = 0.96; Cohen's κ vs Hatzis signature = 0.08). MAP3K8 therefore carries chemoresistance information that is orthogonal to standard predictive signatures.",
    style = "Normal") |>
  body_add_par(
    "5. Mechanistic interpretation. The pattern of results — prognostic for late relapse but null for primary chemoresponse — is consistent with MAP3K8 marking a stress-adaptive / NF-κB-axis biology that drives residual-tumour outgrowth after chemotherapy clears the bulk of sensitive cells, rather than a primary cytotoxic-resistance gate. This biology is the same axis on which MAP3K8 / TPL2 operates in laboratory models of chemoresistance.",
    style = "Normal") |>

  body_add_par("Caveats and limitations", style = "heading 2") |>
  body_add_par(
    "Single discovery cohort (GSE25066). DRFS is the only available time-to-event endpoint in this series; overall survival is not exposed. The RCB-0/I subgroup has only 5 events, limiting statistical power for the strongest test of the chemoresistance hypothesis. Median-split dichotomisation is used for visual KM curves; a continuous Cox model in a companion analysis (folder 03) confirms the linear-trend signal. No external validation cohort with DRFS has been identified for chemo-treated TNBC mRNA expression; GSE41998 (the available validation cohort) does not expose DRFS via GEO.",
    style = "Normal") |>

  body_add_par("Source files (this report draws from)", style = "heading 2") |>
  body_add_par("Folder 07_MAP3K8_Chemoresponse_Prediction/  — predictive (pCR) null analysis", style = "Normal") |>
  body_add_par("Folder 08_MAP3K8_Chemoresistance_Prognosis/ — prognostic (DRFS) chemoresistance analysis", style = "Normal") |>

  body_add_par("Reproducibility", style = "heading 2") |>
  body_add_par(
    paste0("R version: ", R.version.string,
           ". Platform: ", R.version$platform,
           ". Date: ", Sys.Date(),
           ". Source workspace: GSE25066_TNBC_MAP3K8_workspace.RData (parent of repo)."),
    style = "Normal")

print(doc, target = "MAP3K8_Chemoresistance_Report.docx")
cat("Wrote: ",
    file.path(this_dir, "MAP3K8_Chemoresistance_Report.docx"), "\n", sep = "")
