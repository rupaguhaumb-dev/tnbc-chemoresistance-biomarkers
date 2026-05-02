## ============================================================
##  Build_Summary_Doc.R
##  Builds a top-level Summary.docx that integrates the four
##  per-analysis methods documents into a single project-level
##  narrative for MAP3K8 + IL-1beta in chemo-treated TNBC.
##  Reads the CSV outputs of analyses 01-04 (so re-running this
##  script does not re-fit any model).
## ============================================================

suppressPackageStartupMessages({
  library(officer); library(flextable); library(dplyr)
})

## ---- Portable path setup (this script lives at the Analyses/ root) ------
ROOT <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", a, value = TRUE)
  if (length(f)) return(dirname(normalizePath(sub("^--file=", "", f[1]))))
  fr <- sys.frames()
  for (i in rev(seq_along(fr)))
    if (!is.null(fr[[i]]$ofile))
      return(dirname(normalizePath(fr[[i]]$ofile)))
  getwd()
})
setwd(ROOT)

s01 <- read.csv("01_pCR_RCB_Analysis/pCR_RCB_stats.csv")
s02_summary <- read.csv("02_FourGroup_Interaction_KM/interaction_summary.csv")
s02_int     <- read.csv("02_FourGroup_Interaction_KM/interaction_cox_terms.csv")
s02_g4      <- read.csv("02_FourGroup_Interaction_KM/group4_vs_reference_HR.csv")
s03_univ    <- read.csv("03_Continuous_Cox/cox_univariable.csv")
s03_biv     <- read.csv("03_Continuous_Cox/cox_bivariable.csv")
s03_mv      <- read.csv("03_Continuous_Cox/cox_multivariable.csv")
s03_ph      <- read.csv("03_Continuous_Cox/cox_multivariable_PHtest.csv")
s04         <- read.csv("04_GSE41998_Validation/GSE41998_TNBC_pCR_stats.csv")

doc <- read_docx() |>
  body_add_par("MAP3K8 and IL-1beta as chemoresistance biomarkers in TNBC",
               style = "heading 1") |>
  body_add_par("Project-level summary across four pre-specified analyses on the chemo-treated TNBC cohort GSE25066, with an external validation cohort GSE41998. This document integrates the methods and findings of the four per-analysis Word documents stored in the subfolders 01-04.",
               style = "Normal") |>

  ## ---- Cohorts ----------------------------------------------------------
  body_add_par("Cohorts", style = "heading 2") |>
  body_add_par(
    "Discovery cohort: GSE25066 (Hatzis & Pusztai, JAMA 2011; Booser series). 508 women with primary breast cancer treated with neoadjuvant anthracycline-taxane chemotherapy, profiled on Affymetrix HG-U133A (GPL96). The TNBC subset (n = 178), defined as ER-IHC negative AND PR-IHC negative AND HER2 negative, was used for all primary analyses. Endpoints: distant relapse-free survival (DRFS, time `drfs_even_time_years:ch1` x 12 months, event `drfs_1_event_0_censored:ch1`); pathological complete response (pCR vs RD); Residual Cancer Burden class.",
    style = "Normal") |>
  body_add_par(
    "Validation cohort: GSE41998 (Horak et al.; SWOG / BMS-247550-013), 279 patients in a randomised neoadjuvant ixabepilone-vs-paclitaxel trial, GPL96. The TNBC subset (n = 140) and the pCR endpoint were used; DRFS is not publicly provided in GEO for this series.",
    style = "Normal") |>
  body_add_par(
    "Probes on GPL96: MAP3K8 = 205027_s_at (single probe). IL-1beta has two probes (205067_at and 39402_at); we selected the higher-mean probe within the TNBC subset of each cohort, which yielded 205067_at in GSE25066 and 39402_at in GSE41998.",
    style = "Normal") |>

  ## ---- Important direction-of-effect note ------------------------------
  body_add_par("Important note on direction of effect", style = "heading 2") |>
  body_add_par(
    "In the original MAP3K8 KM analysis (parent script and saved workspace), the reported HR = 0.56 (0.34 - 0.92, log-rank P = 0.02) was computed with the factor MAP3K8_group having levels c('high','low'); R's coxph treats the second level as the contrast level, so this HR is for MAP3K8-low versus MAP3K8-high. Inverting it gives HR(MAP3K8-high vs MAP3K8-low) = 1/0.56 ~ 1.80, i.e. MAP3K8-high tumours have higher risk of distant relapse, the same direction as IL-1beta. All analyses in this folder are presented and interpreted in that corrected direction (HR > 1 means the high group is at higher risk).",
    style = "Normal") |>

  ## ---- Analysis 01 ------------------------------------------------------
  body_add_par("Analysis 01 - pCR / RCB", style = "heading 2") |>
  body_add_par(
    "Question: do MAP3K8 or IL-1beta mRNA levels distinguish patients who achieve pCR from those with residual disease, and do they track RCB severity? Tests: Wilcoxon for pCR vs RD, Kruskal-Wallis across RCB classes, univariable logistic regression with continuous expression.",
    style = "Normal") |>
  body_add_flextable(flextable(s01) |> autofit() |> fontsize(size = 9, part = "all")) |>
  body_add_par("Figure A1.1 - pCR vs RD boxplots (analysis 01)", style = "heading 3") |>
  body_add_img(src = file.path(ROOT, "01_pCR_RCB_Analysis/pCR_boxplots.png"),
               width = 6.3, height = 3.2) |>
  body_add_par("Verdict: neither gene predicts pCR. The signal these genes carry is on relapse, not on initial chemoresponse.",
               style = "Normal") |>

  ## ---- Analysis 02 ------------------------------------------------------
  body_add_par("Analysis 02 - MAP3K8 x IL-1beta four-group KM", style = "heading 2") |>
  body_add_par(
    "Question: do the two genes interact, or do they contribute independent prognostic information? Design: median dichotomisation of each gene; combined into a 4-level factor with MAP3K8_low / IL1B_low as reference.",
    style = "Normal") |>
  body_add_par("Group sizes and overall test:", style = "Normal") |>
  body_add_flextable(flextable(s02_summary) |> autofit() |> fontsize(size = 9, part = "all")) |>
  body_add_par("Group versus reference HRs:", style = "Normal") |>
  body_add_flextable(flextable(s02_g4) |> autofit() |> fontsize(size = 9, part = "all")) |>
  body_add_par("Cox interaction parameterisation:", style = "Normal") |>
  body_add_flextable(flextable(s02_int) |> autofit() |> fontsize(size = 9, part = "all")) |>
  body_add_par("Figure A2.1 - Four-group Kaplan-Meier curves (DRFS)", style = "heading 3") |>
  body_add_img(src = file.path(ROOT, "02_FourGroup_Interaction_KM/KMplot_4groups.png"),
               width = 6.3, height = 6.3) |>
  body_add_par("Figure A2.2 - Group-vs-reference HR forest", style = "heading 3") |>
  body_add_img(src = file.path(ROOT, "02_FourGroup_Interaction_KM/forest_group4_vs_reference.png"),
               width = 6.3, height = 3.2) |>
  body_add_par(
    "Verdict: effects are essentially additive (interaction LRT P = 0.86). The doubly-high group has the highest distant-relapse risk; the doubly-low group is the best-prognosis reference. Pairwise log-rank between best and worst groups was borderline after BH adjustment.",
    style = "Normal") |>

  ## ---- Analysis 03 ------------------------------------------------------
  body_add_par("Analysis 03 - Continuous Cox", style = "heading 2") |>
  body_add_par(
    "Question: does the prognostic signal survive the loss of dichotomisation, and does it survive adjustment for age, stage and RCB? Models: univariable, bivariable and multivariable Cox on standardised log2 expression (HR per 1 SD), plus restricted cubic splines.",
    style = "Normal") |>
  body_add_par("Univariable:", style = "Normal") |>
  body_add_flextable(flextable(s03_univ) |> autofit() |> fontsize(size = 9, part = "all")) |>
  body_add_par("Bivariable:", style = "Normal") |>
  body_add_flextable(flextable(s03_biv) |> autofit() |> fontsize(size = 9, part = "all")) |>
  body_add_par("Multivariable (adjusted for age + stage + RCB):", style = "Normal") |>
  body_add_flextable(flextable(s03_mv) |> autofit() |> fontsize(size = 9, part = "all")) |>
  body_add_par("Schoenfeld residual test:", style = "Normal") |>
  body_add_flextable(flextable(s03_ph) |> autofit() |> fontsize(size = 9, part = "all")) |>
  body_add_par("Figure A3.1 - Univariable + bivariable HR forest", style = "heading 3") |>
  body_add_img(src = file.path(ROOT, "03_Continuous_Cox/forest_univariable.png"),
               width = 6.3, height = 3.2) |>
  body_add_par("Figure A3.2 - MAP3K8 RCS hazard-ratio curve", style = "heading 3") |>
  body_add_img(src = file.path(ROOT, "03_Continuous_Cox/spline_MAP3K8.png"),
               width = 5.4, height = 3.8) |>
  body_add_par("Figure A3.3 - IL-1beta RCS hazard-ratio curve", style = "heading 3") |>
  body_add_img(src = file.path(ROOT, "03_Continuous_Cox/spline_IL1B.png"),
               width = 5.4, height = 3.8) |>
  body_add_par(
    "Verdict: both genes consistently trend in the risk direction (HR ~1.2-1.3 per SD), in the same direction across univariable, bivariable and multivariable specifications, and remain so after adjustment for the strongest known prognostic factor in the cohort (RCB-III HR = 2.67, P = 0.0073). IL-1beta violates the proportional-hazards assumption (P = 0.022), consistent with a late-relapse effect that a single HR underestimates. Splines confirm a monotonic dose-response and rule out U-shaped or threshold patterns.",
    style = "Normal") |>

  ## ---- Analysis 04 ------------------------------------------------------
  body_add_par("Analysis 04 - GSE41998 validation", style = "heading 2") |>
  body_add_par(
    "Question: does an independent neoadjuvant chemo-treated TNBC cohort reproduce the GSE25066 findings? Cohort: GSE41998 TNBC subset (n = 140, GPL96), pCR endpoint only (no DRFS in GEO).",
    style = "Normal") |>
  body_add_flextable(flextable(s04) |> autofit() |> fontsize(size = 9, part = "all")) |>
  body_add_par("Figure A4.1 - GSE41998 TNBC pCR boxplots", style = "heading 3") |>
  body_add_img(src = file.path(ROOT, "04_GSE41998_Validation/GSE41998_TNBC_pCR_boxplots.png"),
               width = 6.3, height = 3.2) |>
  body_add_par(
    "Verdict: the pCR null is independently confirmed - neither gene predicts chemoresponse in GSE41998. The DRFS-positive result of GSE25066 is not testable here because GEO does not expose survival times for GSE41998. Direction of the IL-1beta effect (higher in RD than in pCR) is consistent across cohorts.",
    style = "Normal") |>

  ## ---- Bottom-line ------------------------------------------------------
  body_add_par("Bottom-line conclusions for the project", style = "heading 2") |>
  body_add_par(
    "1. MAP3K8 and IL-1beta mRNA are both risk markers for distant relapse in chemo-treated TNBC. Direction of effect is the same for both genes (HR > 1 with high expression), biologically coherent given their shared NF-kB / IRAK / inflammation axis.",
    style = "Normal") |>
  body_add_par(
    "2. The signal is on late relapse (DRFS), not on initial chemoresponse (pCR / RCB). The two endpoints capture different biology, and these genes capture chemoresistance via residual-disease outgrowth rather than via primary cytotoxic insensitivity.",
    style = "Normal") |>
  body_add_par(
    "3. The effects are independent rather than synergistic - the four-group interaction term is non-significant and HRs barely move on mutual adjustment in the bivariable Cox. A patient who is high for both is at greatest risk, but this is approximately the sum of the two univariable signals.",
    style = "Normal") |>
  body_add_par(
    "4. Both effects survive adjustment for age, stage and RCB - the strongest known prognostic factor in the cohort. RCB-III itself confers HR = 2.67; MAP3K8 and IL-1beta carry information beyond this.",
    style = "Normal") |>
  body_add_par(
    "5. IL-1beta acts on late relapse: PH violation and visual late-curve separation both indicate a time-dependent effect that a single HR understates. A time-varying or stratified Cox is the natural next refinement.",
    style = "Normal") |>
  body_add_par(
    "6. Validation in GSE41998 is partial: the pCR null is reproduced; the DRFS finding is untestable because GEO does not expose survival for that series. Stronger external validation requires DRFS-equipped TNBC cohorts (METABRIC TNBC chemo subset; FUSCC TNBC under controlled access).",
    style = "Normal") |>

  ## ---- Recommended next analyses ---------------------------------------
  body_add_par("Recommended next analyses", style = "heading 2") |>
  body_add_par(
    "(a) Time-stratified or time-varying Cox for IL-1beta to formally quantify the late-relapse hazard. (b) Optimal cutpoint analysis (maxstat / surv_cutpoint) on GSE25066 for both genes, to check whether a non-median threshold sharpens the binary KM. (c) METABRIC TNBC chemo-treated subset for an independent OS / RFS validation. (d) Pathway-level signature (IL1B, IL6, IL8, NLRP3, CASP1, MAP3K8) z-score versus DRFS - cytokine signatures usually outperform single genes. (e) For protein-level support, cite published IHC TMA studies for IL-1beta and MAP3K8 rather than attempting to reanalyse mass-spec data in which IL-1beta is not detected and chemo-treated TNBC cohorts are not available.",
    style = "Normal") |>

  body_add_par("Files in this Analyses folder", style = "heading 2") |>
  body_add_par("01_pCR_RCB_Analysis/        - pCR / RCB chemoresponse analysis", style = "Normal") |>
  body_add_par("02_FourGroup_Interaction_KM/- MAP3K8 x IL1B four-group KM with interaction", style = "Normal") |>
  body_add_par("03_Continuous_Cox/          - continuous and multivariable Cox + splines", style = "Normal") |>
  body_add_par("04_GSE41998_Validation/     - independent neoadjuvant TNBC validation", style = "Normal") |>
  body_add_par("Each subfolder contains its own R script, data outputs (CSV), figures (PDF/PNG/TIFF) and a Methods .docx.",
               style = "Normal") |>

  body_add_par("Reproducibility", style = "heading 2") |>
  body_add_par(paste0("R version: ", R.version.string,
                      ". Platform: ", R.version$platform,
                      ". Date: ", Sys.Date(),
                      ". Source workspace: GSE25066_TNBC_MAP3K8_workspace.RData."),
               style = "Normal")

print(doc, target = "Summary.docx")
cat("Wrote: ", file.path(ROOT, "Summary.docx"), "\n", sep = "")
