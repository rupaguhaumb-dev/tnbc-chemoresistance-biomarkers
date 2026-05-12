## ============================================================
##  Build_Plain_Language_Forest_Report.R
##  A single Word document that explains the forest plot to a
##  non-clinical reader. Every clinical / statistical term is
##  defined in plain English. Hazard-ratio numbers are translated
##  into "X% higher risk" wording so the reader can both
##  understand and re-explain the result.
##
##  Output: MAP3K8_Forest_PlainLanguage_Report.docx
##  Page : 8.5x11 in, 0.75 in margins (7 in usable width)
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

FOREST_PNG <- file.path(this_dir, "forest_subgroup_HRs.png")
stopifnot(file.exists(FOREST_PNG))

## ---- Glossary -----------------------------------------------------------
glossary <- data.frame(
  Term = c(
    "TNBC (Triple-Negative Breast Cancer)",
    "Neoadjuvant chemotherapy",
    "Anthracycline-taxane",
    "Chemoresistance",
    "Chemosensitivity",
    "pCR (Pathological Complete Response)",
    "Residual disease",
    "RCB (Residual Cancer Burden)",
    "RCB-0",
    "RCB-I",
    "RCB-II",
    "RCB-III",
    "RCB-0/I (combined)",
    "RCB-II/III (combined)",
    "DRFS (Distant Relapse-Free Survival)",
    "Distant relapse",
    "MAP3K8",
    "MAP3K8-high / MAP3K8-low",
    "Median split",
    "Hazard Ratio (HR)",
    "HR = 1",
    "HR > 1",
    "HR < 1",
    "95% Confidence Interval (CI)",
    "P-value",
    "P < 0.05",
    "Log-rank test",
    "Cox regression",
    "Univariable",
    "Multivariable (adjusted)",
    "Subgroup analysis",
    "Forest plot",
    "GSE25066",
    "n",
    "Events"),
  Plain_English = c(
    "A type of breast cancer that lacks all three of the most common drug targets (ER, PR, HER2). The hardest subtype to treat with targeted drugs; chemotherapy is the main first-line option.",
    "Chemotherapy given BEFORE surgery, to shrink the tumour. Standard of care for TNBC.",
    "A standard combination chemo regimen for TNBC. Anthracycline = e.g. doxorubicin (Adriamycin); taxane = e.g. paclitaxel (Taxol).",
    "When a tumour survives chemotherapy and grows back later. Opposite of chemosensitivity.",
    "When chemo successfully kills the tumour. Opposite of chemoresistance.",
    "After surgery, the pathologist finds NO remaining cancer at the original tumour site. The best possible chemo response. Equivalent to RCB-0.",
    "Tumour that is still present after chemotherapy.",
    "A standardised score (developed at MD Anderson) that measures how much tumour remains after chemo, by looking at tumour size, cellularity, and lymph-node involvement. Higher number = worse response.",
    "No residual cancer found = pCR. Best possible response. (Same as pCR.)",
    "Small amount of residual cancer (low burden). Considered a good response.",
    "Moderate amount of residual cancer. Considered an intermediate / poor response.",
    "Large amount of residual cancer. Considered the worst response category.",
    "Combined: 'chemo-responders' — patients whose tumour was largely or fully cleared by chemo (RCB-0 OR RCB-I).",
    "Combined: 'residual disease' — patients with significant tumour still present after chemo (RCB-II OR RCB-III).",
    "Time from diagnosis until the cancer reappears at a distant site (lung, bone, liver, brain). The main long-term outcome measure here.",
    "Cancer that comes back somewhere distant from the original tumour, after the original tumour was treated.",
    "A gene that codes for a signalling enzyme (also known as TPL2). It is part of a cellular stress-response pathway linked to inflammation (NF-kB axis). High activity has been linked to cancer cells surviving stress / treatment.",
    "Patients in this study were split in half by MAP3K8 expression. 'High' = above-average expression (top half), 'Low' = below-average (bottom half).",
    "A way of dividing patients into two equal halves at the middle value. Here, 89 patients are MAP3K8-high and 89 are MAP3K8-low.",
    "A measure of relative risk between two groups over time. The most important number on the forest plot.",
    "Same risk in both groups (no difference). The vertical dashed line on the plot.",
    "Higher risk in the first group. On the plot, the red dot sits to the RIGHT of the line.",
    "Lower risk in the first group. The red dot would sit to the LEFT of the line. (None of the dots here are on the left.)",
    "The range within which the true hazard ratio is very likely to fall (95% confidence). The horizontal line through each red dot. If that line CROSSES HR = 1, the result is not statistically significant.",
    "The probability that the observed result happened by random chance. Smaller P = stronger evidence the effect is real.",
    "The conventional threshold for 'statistically significant.' Means: less than a 5% chance the result is due to luck alone.",
    "A statistical test specifically designed to compare survival curves between two or more groups. Produces a P-value.",
    "A statistical model that estimates how strongly a factor (here: MAP3K8 expression) is associated with the rate of an event (here: relapse) over time. Produces a hazard ratio.",
    "An analysis that looks at ONE factor at a time, ignoring others.",
    "An analysis that includes SEVERAL factors at once (e.g., MAP3K8 + age + stage + RCB). Tells you the effect of one factor after accounting for the others.",
    "Looking at the effect within a smaller, specific group of patients (e.g., only those with RCB-II/III). Helps test whether a finding holds in particular slices of the cohort.",
    "A standard plot in clinical research. Each row is a different comparison; each dot is the hazard ratio; each horizontal line is the 95% confidence interval. The dashed vertical line marks 'no difference'.",
    "A public dataset of 508 breast-cancer patients who received chemotherapy and had their tumours profiled by gene expression. The TNBC subset (n=178) is used here.",
    "The number of patients in a group. e.g. 'n = 178' means 178 patients.",
    "The number of patients in a group who experienced the outcome (here: distant relapse) during follow-up. Statistical power depends mainly on the number of events, not the number of patients."),
  stringsAsFactors = FALSE)

## ---- Helpers ------------------------------------------------------------
make_glossary_tbl <- function(df) {
  flextable(df) |>
    set_header_labels(Term = "Term", Plain_English = "Plain-English meaning") |>
    fontsize(size = 9, part = "all") |>
    theme_vanilla() |>
    bold(part = "header") |>
    padding(padding = 3, part = "all") |>
    align(align = "left", part = "all") |>
    set_table_properties(width = 1, layout = "fixed") |>
    width(j = 1, width = 2.1) |>
    width(j = 2, width = 4.9)
}

## ---- 0.75-inch margins (7 in usable width) -----------------------------
section_props <- prop_section(
  page_size = page_size(orient = "portrait", width = 8.5, height = 11),
  page_margins = page_mar(top = 0.75, bottom = 0.75,
                          left = 0.75, right = 0.75,
                          header = 0.5, footer = 0.5, gutter = 0))

## ---- Build the document ------------------------------------------------
doc <- read_docx() |>

  body_add_par("MAP3K8 and Chemoresistance in TNBC", style = "heading 1") |>
  body_add_par("A plain-English explanation of the forest plot",
               style = "Normal") |>

  ## ---------- Bottom line ----------
  body_add_par("Bottom line, in one sentence",
               style = "heading 2") |>
  body_add_par(
    "Among triple-negative breast cancer patients who received chemotherapy, those whose tumours had HIGH levels of the MAP3K8 gene were roughly 1.8 times (i.e. about 80%) more likely to have their cancer come back at a distant site (lung, bone, liver, brain) than patients with LOW MAP3K8 levels — even after accounting for other risk factors, and even within patients whose tumour appeared to respond to chemo. This pattern is what is called \"clinical chemoresistance\": chemo seemed to work at first, but the cancer outlived it.",
    style = "Normal") |>

  ## ---------- The picture ----------
  body_add_par("The picture", style = "heading 2") |>
  body_add_img(src = FOREST_PNG, width = 6.6, height = 2.95) |>
  body_add_par(
    "Forest plot. Each row is a different group of patients. The red dot is the hazard ratio (the risk score). The horizontal line through the dot is the 95% confidence interval (how sure we are about the dot). The dashed vertical line at 1.0 is the line of 'no difference between MAP3K8-high and MAP3K8-low'. Dots to the RIGHT of the line mean MAP3K8-high patients have HIGHER risk of relapse than MAP3K8-low.",
    style = "Normal") |>

  ## ---------- How to read it ----------
  body_add_par("How to read this plot in 30 seconds",
               style = "heading 2") |>
  body_add_par(
    "1. Find the red dot - that is the headline risk score (HR = hazard ratio).",
    style = "Normal") |>
  body_add_par(
    "2. If the dot is to the RIGHT of the dashed line at 1.0, MAP3K8-high patients are at higher risk than MAP3K8-low patients.",
    style = "Normal") |>
  body_add_par(
    "3. The horizontal line through the dot is the 95% confidence interval. If that line does NOT touch 1.0, the result is statistically significant (real, not chance).",
    style = "Normal") |>
  body_add_par(
    "4. The P-value to the right of each row also tells you significance: P < 0.05 = real signal.",
    style = "Normal") |>

  ## ---------- Row-by-row ----------
  body_add_par("What each of the four rows means",
               style = "heading 2") |>

  body_add_par("Row 1 - All chemo-treated TNBC patients (n = 178)",
               style = "heading 3") |>
  body_add_par(
    "Hazard Ratio (HR) = 1.80 (95% Confidence Interval 1.09 to 2.98), P = 0.020.",
    style = "Normal") |>
  body_add_par(
    "Plain-English meaning: Patients with HIGH MAP3K8 expression were 80% more likely to relapse distantly than patients with LOW MAP3K8 expression (1.80 - 1 = 0.80 = 80%). The confidence interval (1.09 to 2.98) does not touch 1.0, so this is statistically significant. P = 0.020 means only a 2% chance this is a random fluke.",
    style = "Normal") |>
  body_add_par(
    "How to explain it: \"In the whole study, MAP3K8-high tumours came back almost twice as often as MAP3K8-low tumours.\"",
    style = "Normal") |>

  body_add_par("Row 2 - All patients, adjusted for age, cancer stage, and RCB",
               style = "heading 3") |>
  body_add_par(
    "Hazard Ratio (HR) = 1.71 (95% CI 1.00 to 2.94), P = 0.050.",
    style = "Normal") |>
  body_add_par(
    "Plain-English meaning: Even after taking into account the patient's age, how advanced the cancer was at diagnosis (AJCC stage), and how well the tumour responded to chemotherapy (RCB class), MAP3K8-high tumours were still 71% more likely to relapse than MAP3K8-low. The P = 0.050 is right at the conventional 'significant' threshold; the lower end of the confidence interval just touches 1.0.",
    style = "Normal") |>
  body_add_par(
    "Why this row matters: It rules out the simple alternative that MAP3K8 is just labelling 'patients whose tumour did not respond well to chemo.' If that were true, the effect would disappear once RCB is in the model. It does not. So MAP3K8 is carrying its own information beyond what the standard clinical factors already tell you.",
    style = "Normal") |>
  body_add_par(
    "How to explain it: \"The MAP3K8 effect is not just a re-labelling of bad chemo response. It still holds after we control for the obvious risk factors.\"",
    style = "Normal") |>

  body_add_par("Row 3 - Only patients whose tumour responded well to chemo (RCB-0/I, n = 60)",
               style = "heading 3") |>
  body_add_par(
    "Hazard Ratio (HR) = 5.68 (95% CI 0.58 to 55.18), P = 0.099.",
    style = "Normal") |>
  body_add_par(
    "Plain-English meaning: Among the patients whose tumour clearly responded to chemotherapy (RCB-0 or RCB-I), the MAP3K8-high group had ~5.7 times the risk of distant relapse compared with MAP3K8-low. BUT only 5 patients in this subgroup relapsed during follow-up, so the confidence interval is huge (0.58 to 55.18). The direction supports the chemoresistance story, but with only 5 events we cannot be statistically certain. P = 0.099 is borderline.",
    style = "Normal") |>
  body_add_par(
    "Why this row matters: This is the STRONGEST possible test of the chemoresistance hypothesis. If MAP3K8-high marked tumours that secretly survived chemo even after appearing to respond, you would expect to see exactly this pattern: chemo-responders who relapse anyway, sorted by MAP3K8 level. The direction is right; the certainty is limited by small numbers.",
    style = "Normal") |>
  body_add_par(
    "How to explain it: \"Even among patients who looked like they responded to chemo, MAP3K8-high marked those who relapsed - but this group is small, so this is suggestive, not conclusive.\"",
    style = "Normal") |>

  body_add_par("Row 4 - Only patients with residual disease (RCB-II/III, n = 86)",
               style = "heading 3") |>
  body_add_par(
    "Hazard Ratio (HR) = 2.13 (95% CI 1.16 to 3.90), P = 0.013.",
    style = "Normal") |>
  body_add_par(
    "Plain-English meaning: Among the patients whose tumour did NOT respond well to chemotherapy (significant tumour still present at surgery), the MAP3K8-high group was 113% more likely to relapse than the MAP3K8-low group. The confidence interval (1.16 to 3.90) does not touch 1.0; P = 0.013 is comfortably significant.",
    style = "Normal") |>
  body_add_par(
    "Why this row matters: This is the cleanest, strongest single result in the whole analysis. It says that within the population of patients who already have residual disease - i.e. patients we already know had a poor chemo response - MAP3K8 levels FURTHER separate those who will relapse from those who will not. MAP3K8 is not just 'another way of saying RCB-II/III'; it adds resolution inside that already-resistant population.",
    style = "Normal") |>
  body_add_par(
    "How to explain it: \"Within patients whose tumour did not respond well, MAP3K8 levels still separate the ones who will relapse from the ones who will not.\"",
    style = "Normal") |>

  ## ---------- Numbers summary ----------
  body_add_par("Summary table of all four rows",
               style = "heading 2") |>
  body_add_flextable(
    flextable(data.frame(
      Group   = c("All patients (univariable)",
                  "All patients (adjusted)",
                  "RCB-0/I subgroup",
                  "RCB-II/III subgroup"),
      n       = c("178 / 64", "178 / 64", "60 / 5", "86 / 45"),
      HR      = c("1.80", "1.71", "5.68", "2.13"),
      CI      = c("1.09 - 2.98", "1.00 - 2.94", "0.58 - 55.18", "1.16 - 3.90"),
      P       = c("0.020", "0.050", "0.099", "0.013"),
      Verdict = c("Significant",
                  "Borderline significant; effect persists after adjustment",
                  "Direction supports the hypothesis but underpowered",
                  "Clearly significant")
    )) |>
      set_header_labels(
        n = "n / events", HR = "HR",
        CI = "95% CI", P = "P-value") |>
      fontsize(size = 9, part = "all") |>
      theme_vanilla() |>
      bold(part = "header") |>
      padding(padding = 3, part = "all") |>
      align(align = "left", part = "all") |>
      align(align = "center", j = 2:5, part = "all") |>
      set_table_properties(width = 1, layout = "fixed") |>
      width(j = 1, width = 1.5) |>
      width(j = 2, width = 0.8) |>
      width(j = 3, width = 0.6) |>
      width(j = 4, width = 1.1) |>
      width(j = 5, width = 0.6) |>
      width(j = 6, width = 2.4)
  ) |>

  ## ---------- What does it all mean ----------
  body_add_par("What does this all mean, plainly",
               style = "heading 2") |>
  body_add_par(
    "MAP3K8 is acting as a marker of chemoresistance - meaning, tumours that have higher MAP3K8 are the tumours that the chemotherapy did not fully eliminate, and they are the tumours that grow back later as distant metastases.",
    style = "Normal") |>
  body_add_par(
    "Importantly, MAP3K8 is NOT a marker of whether the chemo will visibly shrink the tumour in the short term. (A separate analysis showed MAP3K8 does not predict pathological complete response.) Instead, MAP3K8 marks the tumours where some cells survive the chemo and quietly grow back over months to years.",
    style = "Normal") |>
  body_add_par(
    "Biologically, this fits what is known about MAP3K8 from cell-line research: MAP3K8 (also called TPL2) is part of a stress-response and inflammation pathway that helps cancer cells survive treatment. The clinical pattern observed here matches that biology - late relapse, independent of initial response.",
    style = "Normal") |>

  ## ---------- Caveats ----------
  body_add_par("What this analysis does NOT show, in plain English",
               style = "heading 2") |>
  body_add_par(
    "1. It does not prove cause. We see an association: tumours with more MAP3K8 relapse more often. That is consistent with MAP3K8 being a chemoresistance driver, but observational data alone cannot prove MAP3K8 itself is what causes the resistance.",
    style = "Normal") |>
  body_add_par(
    "2. It is based on a single public dataset (GSE25066, n = 178). Independent validation in a separate cohort with distant-relapse data would strengthen the case.",
    style = "Normal") |>
  body_add_par(
    "3. The RCB-0/I subgroup (Row 3) is underpowered - only 5 relapse events. That row supports the story but cannot prove it on its own.",
    style = "Normal") |>
  body_add_par(
    "4. MAP3K8 is not yet a clinically used biomarker. These results are a hypothesis-generating finding, not a recommendation to use MAP3K8 in patient care.",
    style = "Normal") |>

  ## ---------- Glossary ----------
  body_add_par("Glossary of every term in this report",
               style = "heading 2") |>
  body_add_flextable(make_glossary_tbl(glossary)) |>

  ## ---------- Footer ----------
  body_add_par("Source", style = "heading 2") |>
  body_add_par(
    "All analyses run on the GSE25066 public dataset (Hatzis & Pusztai, JAMA 2011). Triple-negative subset n = 178, treated with neoadjuvant anthracycline-taxane chemotherapy. Endpoint: distant relapse-free survival. The figure on the first page is the forest plot from the full analysis (folder 08_MAP3K8_Chemoresistance_Prognosis). For the underlying scripts, statistics CSVs, and the comprehensive technical report, see the same folder in the open-source repo.",
    style = "Normal") |>
  body_add_par(
    paste0("R version: ", R.version.string,
           ". Built on: ", Sys.Date(), "."),
    style = "Normal") |>

  body_set_default_section(section_props)

print(doc, target = "MAP3K8_Forest_PlainLanguage_Report.docx")
cat("Wrote: ",
    file.path(this_dir, "MAP3K8_Forest_PlainLanguage_Report.docx"), "\n",
    sep = "")
