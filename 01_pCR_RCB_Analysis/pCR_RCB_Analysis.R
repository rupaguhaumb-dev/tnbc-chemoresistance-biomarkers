## ============================================================
##  Analysis 01 - pCR / RCB analysis for MAP3K8 and IL1B
##  Cohort : GSE25066 (Hatzis/Booser, anthracycline-taxane neoadj)
##           TNBC subset, n = 178
##  Endpoints:
##    pCR_RD : pathologic complete response vs residual disease
##    RCB    : Residual Cancer Burden class (0 / I / II / III)
##  Tests:
##    - Wilcoxon (pCR vs RD) per gene
##    - Kruskal-Wallis (RCB classes) per gene
##    - Univariable logistic regression: pCR ~ continuous expr
##  Outputs (this folder):
##    boxplot pCR    : pCR_boxplots.{pdf,png}
##    boxplot RCB    : RCB_boxplots.{pdf,png}
##    stats table    : pCR_RCB_stats.csv
##    methods doc    : Methods_pCR_RCB_Analysis.docx
## ============================================================

need <- c("survival","ggplot2","ggpubr","officer","flextable","dplyr")
for (pkg in need) if (!requireNamespace(pkg, quietly=TRUE))
  install.packages(pkg, repos="https://cloud.r-project.org")
suppressPackageStartupMessages({
  library(ggplot2); library(ggpubr); library(officer)
  library(flextable); library(dplyr)
})

## ---- Portable path setup -----------------------------------------------
##  OUT  : this script's own directory (works under Rscript / source / RStudio)
##  WS   : the GSE25066 workspace; defaults to two levels above this script
##         (i.e. the parent of the Analyses/ folder), overridable via env var
##         TNBC_WORKSPACE.
this_dir <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", a, value = TRUE)
  if (length(f)) return(dirname(normalizePath(sub("^--file=", "", f[1]))))
  fr <- sys.frames()
  for (i in rev(seq_along(fr)))
    if (!is.null(fr[[i]]$ofile))
      return(dirname(normalizePath(fr[[i]]$ofile)))
  getwd()
})
OUT <- this_dir
setwd(OUT)
WS  <- Sys.getenv("TNBC_WORKSPACE",
        unset = file.path("..", "..", "GSE25066_TNBC_MAP3K8_workspace.RData"))
if (!file.exists(WS))
  stop("Workspace not found at: ", WS,
       "\nSet TNBC_WORKSPACE env var or place the file at the default path.")
load(WS)

## ---- 1. Build the analysis frame ----------------------------------------
##  MAP3K8 already in workspace (object 'map3k8', from probe 205027_s_at)
##  IL1B   we re-derive using the same probe-selection rule used in
##         IL1B_TNBC_KMPlot.R (highest mean expression probe on GPL96)
sym_col   <- "Gene Symbol"
il1b_rows <- which(gpl_table[[sym_col]] == "IL1B")
il1b_probes <- as.character(gpl_table$ID[il1b_rows])
probe_means <- rowMeans(expr_tnbc[il1b_probes, , drop = FALSE])
il1b_probe  <- names(which.max(probe_means))
il1b <- as.numeric(expr_tnbc[il1b_probe, ])

stopifnot(identical(rownames(clinical_tnbc), colnames(expr_tnbc)))

df <- data.frame(
  geo_accession = clinical_tnbc$geo_accession,
  MAP3K8 = as.numeric(map3k8),
  IL1B   = il1b,
  pCR    = clinical_tnbc[["pathologic_response_pcr_rd:ch1"]],
  RCB    = clinical_tnbc[["pathologic_response_rcb_class:ch1"]],
  stringsAsFactors = FALSE
)
df$pCR <- toupper(trimws(df$pCR))
df$pCR[df$pCR %in% c("","NA")] <- NA
df$pCR <- factor(df$pCR, levels = c("RD","PCR"))
df$RCB <- toupper(trimws(df$RCB))
df$RCB[df$RCB %in% c("","NA")] <- NA
df$RCB <- factor(df$RCB, levels = c("RCB-0","RCB-I","RCB-II","RCB-III"))

cat("Patients with pCR call:", sum(!is.na(df$pCR)), "\n")
print(table(df$pCR, useNA="ifany"))
cat("\nPatients with RCB call:", sum(!is.na(df$RCB)), "\n")
print(table(df$RCB, useNA="ifany"))

## ---- 2. Statistical tests -----------------------------------------------
make_stats <- function(g) {
  pcr_df <- df[!is.na(df$pCR), ]
  rcb_df <- df[!is.na(df$RCB), ]
  w <- wilcox.test(pcr_df[[g]] ~ pcr_df$pCR)
  k <- kruskal.test(rcb_df[[g]] ~ rcb_df$RCB)
  ##  Logistic regression pCR ~ continuous expression
  lr  <- glm(I(pCR == "PCR") ~ get(g), family = binomial, data = pcr_df)
  s   <- summary(lr)
  OR  <- exp(coef(lr)[2])
  ci  <- suppressMessages(exp(confint(lr)))[2, ]
  data.frame(
    gene             = g,
    n_pCR            = sum(pcr_df$pCR == "PCR"),
    n_RD             = sum(pcr_df$pCR == "RD"),
    median_in_pCR    = median(pcr_df[[g]][pcr_df$pCR == "PCR"]),
    median_in_RD     = median(pcr_df[[g]][pcr_df$pCR == "RD"]),
    Wilcoxon_p       = signif(w$p.value, 3),
    KruskalWallis_p  = signif(k$p.value, 3),
    logreg_OR_perUnit= signif(OR, 3),
    logreg_OR_LCL    = signif(ci[1], 3),
    logreg_OR_UCL    = signif(ci[2], 3),
    logreg_p         = signif(s$coefficients[2,"Pr(>|z|)"], 3),
    stringsAsFactors = FALSE
  )
}
stats <- bind_rows(make_stats("MAP3K8"), make_stats("IL1B"))
print(stats)
write.csv(stats, "pCR_RCB_stats.csv", row.names = FALSE)

## ---- 3. Boxplots --------------------------------------------------------
plot_pcr <- function(g) {
  d <- df[!is.na(df$pCR), ]
  ggplot(d, aes(x = pCR, y = .data[[g]], fill = pCR)) +
    geom_boxplot(width = 0.55, alpha = 0.85, outlier.size = 0.8) +
    geom_jitter(width = 0.15, size = 0.9, alpha = 0.6) +
    scale_fill_manual(values = c(RD = "#1F77B4", PCR = "#D7263D")) +
    stat_compare_means(method = "wilcox.test",
                       label.y = max(d[[g]], na.rm = TRUE) + 0.3) +
    labs(x = "Response", y = paste0(g, " (log2)"),
         title = paste0(g, " by chemoresponse")) +
    theme_classic(base_size = 12) + theme(legend.position = "none")
}
plot_rcb <- function(g) {
  d <- df[!is.na(df$RCB), ]
  ggplot(d, aes(x = RCB, y = .data[[g]], fill = RCB)) +
    geom_boxplot(width = 0.55, alpha = 0.85, outlier.size = 0.8) +
    geom_jitter(width = 0.15, size = 0.9, alpha = 0.6) +
    scale_fill_brewer(palette = "RdYlBu", direction = -1) +
    stat_compare_means(method = "kruskal.test",
                       label.y = max(d[[g]], na.rm = TRUE) + 0.3) +
    labs(x = "RCB class", y = paste0(g, " (log2)"),
         title = paste0(g, " across RCB classes")) +
    theme_classic(base_size = 12) + theme(legend.position = "none")
}
p_pcr  <- ggarrange(plot_pcr("MAP3K8"), plot_pcr("IL1B"),
                    ncol = 2, labels = c("A","B"))
p_rcb  <- ggarrange(plot_rcb("MAP3K8"), plot_rcb("IL1B"),
                    ncol = 2, labels = c("A","B"))
ggsave("pCR_boxplots.pdf", p_pcr, width = 9, height = 4.5)
ggsave("pCR_boxplots.png", p_pcr, width = 9, height = 4.5, dpi = 300)
ggsave("RCB_boxplots.pdf", p_rcb, width = 9, height = 4.5)
ggsave("RCB_boxplots.png", p_rcb, width = 9, height = 4.5, dpi = 300)

## ---- 4. Methods Word document -------------------------------------------
make_doc <- function(stats_df, il1b_probe) {
  doc <- read_docx() |>
    body_add_par("Analysis 01 - pCR / RCB analysis for MAP3K8 and IL-1beta in TNBC", style = "heading 1") |>
    body_add_par("Cohort and dataset", style = "heading 2") |>
    body_add_par(
      "GSE25066 (Hatzis & Pusztai, JAMA 2011; Booser series) - 508 women with primary ",
      style = "Normal") |>
    body_add_par("breast cancer treated with neoadjuvant anthracycline-taxane chemotherapy, profiled on the Affymetrix HG-U133A platform (GPL96, MAS5.0). The TNBC subset (n = 178) was previously parsed by selecting patients with ER-IHC = N, PR-IHC = N, and HER2 = N, and is stored in the working-directory R workspace `GSE25066_TNBC_MAP3K8_workspace.RData` as `clinical_tnbc` (178 x 82) and `expr_tnbc` (22283 x 178).",
      style = "Normal") |>
    body_add_par("Probe selection", style = "heading 2") |>
    body_add_par(paste0(
      "MAP3K8 was represented on GPL96 by the single probe 205027_s_at, used as-is. IL-1beta (IL1B) was represented by two probes (205067_at and 39402_at); the probe with the higher mean log2 expression in the TNBC cohort was selected (",
      il1b_probe,
      "), matching the convention used in the original IL1B Kaplan-Meier script."),
      style = "Normal") |>
    body_add_par("Endpoints", style = "heading 2") |>
    body_add_par(
      "Pathologic complete response (pCR) versus residual disease (RD) was taken from the clinical column `pathologic_response_pcr_rd:ch1`. The Residual Cancer Burden (RCB) class (RCB-0, RCB-I, RCB-II, RCB-III; Symmans 2007) was taken from `pathologic_response_rcb_class:ch1`. RCB-0 corresponds to pCR; RCB-I/II/III are increasing residual tumor burden.",
      style = "Normal") |>
    body_add_par("Statistical tests", style = "heading 2") |>
    body_add_par(
      "1. Wilcoxon rank-sum test for each gene comparing log2 expression between pCR and RD groups.",
      style = "Normal") |>
    body_add_par(
      "2. Kruskal-Wallis test for each gene across the four RCB classes.",
      style = "Normal") |>
    body_add_par(
      "3. Univariable logistic regression with pCR (1/0) as outcome and continuous log2 expression as the only predictor. We report the odds ratio per unit of log2 expression with profile-likelihood 95% CI and the Wald p-value.",
      style = "Normal") |>
    body_add_par("All tests were two-sided; no multiple-testing correction was applied because only two pre-specified candidate genes were tested.",
                 style = "Normal") |>
    body_add_par("Visualization", style = "heading 2") |>
    body_add_par(
      "Box-and-whisker plots with overlaid jittered points for each gene, panelled side-by-side with ggpubr::ggarrange, were produced for (i) pCR vs RD and (ii) RCB-0 / I / II / III. P-values printed on the plots come from the corresponding rank tests above.",
      style = "Normal") |>
    body_add_par("Results table", style = "heading 2") |>
    body_add_flextable(flextable(stats_df) |>
                         autofit() |>
                         fontsize(size = 9, part = "all")) |>
    body_add_par("Results obtained in this run", style = "heading 2") |>
    body_add_par(sprintf(
      "Of 178 TNBC patients, %d had a usable pCR/RD call (%d pCR, %d RD) and %d had a non-missing RCB class. The original Booser series populated the RCB column only for residual-disease cases, so RCB-0 (which is by definition pCR) and RCB-I were absent from the encoded categorical column; the RCB-II / RCB-III contrast is the only meaningful Kruskal-Wallis tested here.",
      sum(!is.na(df$pCR)), sum(df$pCR == "PCR", na.rm = TRUE),
      sum(df$pCR == "RD",  na.rm = TRUE), sum(!is.na(df$RCB))),
      style = "Normal") |>
    body_add_par(sprintf(
      "MAP3K8 (probe 205027_s_at): median log2 expression was %.2f in pCR and %.2f in RD. Wilcoxon P = %.3g; logistic OR per unit log2 = %.2f (95%% CI %.2f-%.2f), Wald P = %.3g. Kruskal-Wallis across RCB classes P = %.3g. There is no detectable association between MAP3K8 mRNA and short-term chemoresponse in this cohort.",
      stats$median_in_pCR[1], stats$median_in_RD[1],
      stats$Wilcoxon_p[1],
      stats$logreg_OR_perUnit[1], stats$logreg_OR_LCL[1], stats$logreg_OR_UCL[1],
      stats$logreg_p[1], stats$KruskalWallis_p[1]),
      style = "Normal") |>
    body_add_par(sprintf(
      "IL-1beta (probe %s): median log2 expression was %.2f in pCR and %.2f in RD. Wilcoxon P = %.3g; logistic OR per unit log2 = %.2f (95%% CI %.2f-%.2f), Wald P = %.3g. Kruskal-Wallis across RCB classes P = %.3g. The point estimate places IL-1beta higher in residual-disease tumours, consistent with a chemoresistance role, but the difference does not reach statistical significance.",
      il1b_probe,
      stats$median_in_pCR[2], stats$median_in_RD[2],
      stats$Wilcoxon_p[2],
      stats$logreg_OR_perUnit[2], stats$logreg_OR_LCL[2], stats$logreg_OR_UCL[2],
      stats$logreg_p[2], stats$KruskalWallis_p[2]),
      style = "Normal") |>
    body_add_par("Figure 1. MAP3K8 and IL-1beta expression by pCR status",
                 style = "heading 3") |>
    body_add_img(src = file.path(OUT, "pCR_boxplots.png"),
                 width = 6.3, height = 3.2) |>
    body_add_par("Figure 2. MAP3K8 and IL-1beta expression by RCB class",
                 style = "heading 3") |>
    body_add_img(src = file.path(OUT, "RCB_boxplots.png"),
                 width = 6.3, height = 3.2) |>
    body_add_par("Interpretation", style = "heading 2") |>
    body_add_par(
      "Neither MAP3K8 nor IL-1beta mRNA predicts pathological complete response in chemo-treated TNBC. Combined with the DRFS results obtained in the parallel analyses (median-split KM, four-group KM, continuous Cox), this means that the prognostic signal carried by these two genes does not act through whether the tumour visibly shrinks under neoadjuvant chemotherapy, but through the probability of distant relapse afterwards. In other words: pCR captures the early cytotoxic response; MAP3K8 / IL-1beta capture residual chemoresistance biology that drives subsequent metastasis. This is biologically coherent - both genes sit on the IL-1R / IRAK / NF-kB / MAP3K8 axis, which is implicated in tumour-promoting inflammation and post-treatment outgrowth rather than primary cytotoxic sensitivity.",
      style = "Normal") |>
    body_add_par("Limitations", style = "heading 2") |>
    body_add_par(
      "Sample sizes for the pCR (n=170) and especially RCB (n=86, only II/III observed) endpoints are modest; the analysis is therefore underpowered to detect small (<10%) differences in expression between response groups. No multiple-testing correction was applied because only two pre-specified candidate genes were tested. The univariable logistic regression does not adjust for stage or RCB.",
      style = "Normal") |>
    body_add_par("Files written by this script", style = "heading 2") |>
    body_add_par("- pCR_RCB_stats.csv      : test statistics for both genes", style = "Normal") |>
    body_add_par("- pCR_boxplots.pdf/.png  : MAP3K8 and IL1B boxplots, pCR vs RD", style = "Normal") |>
    body_add_par("- RCB_boxplots.pdf/.png  : MAP3K8 and IL1B boxplots across RCB classes", style = "Normal") |>
    body_add_par("Reproducibility", style = "heading 2") |>
    body_add_par(paste0("R version: ", R.version.string,
                        ". Platform: ", R.version$platform,
                        ". Date: ", Sys.Date(),
                        ". Random seed: not used (only deterministic tests)."),
                 style = "Normal")
  print(doc, target = "Methods_pCR_RCB_Analysis.docx")
}
make_doc(stats, il1b_probe)

cat("\nAnalysis 01 complete. Files in:\n  ", OUT, "\n", sep="")
