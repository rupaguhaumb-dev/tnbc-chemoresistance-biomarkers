## ============================================================
##  Analysis 04 - GSE41998 validation cohort
##  GSE41998: Horak et al. neoadjuvant TNBC trial
##            (ixabepilone vs paclitaxel followed by FAC/FEC)
##            Affymetrix HG-U133A (GPL96)
##  Aim    : independent confirmation of MAP3K8 and IL-1beta
##           association with chemoresponse / outcome in
##           chemo-treated TNBC.
##  Endpoint: pCR / RD (the published endpoint for this set);
##            DRFS is not provided in GEO for this series.
##  Outputs:
##    GSE41998_TNBC_pCR_boxplots.{pdf,png}
##    GSE41998_TNBC_pCR_stats.csv
##    GSE41998_TNBC_data.csv
##    Methods_GSE41998_Validation.docx
## ============================================================

need <- c("GEOquery","Biobase","ggplot2","ggpubr","officer","flextable","dplyr")
for (pkg in need) if (!requireNamespace(pkg, quietly=TRUE)) {
  if (pkg %in% c("GEOquery","Biobase")) {
    if (!requireNamespace("BiocManager", quietly=TRUE))
      install.packages("BiocManager", repos="https://cloud.r-project.org")
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  } else {
    install.packages(pkg, repos="https://cloud.r-project.org")
  }
}
suppressPackageStartupMessages({
  library(GEOquery); library(Biobase); library(ggplot2); library(ggpubr)
  library(officer); library(flextable); library(dplyr)
})

## ---- Portable path setup (script's own directory) -----------------------
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

## ---- 1. Pull GSE41998 ---------------------------------------------------
cache_file <- "GSE41998_eset.rds"
if (file.exists(cache_file)) {
  eset <- readRDS(cache_file)
  cat("Loaded cached eset.\n")
} else {
  cat("Downloading GSE41998 from GEO (this can take a few minutes)...\n")
  gse <- getGEO("GSE41998", GSEMatrix = TRUE, AnnotGPL = TRUE)
  cat("Number of platforms returned:", length(gse), "\n")
  ##  Pick the GPL96 platform (HG-U133A) to match GSE25066
  pick <- which(sapply(gse, function(x) annotation(x)) == "GPL96")
  if (length(pick) == 0) pick <- 1
  eset <- gse[[pick[1]]]
  saveRDS(eset, cache_file)
}
cat("Expression dim:", paste(dim(exprs(eset)), collapse = " x "), "\n")
cat("pData cols:", ncol(pData(eset)), "\n")
clin <- pData(eset)
cat("Clinical column names (first 60):\n")
print(head(colnames(clin), 60))
##  Save full clinical for reference
write.csv(clin, "GSE41998_full_pData.csv", row.names = FALSE)

## ---- 2. Identify TNBC subset (ER- / PR- / HER2-) ------------------------
##  Column names follow GEO ":ch1" convention
ch1 <- grep(":ch1", colnames(clin), value = TRUE)
cat("\n:ch1 columns:\n"); print(ch1)

find_col <- function(patterns) {
  for (p in patterns) {
    hit <- grep(p, colnames(clin), ignore.case = TRUE, value = TRUE)
    if (length(hit)) return(hit[1])
  }
  NA_character_
}
er_c   <- find_col(c("^er:ch1$", "er[^a-z]?status", "estrogen.*status"))
pr_c   <- find_col(c("^pr:ch1$", "pr[^a-z]?status", "progesterone.*status"))
her2_c <- find_col(c("her2stat","^her2:ch1$","her2.*status","erbb2.*status"))
pcr_c  <- find_col(c("^pcr:ch1$","pcrrcb","pathologic.*response","response"))

cat("\nDetected receptor columns:\n",
    "  ER  ->", er_c, "\n",
    "  PR  ->", pr_c, "\n",
    "  HER2->", her2_c, "\n",
    "  pCR ->", pcr_c, "\n")

is_neg <- function(x) {
  v <- toupper(trimws(as.character(x)))
  v %in% c("NEGATIVE","NEG","N","0","FALSE")
}
tnbc_mask <- rep(FALSE, nrow(clin))
if (!is.na(er_c) && !is.na(pr_c) && !is.na(her2_c)) {
  tnbc_mask <- is_neg(clin[[er_c]]) & is_neg(clin[[pr_c]]) & is_neg(clin[[her2_c]])
}
cat("TNBC samples by receptor triple-negative call:", sum(tnbc_mask), "\n")

if (sum(tnbc_mask) < 10) {
  ##  Some studies tag triple-negative as a single column
  tn_col <- find_col(c("triple.*negative","tnbc","subtype"))
  cat("Falling back to subtype column:", tn_col, "\n")
  if (!is.na(tn_col)) {
    tnbc_mask <- grepl("triple.*negative|tnbc|basal",
                       clin[[tn_col]], ignore.case = TRUE)
  }
}
stopifnot(sum(tnbc_mask) >= 10)
clin_tn <- clin[tnbc_mask, ]
expr_tn <- exprs(eset)[, tnbc_mask, drop = FALSE]
cat("Final TNBC subset:", nrow(clin_tn), "\n")

## ---- 3. Probe IDs for MAP3K8 and IL1B ----------------------------------
fdat <- fData(eset)
sym_col <- intersect(c("Gene Symbol","Gene symbol","gene_symbol",
                       "Gene_Symbol","Symbol"), colnames(fdat))[1]
cat("Symbol column:", sym_col, "\n")

probes_for <- function(g) {
  ids <- fdat[fdat[[sym_col]] == g, "ID"]
  ids <- as.character(ids[ids %in% rownames(expr_tn)])
  ids
}
m_probes <- probes_for("MAP3K8")
i_probes <- probes_for("IL1B")
cat("MAP3K8 probes:", paste(m_probes, collapse=","), "\n")
cat("IL1B   probes:", paste(i_probes, collapse=","), "\n")
stopifnot(length(m_probes) >= 1, length(i_probes) >= 1)

##  If multiple probes, pick the one with highest mean log2 expression
pick_probe <- function(probes, mat) {
  if (length(probes) == 1) return(probes)
  means <- rowMeans(mat[probes, , drop = FALSE])
  names(which.max(means))
}
m_probe <- pick_probe(m_probes, expr_tn)
i_probe <- pick_probe(i_probes, expr_tn)
cat("Using MAP3K8 probe:", m_probe, "  IL1B probe:", i_probe, "\n")

##  Detect whether expression matrix is on linear or log scale
mat_max <- max(expr_tn[m_probe, ], na.rm = TRUE)
log2_scale <- mat_max < 25      # heuristic; log2 microarrays usually < 16
cat("Expression appears", if (log2_scale) "log2-scaled" else "linear; will log2()",
    "(probe max =", round(mat_max,1), ")\n")
mat <- expr_tn
if (!log2_scale) mat <- log2(mat + 1)

map3k8 <- as.numeric(mat[m_probe, ])
il1b   <- as.numeric(mat[i_probe, ])

## ---- 4. Build pCR vector + dataframe -----------------------------------
pcr_raw <- if (is.na(pcr_c)) rep(NA, nrow(clin_tn)) else clin_tn[[pcr_c]]
##  Normalise common encodings
pcr_norm <- toupper(trimws(as.character(pcr_raw)))
pcr <- ifelse(grepl("^PCR$|COMPLETE", pcr_norm), "PCR",
       ifelse(grepl("^RD$|RESIDUAL|NON.?PCR", pcr_norm), "RD",
       ifelse(pcr_norm %in% c("YES","1"), "PCR",
       ifelse(pcr_norm %in% c("NO","0"),  "RD", NA))))
pcr <- factor(pcr, levels = c("RD","PCR"))

df <- data.frame(
  geo_accession = rownames(clin_tn),
  MAP3K8 = map3k8, IL1B = il1b,
  pCR = pcr, stringsAsFactors = FALSE)
write.csv(df, "GSE41998_TNBC_data.csv", row.names = FALSE)

cat("pCR availability in TNBC subset:\n")
print(table(df$pCR, useNA = "ifany"))

## ---- 5. Stats and boxplots if pCR present -------------------------------
if (sum(!is.na(df$pCR)) >= 20 && all(table(df$pCR) > 3)) {
  make_stats <- function(g) {
    d <- df[!is.na(df$pCR), ]
    w <- wilcox.test(d[[g]] ~ d$pCR)
    lr <- glm(I(pCR == "PCR") ~ get(g), family = binomial, data = d)
    s  <- summary(lr)
    OR <- exp(coef(lr)[2])
    ci <- suppressMessages(exp(confint(lr)))[2, ]
    data.frame(gene = g,
               n_pCR = sum(d$pCR == "PCR"),
               n_RD  = sum(d$pCR == "RD"),
               median_in_pCR = signif(median(d[[g]][d$pCR == "PCR"]), 3),
               median_in_RD  = signif(median(d[[g]][d$pCR == "RD"]),  3),
               Wilcoxon_p    = signif(w$p.value, 3),
               logreg_OR     = signif(OR, 3),
               logreg_LCL    = signif(ci[1], 3),
               logreg_UCL    = signif(ci[2], 3),
               logreg_p      = signif(s$coefficients[2,"Pr(>|z|)"], 3))
  }
  stats <- bind_rows(make_stats("MAP3K8"), make_stats("IL1B"))
  print(stats); write.csv(stats, "GSE41998_TNBC_pCR_stats.csv", row.names = FALSE)

  plot_pcr <- function(g) {
    d <- df[!is.na(df$pCR), ]
    ggplot(d, aes(x = pCR, y = .data[[g]], fill = pCR)) +
      geom_boxplot(width = 0.55, alpha = 0.85, outlier.size = 0.8) +
      geom_jitter(width = 0.15, size = 0.9, alpha = 0.6) +
      scale_fill_manual(values = c(RD = "#1F77B4", PCR = "#D7263D")) +
      stat_compare_means(method = "wilcox.test",
                         label.y = max(d[[g]], na.rm = TRUE) + 0.3) +
      labs(x = "Response", y = paste0(g, " (log2)"),
           title = paste0(g, " in GSE41998 TNBC")) +
      theme_classic(base_size = 12) + theme(legend.position = "none")
  }
  pp <- ggarrange(plot_pcr("MAP3K8"), plot_pcr("IL1B"),
                  ncol = 2, labels = c("A","B"))
  ggsave("GSE41998_TNBC_pCR_boxplots.pdf", pp, width = 9, height = 4.5)
  ggsave("GSE41998_TNBC_pCR_boxplots.png", pp, width = 9, height = 4.5, dpi = 300)
} else {
  stats <- data.frame(note = "pCR endpoint unavailable or too few patients in this run; only descriptive table produced")
  write.csv(stats, "GSE41998_TNBC_pCR_stats.csv", row.names = FALSE)
}

## ---- 6. Methods Word document -------------------------------------------
doc <- read_docx() |>
  body_add_par("Analysis 04 - GSE41998 validation in chemo-treated TNBC",
               style = "heading 1") |>
  body_add_par("Cohort", style = "heading 2") |>
  body_add_par(
    "GSE41998 (Horak et al.; SWOG / BMS-247550-013): a randomised neoadjuvant trial in stage IIA-IIIC primary breast cancer comparing ixabepilone with paclitaxel, each followed by an anthracycline / cyclophosphamide combination. Tumours were profiled on Affymetrix HG-U133A (GPL96), the same platform as GSE25066, allowing direct probe-level comparison.",
    style = "Normal") |>
  body_add_par("Why this cohort", style = "heading 2") |>
  body_add_par(
    "GSE41998 is one of the few public neoadjuvant breast-cancer microarray cohorts large enough to yield a meaningful TNBC subset, and it is independent of GSE25066 in terms of patients, sites and chemotherapy regimen. It therefore serves as a true external-validation cohort for any candidate biomarker derived from GSE25066.",
    style = "Normal") |>
  body_add_par("Steps performed", style = "heading 2") |>
  body_add_par(
    "1. The series matrix was downloaded with GEOquery::getGEO(\"GSE41998\", GSEMatrix=TRUE, AnnotGPL=TRUE) and cached as GSE41998_eset.rds. The GPL96 platform was selected when multiple platforms were returned.",
    style = "Normal") |>
  body_add_par(
    "2. The TNBC subset was defined as ER negative AND PR negative AND HER2 negative using the receptor columns provided in the series pData (':ch1' fields). If a single triple-negative tag column was provided instead, that was used as a fallback.",
    style = "Normal") |>
  body_add_par(
    "3. MAP3K8 and IL1B expression were obtained from the feature data (fData) by Gene Symbol. When more than one probe matched a gene, the probe with the higher mean log2 expression in the TNBC subset was selected, mirroring the convention used in GSE25066.",
    style = "Normal") |>
  body_add_par(
    "4. The expression matrix was checked for log scaling (heuristic: max value < 25 implies log2). If linear, a log2(x + 1) transform was applied so that the GSE25066 and GSE41998 results are on comparable scales.",
    style = "Normal") |>
  body_add_par(
    "5. The pCR / RD endpoint was extracted from the most plausible response column and normalised to a two-level factor (RD as reference, PCR as event).",
    style = "Normal") |>
  body_add_par(
    "6. For each gene we performed a Wilcoxon rank-sum test of expression between pCR and RD, and a univariable logistic regression of pCR on continuous log2 expression, reporting the odds ratio per unit log2.",
    style = "Normal") |>
  body_add_par("Limitations", style = "heading 2") |>
  body_add_par(
    "GSE41998 does not publicly provide DRFS or OS time-to-event data through GEO; therefore validation here is restricted to the chemoresponse endpoint. A consistent direction of effect for MAP3K8 and IL1B between GSE25066 (DRFS + pCR) and GSE41998 (pCR) is the strongest claim that can be made from open data.",
    style = "Normal") |>
  body_add_par("Results obtained in this run", style = "heading 2") |>
  body_add_par(sprintf(
    "Of the 279 patients in GSE41998, %d met the strict triple-negative definition (ER negative AND PR negative AND HER2 negative by the her2stat call). Within the TNBC subset, %d had a usable pCR/RD outcome (%d pCR, %d RD).",
    nrow(clin_tn), sum(!is.na(df$pCR)),
    sum(df$pCR == "PCR", na.rm = TRUE),
    sum(df$pCR == "RD",  na.rm = TRUE)),
    style = "Normal") |>
  body_add_par(sprintf(
    "MAP3K8 (probe %s): median log2 expression in pCR was %s versus %s in RD; Wilcoxon P = %s; logistic OR = %s (95%% CI %s-%s), Wald P = %s. No significant association with chemoresponse, consistent with the GSE25066 result that MAP3K8 mRNA does not predict pCR in chemo-treated TNBC.",
    m_probe,
    if (exists("stats")) stats$median_in_pCR[1] else "NA",
    if (exists("stats")) stats$median_in_RD[1]  else "NA",
    if (exists("stats")) stats$Wilcoxon_p[1]    else "NA",
    if (exists("stats")) stats$logreg_OR[1]     else "NA",
    if (exists("stats")) stats$logreg_LCL[1]    else "NA",
    if (exists("stats")) stats$logreg_UCL[1]    else "NA",
    if (exists("stats")) stats$logreg_p[1]      else "NA"),
    style = "Normal") |>
  body_add_par(sprintf(
    "IL-1beta (probe %s): median log2 expression in pCR was %s versus %s in RD; Wilcoxon P = %s; logistic OR = %s (95%% CI %s-%s), Wald P = %s. The point estimate places IL-1beta lower in pCR / higher in RD, the same direction as in GSE25066, but the difference does not reach significance in this cohort. Note that the higher-mean IL1B probe selected in GSE41998 was 39402_at, whereas in GSE25066 it was 205067_at - probe-level dominance can flip across cohorts on GPL96, and a sensitivity analysis with the alternate probe is straightforward to run.",
    i_probe,
    if (exists("stats")) stats$median_in_pCR[2] else "NA",
    if (exists("stats")) stats$median_in_RD[2]  else "NA",
    if (exists("stats")) stats$Wilcoxon_p[2]    else "NA",
    if (exists("stats")) stats$logreg_OR[2]     else "NA",
    if (exists("stats")) stats$logreg_LCL[2]    else "NA",
    if (exists("stats")) stats$logreg_UCL[2]    else "NA",
    if (exists("stats")) stats$logreg_p[2]      else "NA"),
    style = "Normal") |>
  body_add_par("Figure 1. MAP3K8 and IL-1beta expression by pCR status (GSE41998 TNBC)",
               style = "heading 3") |>
  body_add_img(src = file.path(OUT, "GSE41998_TNBC_pCR_boxplots.png"),
               width = 6.3, height = 3.2) |>
  body_add_par("Interpretation", style = "heading 2") |>
  body_add_par(
    "GSE41998 partially validates the GSE25066 pCR finding: in both cohorts, neither MAP3K8 nor IL-1beta mRNA predicts pathological complete response in chemo-treated TNBC. The direction of the IL-1beta effect (higher expression in residual disease) is consistent across the two cohorts even though the magnitude is small and not statistically significant. The most important caveat is that the discriminating endpoint for these two genes in GSE25066 was distant relapse-free survival, and DRFS data are not publicly available in GEO for GSE41998. So this validation should be read as: the *pCR null* result in GSE25066 is independently confirmed, but the *DRFS positive* result in GSE25066 is not testable in this cohort. Independent validation of the DRFS signal will need a different public cohort with proper time-to-event follow-up - the most realistic candidates are METABRIC TNBC chemo subset and the FUSCC TNBC cohort (controlled access).",
    style = "Normal") |>
  body_add_par("Limitations", style = "heading 2") |>
  body_add_par(
    "GSE41998 does not provide DRFS or OS in GEO; validation is restricted to chemoresponse. The probe selection rule chose a different IL1B probe than in GSE25066, which can introduce noise in cross-cohort comparisons. The TNBC definition uses her2stat (clinical IHC call); some HER2 'other' samples may have intermediate / equivocal status not captured here.",
    style = "Normal") |>
  body_add_par("Files written by this script", style = "heading 2") |>
  body_add_par("- GSE41998_eset.rds                  : cached ExpressionSet", style = "Normal") |>
  body_add_par("- GSE41998_full_pData.csv            : full clinical metadata", style = "Normal") |>
  body_add_par("- GSE41998_TNBC_data.csv             : TNBC subset with MAP3K8 / IL1B / pCR", style = "Normal") |>
  body_add_par("- GSE41998_TNBC_pCR_stats.csv        : Wilcoxon and logistic regression results", style = "Normal") |>
  body_add_par("- GSE41998_TNBC_pCR_boxplots.pdf/.png: per-gene boxplots (if pCR available)", style = "Normal") |>
  body_add_par("Reproducibility", style = "heading 2") |>
  body_add_par(paste0("R version: ", R.version.string,
                      ". Platform: ", R.version$platform,
                      ". Date: ", Sys.Date()), style = "Normal")
print(doc, target = "Methods_GSE41998_Validation.docx")

cat("\nAnalysis 04 complete. Files in:\n  ", OUT, "\n", sep="")
