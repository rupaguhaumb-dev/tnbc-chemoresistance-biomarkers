## ============================================================
##  Analysis 07 - MAP3K8 as a chemoresistance / chemosensitivity
##                predictor in chemo-treated TNBC
##  Cohort : GSE25066 TNBC (n = 178), neoadjuvant anthracycline-
##           taxane chemotherapy.
##
##  Endpoints (all from clinical_tnbc):
##    A. Actual pCR vs RD          (pathologic_response_pcr_rd:ch1)
##    B. Actual RCB class          (pathologic_response_rcb_class:ch1)
##    C. Hatzis genomic signature  (chemosensitivity_prediction:ch1)
##                                 -> Rx Sensitive vs Rx Insensitive
##    D. DLDA30 prediction         (dlda30_prediction:ch1)
##    E. RCB-0/I prediction        (rcb_0_i_prediction:ch1)
##
##  Tests / outputs:
##    1. Distribution boxplots of MAP3K8 by each prediction class
##    2. Wilcoxon (and Kruskal-Wallis for RCB class)
##    3. Univariable logistic regression: pCR ~ MAP3K8 (continuous)
##    4. ROC curves and AUC for MAP3K8 -> pCR, vs the Hatzis
##       signature, plus a combined logistic model
##    5. Optimal cutpoint by Youden's J; sensitivity, specificity,
##       PPV, NPV at median and at optimal cutpoints
##    6. 2x2 concordance table and Cohen's kappa between MAP3K8
##       median-split and the Hatzis chemosensitivity prediction
##
##  Outputs (this folder):
##    boxplots_MAP3K8_by_class.{pdf,png}
##    ROC_MAP3K8_vs_signature.{pdf,png}
##    MAP3K8_chemoresponse_stats.csv
##    MAP3K8_classifier_metrics.csv
##    MAP3K8_vs_Hatzis_concordance.csv
## ============================================================

need <- c("survival","ggplot2","ggpubr","pROC","dplyr","cutpointr")
for (pkg in need) if (!requireNamespace(pkg, quietly=TRUE))
  install.packages(pkg, repos="https://cloud.r-project.org")
suppressPackageStartupMessages({
  library(ggplot2); library(ggpubr); library(pROC)
  library(dplyr)
  ok_cutpointr <- requireNamespace("cutpointr", quietly = TRUE)
  if (ok_cutpointr) library(cutpointr)
})

## ---- Portable path setup ------------------------------------------------
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
WS <- Sys.getenv("TNBC_WORKSPACE",
       unset = file.path("..", "..", "GSE25066_TNBC_MAP3K8_workspace.RData"))
if (!file.exists(WS))
  stop("Workspace not found at: ", WS)
load(WS)

## ---- Build analysis frame ----------------------------------------------
stopifnot(identical(rownames(clinical_tnbc), colnames(expr_tnbc)))
df <- data.frame(
  geo_accession = clinical_tnbc$geo_accession,
  MAP3K8   = as.numeric(map3k8),
  pCR      = clinical_tnbc[["pathologic_response_pcr_rd:ch1"]],
  RCB      = clinical_tnbc[["pathologic_response_rcb_class:ch1"]],
  Hatzis   = clinical_tnbc[["chemosensitivity_prediction:ch1"]],
  DLDA30   = clinical_tnbc[["dlda30_prediction:ch1"]],
  RCB0I    = clinical_tnbc[["rcb_0_i_prediction:ch1"]],
  stringsAsFactors = FALSE)

clean_factor <- function(x, lvls) {
  v <- toupper(trimws(as.character(x)))
  v[v %in% c("", "NA")] <- NA
  factor(v, levels = toupper(lvls))
}
df$pCR    <- clean_factor(df$pCR,    c("RD","PCR"))
df$RCB    <- clean_factor(df$RCB,    c("RCB-0/I","RCB-II","RCB-III"))
df$Hatzis <- clean_factor(df$Hatzis, c("RX INSENSITIVE","RX SENSITIVE"))
df$DLDA30 <- clean_factor(df$DLDA30, c("RD","PCR"))
df$RCB0I  <- clean_factor(df$RCB0I,  c("RCB-II/III","RCB-0/I"))

cat("Cohort sizes per prediction call:\n")
cat("  Actual pCR/RD       :", sum(!is.na(df$pCR)),    "\n")
cat("  Actual RCB class    :", sum(!is.na(df$RCB)),    "\n")
cat("  Hatzis signature    :", sum(!is.na(df$Hatzis)), "\n")
cat("  DLDA30 prediction   :", sum(!is.na(df$DLDA30)), "\n")
cat("  RCB-0/I prediction  :", sum(!is.na(df$RCB0I)),  "\n")

## ---- 1. Distribution tests ---------------------------------------------
do_wilcox <- function(g) {
  d <- df[!is.na(df[[g]]), c("MAP3K8", g)]
  if (length(unique(d[[g]])) != 2) return(NULL)
  w  <- wilcox.test(d$MAP3K8 ~ d[[g]])
  lv <- levels(droplevels(d[[g]]))
  data.frame(
    endpoint = g,
    group1 = lv[1], n1 = sum(d[[g]] == lv[1]),
    median1 = signif(median(d$MAP3K8[d[[g]] == lv[1]]), 4),
    group2 = lv[2], n2 = sum(d[[g]] == lv[2]),
    median2 = signif(median(d$MAP3K8[d[[g]] == lv[2]]), 4),
    Wilcoxon_p = signif(w$p.value, 3))
}
stats_2cls <- bind_rows(lapply(c("pCR","Hatzis","DLDA30","RCB0I"), do_wilcox))
print(stats_2cls)

##  RCB class is 3-level - use Kruskal-Wallis
rcb_d <- df[!is.na(df$RCB), c("MAP3K8","RCB")]
k_rcb <- kruskal.test(MAP3K8 ~ RCB, data = rcb_d)
stats_3cls <- data.frame(
  endpoint = "RCB_class",
  n_RCB_0I  = sum(rcb_d$RCB == "RCB-0/I"),
  median_RCB_0I  = signif(median(rcb_d$MAP3K8[rcb_d$RCB == "RCB-0/I"]),  4),
  n_RCB_II  = sum(rcb_d$RCB == "RCB-II"),
  median_RCB_II  = signif(median(rcb_d$MAP3K8[rcb_d$RCB == "RCB-II"]),  4),
  n_RCB_III = sum(rcb_d$RCB == "RCB-III"),
  median_RCB_III = signif(median(rcb_d$MAP3K8[rcb_d$RCB == "RCB-III"]), 4),
  KruskalWallis_p = signif(k_rcb$p.value, 3))
print(stats_3cls)
write.csv(stats_2cls, "MAP3K8_chemoresponse_stats_2class.csv", row.names = FALSE)
write.csv(stats_3cls, "MAP3K8_chemoresponse_stats_RCB.csv",    row.names = FALSE)

## ---- 2. Boxplots --------------------------------------------------------
make_box <- function(g, title) {
  d <- df[!is.na(df[[g]]), c("MAP3K8", g)]
  ggplot(d, aes(x = .data[[g]], y = MAP3K8, fill = .data[[g]])) +
    geom_boxplot(width = 0.55, alpha = 0.85, outlier.size = 0.7) +
    geom_jitter(width = 0.15, size = 0.8, alpha = 0.55) +
    stat_compare_means(method = "wilcox.test", label.y = max(d$MAP3K8) + 0.3) +
    labs(x = NULL, y = "MAP3K8 log2 expression", title = title) +
    theme_classic(base_size = 11) +
    theme(legend.position = "none",
          plot.title = element_text(size = 11, face = "bold"))
}
make_box_kw <- function(g, title) {
  d <- df[!is.na(df[[g]]), c("MAP3K8", g)]
  ggplot(d, aes(x = .data[[g]], y = MAP3K8, fill = .data[[g]])) +
    geom_boxplot(width = 0.55, alpha = 0.85, outlier.size = 0.7) +
    geom_jitter(width = 0.15, size = 0.8, alpha = 0.55) +
    stat_compare_means(method = "kruskal.test", label.y = max(d$MAP3K8) + 0.3) +
    scale_fill_brewer(palette = "RdYlBu", direction = -1) +
    labs(x = NULL, y = "MAP3K8 log2 expression", title = title) +
    theme_classic(base_size = 11) +
    theme(legend.position = "none",
          plot.title = element_text(size = 11, face = "bold"))
}
pal_pair <- function() scale_fill_manual(values = c("#D7263D", "#1F77B4"))
p1 <- make_box("pCR",    "MAP3K8 by actual pCR vs RD")    + pal_pair()
p2 <- make_box("Hatzis", "MAP3K8 by Hatzis chemosensitivity prediction") + pal_pair()
p3 <- make_box("DLDA30", "MAP3K8 by DLDA30 prediction")   + pal_pair()
p4 <- make_box("RCB0I",  "MAP3K8 by RCB-0/I prediction")  + pal_pair()
p5 <- make_box_kw("RCB", "MAP3K8 by actual RCB class")

panel <- ggarrange(p1, p2, p3, p4, p5,
                   ncol = 2, nrow = 3,
                   labels = c("A","B","C","D","E"))
ggsave("boxplots_MAP3K8_by_class.pdf", panel, width = 11, height = 12)
ggsave("boxplots_MAP3K8_by_class.png", panel, width = 11, height = 12, dpi = 300)

## ---- 3. Logistic regression: pCR ~ MAP3K8 (continuous) -----------------
d_pcr <- df[!is.na(df$pCR), ]
d_pcr$y <- as.integer(d_pcr$pCR == "PCR")
lr <- glm(y ~ MAP3K8, data = d_pcr, family = binomial)
s  <- summary(lr)
OR <- exp(coef(lr))
ci <- suppressMessages(exp(confint(lr)))
log_tab <- data.frame(
  term = c("Intercept","MAP3K8 (per unit log2)"),
  OR   = signif(OR, 3),
  LCL95 = signif(ci[,1], 3),
  UCL95 = signif(ci[,2], 3),
  p    = signif(s$coefficients[,"Pr(>|z|)"], 3))
print(log_tab)
write.csv(log_tab, "MAP3K8_logistic_pCR.csv", row.names = FALSE)

## ---- 4. ROC curves ------------------------------------------------------
##  (a) MAP3K8 alone (continuous)
roc_map  <- suppressMessages(pROC::roc(d_pcr$y, d_pcr$MAP3K8))
auc_map  <- as.numeric(pROC::auc(roc_map))
ci_map   <- as.numeric(pROC::ci.auc(roc_map))

##  (b) Hatzis signature alone (binary -> use as probability of pCR)
d_hat <- df[!is.na(df$pCR) & !is.na(df$Hatzis), ]
d_hat$y    <- as.integer(d_hat$pCR == "PCR")
d_hat$hatz <- as.integer(d_hat$Hatzis == "RX SENSITIVE")
roc_hat <- suppressMessages(pROC::roc(d_hat$y, d_hat$hatz))
auc_hat <- as.numeric(pROC::auc(roc_hat))
ci_hat  <- as.numeric(pROC::ci.auc(roc_hat))

##  (c) Combined (Hatzis + MAP3K8 continuous)
lr2 <- glm(y ~ hatz + MAP3K8, data = d_hat, family = binomial)
d_hat$pred_combined <- predict(lr2, type = "response")
roc_combo <- suppressMessages(pROC::roc(d_hat$y, d_hat$pred_combined))
auc_combo <- as.numeric(pROC::auc(roc_combo))
ci_combo  <- as.numeric(pROC::ci.auc(roc_combo))

##  DeLong test: combined vs Hatzis alone
dl <- suppressMessages(pROC::roc.test(roc_hat, roc_combo, method = "delong"))

cat(sprintf("AUC MAP3K8 alone     : %.3f (95%% CI %.3f-%.3f)\n",
            auc_map, ci_map[1], ci_map[3]))
cat(sprintf("AUC Hatzis alone     : %.3f (95%% CI %.3f-%.3f)\n",
            auc_hat, ci_hat[1], ci_hat[3]))
cat(sprintf("AUC Combined         : %.3f (95%% CI %.3f-%.3f)\n",
            auc_combo, ci_combo[1], ci_combo[3]))
cat(sprintf("DeLong test (combined vs Hatzis alone): P = %.3g\n",
            dl$p.value))

##  ROC plot
roc_df <- rbind(
  data.frame(spec = 1 - roc_map$specificities, sens = roc_map$sensitivities,
             model = sprintf("MAP3K8 alone (AUC %.2f)", auc_map)),
  data.frame(spec = 1 - roc_hat$specificities, sens = roc_hat$sensitivities,
             model = sprintf("Hatzis signature (AUC %.2f)", auc_hat)),
  data.frame(spec = 1 - roc_combo$specificities, sens = roc_combo$sensitivities,
             model = sprintf("Hatzis + MAP3K8 (AUC %.2f)", auc_combo)))
roc_df$model <- factor(roc_df$model, levels = unique(roc_df$model))
roc_plot <- ggplot(roc_df, aes(spec, sens, color = model)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey60") +
  geom_step(linewidth = 0.9) +
  scale_color_manual(values = c("#D7263D","#1F77B4","#2CA02C")) +
  coord_equal() +
  labs(x = "1 - Specificity", y = "Sensitivity",
       title = "ROC: MAP3K8 vs Hatzis chemosensitivity signature for pCR",
       color = NULL) +
  theme_classic(base_size = 11) +
  theme(legend.position = c(0.62, 0.18),
        legend.background = element_rect(fill = alpha("white", 0.7), color = NA))
ggsave("ROC_MAP3K8_vs_signature.pdf", roc_plot, width = 6.5, height = 6.5)
ggsave("ROC_MAP3K8_vs_signature.png", roc_plot, width = 6.5, height = 6.5, dpi = 300)

## ---- 5. Cut-point performance metrics ----------------------------------
classify <- function(thr) {
  pred <- ifelse(d_pcr$MAP3K8 >= thr, "HIGH", "LOW")
  ##  By the DRFS direction, MAP3K8-high = worse outcome => predicts RD
  ##  Sensitivity for RD = TP / (TP + FN), positive = RD, negative = pCR
  TP <- sum(pred == "HIGH" & d_pcr$pCR == "RD")
  TN <- sum(pred == "LOW"  & d_pcr$pCR == "PCR")
  FP <- sum(pred == "HIGH" & d_pcr$pCR == "PCR")
  FN <- sum(pred == "LOW"  & d_pcr$pCR == "RD")
  sens <- TP / (TP + FN)
  spec <- TN / (TN + FP)
  ppv  <- TP / (TP + FP)
  npv  <- TN / (TN + FN)
  data.frame(cutoff = signif(thr, 4),
             TP = TP, TN = TN, FP = FP, FN = FN,
             sensitivity = signif(sens, 3),
             specificity = signif(spec, 3),
             PPV = signif(ppv, 3),
             NPV = signif(npv, 3))
}
##  (i) Median cutpoint
med_cut <- median(d_pcr$MAP3K8)
m_med   <- cbind(label = "median", classify(med_cut))
##  (ii) Optimal cutpoint by Youden (use cutpointr if available, else
##       compute from pROC coordinates)
if (ok_cutpointr) {
  opt <- cutpointr(d_pcr, MAP3K8, pCR,
                   pos_class = "RD", neg_class = "PCR",
                   method = maximize_metric, metric = youden,
                   silent = TRUE)
  opt_cut <- as.numeric(opt$optimal_cutpoint)
} else {
  coords_df <- coords(roc_map, x = "best", best.method = "youden",
                      ret = c("threshold","sensitivity","specificity"),
                      transpose = FALSE)
  opt_cut <- coords_df$threshold[1]
}
m_opt <- cbind(label = "optimal_Youden", classify(opt_cut))

metrics <- rbind(m_med, m_opt)
print(metrics)
write.csv(metrics, "MAP3K8_classifier_metrics.csv", row.names = FALSE)

## ---- 6. Concordance with Hatzis signature ------------------------------
d_conc <- df[!is.na(df$Hatzis), ]
d_conc$MAP3K8_grp <- ifelse(d_conc$MAP3K8 >= median(d_conc$MAP3K8),
                            "MAP3K8-high", "MAP3K8-low")
tbl <- table(d_conc$MAP3K8_grp, d_conc$Hatzis)
print(tbl)
##  Cohen's kappa (manual, since psych::cohen.kappa requires extra dep)
agree <- function(t) {
  n   <- sum(t)
  po  <- sum(diag(t)) / n
  pe  <- sum(rowSums(t) * colSums(t)) / n^2
  (po - pe) / (1 - pe)
}
##  Conceptually MAP3K8-high should align with Rx Insensitive (poor responders)
##  Re-arrange so that MAP3K8-high vs Rx Insensitive is the diagonal
tbl_aligned <- tbl[c("MAP3K8-high","MAP3K8-low"),
                   c("RX INSENSITIVE","RX SENSITIVE")]
kappa <- agree(tbl_aligned)
cat(sprintf("Cohen's kappa (MAP3K8 dichot vs Hatzis): %.3f\n", kappa))

conc_out <- as.data.frame.matrix(tbl_aligned)
conc_out$kappa <- ""
conc_out[1, "kappa"] <- signif(kappa, 3)
write.csv(conc_out, "MAP3K8_vs_Hatzis_concordance.csv")

## ---- Summary --------------------------------------------------------
summary_df <- data.frame(
  metric = c("AUC MAP3K8 alone (95% CI)",
             "AUC Hatzis signature (95% CI)",
             "AUC Combined (95% CI)",
             "DeLong P (combined vs Hatzis alone)",
             "Median cutoff (log2)",
             "Optimal Youden cutoff (log2)",
             "Cohen's kappa MAP3K8 dichot vs Hatzis",
             "Logistic OR pCR per unit MAP3K8",
             "Logistic P pCR ~ MAP3K8"),
  value = c(sprintf("%.3f (%.3f-%.3f)", auc_map,   ci_map[1],   ci_map[3]),
            sprintf("%.3f (%.3f-%.3f)", auc_hat,   ci_hat[1],   ci_hat[3]),
            sprintf("%.3f (%.3f-%.3f)", auc_combo, ci_combo[1], ci_combo[3]),
            signif(dl$p.value, 3),
            signif(med_cut, 4),
            signif(opt_cut, 4),
            signif(kappa, 3),
            signif(OR["MAP3K8"], 3),
            signif(s$coefficients["MAP3K8","Pr(>|z|)"], 3)))
print(summary_df)
write.csv(summary_df, "MAP3K8_chemoresponse_SUMMARY.csv", row.names = FALSE)

cat("\nAnalysis 07 complete. Files in:\n  ", this_dir, "\n", sep="")
