## ============================================================
##  Analysis 06a - TROP2 (TACSTD2) single-gene DRFS Kaplan-Meier
##  Cohort : GSE25066 TNBC (n = 178), DRFS endpoint
##  Design : median split of TROP2 expression (probe selected as
##           the higher-mean of the three GPL96 TACSTD2 probes
##           202285_s_at / 202286_s_at / 202287_s_at within the
##           TNBC subset). TROP2-low is the reference, so the
##           HR is reported for high vs low.
##  Outputs (this folder):
##    KMplot_TROP2.{pdf,png,tiff}
##    TROP2_KM_stats.csv
## ============================================================

need <- c("survival","survminer","ggplot2")
for (pkg in need) if (!requireNamespace(pkg, quietly=TRUE))
  install.packages(pkg, repos="https://cloud.r-project.org")
suppressPackageStartupMessages({
  library(survival); library(survminer); library(ggplot2)
})

## ---- Portable path setup ------------------------------------------------
this_dir <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", a, value = TRUE)
  if (length(f)) {
    p <- sub("^--file=", "", f[1])
    p <- gsub("~\\+~", " ", p)   # Rscript URL-encodes spaces as ~+~
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
  stop("Workspace not found at: ", WS,
       "\nSet TNBC_WORKSPACE env var or place the file at the default path.")
load(WS)

## ---- Pick higher-mean TROP2 probe --------------------------------------
sym_col      <- "Gene Symbol"
trop2_rows   <- which(gpl_table[[sym_col]] == "TACSTD2")
trop2_probes <- as.character(gpl_table$ID[trop2_rows])
stopifnot(length(trop2_probes) >= 1)
trop2_probe  <- names(which.max(rowMeans(expr_tnbc[trop2_probes, , drop = FALSE])))
cat("TROP2 probes available:", paste(trop2_probes, collapse = ", "), "\n")
cat("TROP2 probe used:", trop2_probe,
    "(highest mean log2 in TNBC cohort)\n")
trop2 <- as.numeric(expr_tnbc[trop2_probe, ])

## ---- Build vectors ------------------------------------------------------
stopifnot(identical(rownames(clinical_tnbc), colnames(expr_tnbc)))
t_months <- as.numeric(clinical_tnbc[["drfs_even_time_years:ch1"]]) * 12
event    <- as.numeric(clinical_tnbc[["drfs_1_event_0_censored:ch1"]])

df <- data.frame(geo_accession = clinical_tnbc$geo_accession,
                 TROP2 = trop2, time_months = t_months, event = event)
df <- df[is.finite(df$TROP2) & is.finite(df$time_months) &
         is.finite(df$event), ]
cut <- median(df$TROP2)
df$TROP2_grp <- factor(ifelse(df$TROP2 >= cut, "high", "low"),
                       levels = c("low","high"))   # low = reference

## ---- Models -------------------------------------------------------------
so   <- Surv(df$time_months, df$event)
fit  <- survfit(so ~ TROP2_grp, data = df)
sd   <- survdiff(so ~ TROP2_grp, data = df)
p_lr <- 1 - pchisq(sd$chisq, length(sd$n)-1)
cox  <- coxph(so ~ TROP2_grp, data = df)
s    <- summary(cox)
HR   <- s$conf.int[1, "exp(coef)"]
LCL  <- s$conf.int[1, "lower .95"]
UCL  <- s$conf.int[1, "upper .95"]
cat(sprintf("TROP2 high vs low — HR = %.2f (%.2f-%.2f), log-rank P = %.3g\n",
            HR, LCL, UCL, p_lr))

stats <- data.frame(
  metric = c("n_total","events","n_high","n_low","cutoff_log2",
             "HR_high_vs_low","LCL95","UCL95","logrank_P","probe_used"),
  value  = c(nrow(df), sum(df$event),
             sum(df$TROP2_grp == "high"), sum(df$TROP2_grp == "low"),
             signif(cut,4), signif(HR,3), signif(LCL,3),
             signif(UCL,3), signif(p_lr,3), trop2_probe))
write.csv(stats, "TROP2_KM_stats.csv", row.names = FALSE)

## ---- Plot ---------------------------------------------------------------
title_block <- sprintf(
  "GSE25066 - TNBC (n = %d)\nDRFS by TROP2 / TACSTD2 (median split, probe %s)",
  nrow(df), trop2_probe)
annot <- sprintf("HR = %.2f (%.2f-%.2f)\nLog-rank P = %.3g", HR, LCL, UCL, p_lr)

p <- ggsurvplot(
  fit, data = df,
  risk.table = TRUE, pval = FALSE, conf.int = FALSE,
  palette = c("#1F77B4", "#D7263D"),
  legend = "top",
  legend.title = "TROP2",
  legend.labs = c(paste0("low  (n=", sum(df$TROP2_grp == "low"),  ")"),
                  paste0("high (n=", sum(df$TROP2_grp == "high"), ")")),
  xlab = "Time (months)", ylab = "DRFS probability",
  title = title_block,
  risk.table.title = "Number at risk",
  risk.table.height = 0.27,
  ggtheme = theme_classic(base_size = 12),
  font.title = c(11, "bold"))
p$plot <- p$plot +
  annotate("text",
           x = max(df$time_months, na.rm = TRUE) * 0.55, y = 0.92,
           label = annot, hjust = 0, size = 4)

pdf("KMplot_TROP2.pdf", width = 7.5, height = 7.5, useDingbats = FALSE)
print(p); dev.off()
png("KMplot_TROP2.png", width = 7.5, height = 7.5, units = "in", res = 300)
print(p); dev.off()
tiff("KMplot_TROP2.tiff", width = 7.5, height = 7.5, units = "in", res = 600)
print(p); dev.off()
cat("TROP2 KM done. Files in:", this_dir, "\n")
