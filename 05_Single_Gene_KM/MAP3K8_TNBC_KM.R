## ============================================================
##  Analysis 05a - MAP3K8 single-gene DRFS Kaplan-Meier
##  Cohort : GSE25066 TNBC (n = 178), DRFS endpoint
##  Design : median split of MAP3K8 (probe 205027_s_at)
##           with MAP3K8-low as the reference group, so that
##           the reported HR matches the convention "high vs low".
##  Outputs (this folder):
##    KMplot_MAP3K8.{pdf,png,tiff}
##    MAP3K8_KM_stats.csv
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
  if (length(f)) return(dirname(normalizePath(sub("^--file=", "", f[1]))))
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

## ---- Build vectors ------------------------------------------------------
stopifnot(identical(rownames(clinical_tnbc), colnames(expr_tnbc)))
m         <- as.numeric(map3k8)
t_months  <- as.numeric(clinical_tnbc[["drfs_even_time_years:ch1"]]) * 12
event     <- as.numeric(clinical_tnbc[["drfs_1_event_0_censored:ch1"]])

df <- data.frame(geo_accession = clinical_tnbc$geo_accession,
                 MAP3K8 = m, time_months = t_months, event = event)
df <- df[is.finite(df$MAP3K8) & is.finite(df$time_months) &
         is.finite(df$event), ]
cut <- median(df$MAP3K8)
df$MAP3K8_grp <- factor(ifelse(df$MAP3K8 >= cut, "high", "low"),
                        levels = c("low","high"))   # low = reference

## ---- Models -------------------------------------------------------------
so   <- Surv(df$time_months, df$event)
fit  <- survfit(so ~ MAP3K8_grp, data = df)
sd   <- survdiff(so ~ MAP3K8_grp, data = df)
p_lr <- 1 - pchisq(sd$chisq, length(sd$n)-1)
cox  <- coxph(so ~ MAP3K8_grp, data = df)
s    <- summary(cox)
HR   <- s$conf.int[1, "exp(coef)"]
LCL  <- s$conf.int[1, "lower .95"]
UCL  <- s$conf.int[1, "upper .95"]
cat(sprintf("MAP3K8 high vs low — HR = %.2f (%.2f-%.2f), log-rank P = %.3g\n",
            HR, LCL, UCL, p_lr))

stats <- data.frame(
  metric = c("n_total","events","n_high","n_low","cutoff_log2",
             "HR_high_vs_low","LCL95","UCL95","logrank_P"),
  value  = c(nrow(df), sum(df$event),
             sum(df$MAP3K8_grp == "high"), sum(df$MAP3K8_grp == "low"),
             signif(cut,4), signif(HR,3), signif(LCL,3),
             signif(UCL,3), signif(p_lr,3)))
write.csv(stats, "MAP3K8_KM_stats.csv", row.names = FALSE)

## ---- Plot ---------------------------------------------------------------
title_block <- sprintf(
  "GSE25066 - TNBC (n = %d)\nDRFS by MAP3K8 (median split, probe 205027_s_at)",
  nrow(df))
annot <- sprintf("HR = %.2f (%.2f-%.2f)\nLog-rank P = %.3g", HR, LCL, UCL, p_lr)

p <- ggsurvplot(
  fit, data = df,
  risk.table = TRUE, pval = FALSE, conf.int = FALSE,
  palette = c("#1F77B4", "#D7263D"),
  legend = "top",
  legend.title = "MAP3K8",
  legend.labs = c(paste0("low  (n=", sum(df$MAP3K8_grp == "low"),  ")"),
                  paste0("high (n=", sum(df$MAP3K8_grp == "high"), ")")),
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

pdf("KMplot_MAP3K8.pdf", width = 7.5, height = 7.5, useDingbats = FALSE)
print(p); dev.off()
png("KMplot_MAP3K8.png", width = 7.5, height = 7.5, units = "in", res = 300)
print(p); dev.off()
tiff("KMplot_MAP3K8.tiff", width = 7.5, height = 7.5, units = "in", res = 600)
print(p); dev.off()
cat("MAP3K8 KM done. Files in:", this_dir, "\n")
