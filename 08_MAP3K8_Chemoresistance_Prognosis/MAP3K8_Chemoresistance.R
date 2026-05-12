## ============================================================
##  Analysis 08 - MAP3K8 marks chemoresistance in chemo-treated TNBC
##
##  Hypothesis tested:
##    Among triple-negative breast cancer patients treated with
##    neoadjuvant anthracycline-taxane chemotherapy, MAP3K8-high
##    tumours carry a chemoresistance phenotype — they relapse
##    distantly after chemotherapy at a higher rate than MAP3K8-low
##    tumours, regardless of their initial pathological response.
##
##  Strategy (three converging lines of evidence):
##    (A) DRFS in the full cohort: MAP3K8-high vs MAP3K8-low
##    (B) Multivariable Cox: DRFS ~ MAP3K8 + age + AJCC stage + RCB,
##        showing the MAP3K8 effect is independent of initial
##        chemoresponse (RCB).
##    (C) Pre-specified subgroup KMs:
##        - RCB-0/I  (chemo-responders): does MAP3K8-high still
##          mark relapse even in patients whose tumour visibly
##          responded?  (STRONGEST chemoresistance test)
##        - RCB-II/III (residual disease): does MAP3K8 stratify
##          further within already-resistant tumours?
##    (D) Forest plot of HRs across the overall + subgroup analyses.
##
##  Cohort : GSE25066 TNBC (n = 178), DRFS endpoint (months)
##
##  Outputs (this folder):
##    KM_MAP3K8_overall.{pdf,png,tiff}
##    KM_MAP3K8_RCB0I.{pdf,png}
##    KM_MAP3K8_RCBII_III.{pdf,png}
##    forest_subgroup_HRs.{pdf,png}
##    cumulative_relapse_rates.csv
##    multivariable_Cox.csv
##    subgroup_HR_table.csv
##    SUMMARY.csv
## ============================================================

need <- c("survival","survminer","ggplot2","dplyr")
for (pkg in need) if (!requireNamespace(pkg, quietly=TRUE))
  install.packages(pkg, repos="https://cloud.r-project.org")
suppressPackageStartupMessages({
  library(survival); library(survminer); library(ggplot2); library(dplyr)
})

## ---- Portable path setup -----------------------------------------------
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
t_months <- as.numeric(clinical_tnbc[["drfs_even_time_years:ch1"]]) * 12
event    <- as.numeric(clinical_tnbc[["drfs_1_event_0_censored:ch1"]])
age      <- suppressWarnings(as.numeric(clinical_tnbc[["age_years:ch1"]]))
stage    <- factor(clinical_tnbc[["clinical_ajcc_stage:ch1"]])
rcb_raw  <- clinical_tnbc[["pathologic_response_rcb_class:ch1"]]

df <- data.frame(
  geo_accession = clinical_tnbc$geo_accession,
  MAP3K8 = as.numeric(map3k8),
  time_months = t_months, event = event,
  age = age, stage = stage,
  RCB_full = factor(rcb_raw),
  stringsAsFactors = FALSE)
df <- df[is.finite(df$MAP3K8) & is.finite(df$time_months) &
         is.finite(df$event), ]

##  Collapse RCB into chemo-responders (RCB-0/I) vs residual (II/III).
df$RCB_cat <- factor(
  ifelse(df$RCB_full == "RCB-0/I",  "Chemo-responders (RCB-0/I)",
  ifelse(df$RCB_full %in% c("RCB-II","RCB-III"),
                                    "Residual disease (RCB-II/III)",
                                    NA)),
  levels = c("Chemo-responders (RCB-0/I)",
             "Residual disease (RCB-II/III)"))

##  MAP3K8 median split, MAP3K8-low as reference
cut <- median(df$MAP3K8)
df$MAP3K8_grp <- factor(ifelse(df$MAP3K8 >= cut, "MAP3K8-high", "MAP3K8-low"),
                        levels = c("MAP3K8-low","MAP3K8-high"))

cat("Cohort n =", nrow(df), "with", sum(df$event), "DRFS events\n")
cat("MAP3K8 median log2 cutoff =", round(cut, 3), "\n")
cat("MAP3K8 group sizes:\n"); print(table(df$MAP3K8_grp))
cat("\nRCB subgroup sizes:\n"); print(table(df$RCB_cat, df$MAP3K8_grp, useNA = "ifany"))

## ============================================================
##  (A) DRFS in the full cohort
## ============================================================
so <- Surv(df$time_months, df$event)
fit_all   <- survfit(so ~ MAP3K8_grp, data = df)
sdiff_all <- survdiff(so ~ MAP3K8_grp, data = df)
p_all_lr  <- 1 - pchisq(sdiff_all$chisq, length(sdiff_all$n) - 1)
cox_all   <- coxph(so ~ MAP3K8_grp, data = df)
s_all     <- summary(cox_all)
HR_all    <- s_all$conf.int[1, "exp(coef)"]
LCL_all   <- s_all$conf.int[1, "lower .95"]
UCL_all   <- s_all$conf.int[1, "upper .95"]
p_all_cox <- s_all$coefficients[1, "Pr(>|z|)"]

cat(sprintf("\n[A] Full cohort: MAP3K8-high vs -low\n"))
cat(sprintf("    HR = %.2f (%.2f-%.2f), log-rank P = %.3g, Cox P = %.3g\n",
            HR_all, LCL_all, UCL_all, p_all_lr, p_all_cox))

##  Cumulative incidence of relapse at 2, 5, 7 years
sm <- summary(fit_all, times = c(24, 60, 84), extend = TRUE)
cum_df <- data.frame(
  group = as.character(sm$strata),
  time_months = sm$time,
  S_t = round(sm$surv, 3),
  CI_lower = round(sm$lower, 3),
  CI_upper = round(sm$upper, 3),
  relapse_pct = round((1 - sm$surv) * 100, 1))
write.csv(cum_df, "cumulative_relapse_rates.csv", row.names = FALSE)
print(cum_df)

## ============================================================
##  (B) Multivariable Cox: MAP3K8 effect independent of RCB
## ============================================================
df_mv <- df[!is.na(df$age) & !is.na(df$stage) & !is.na(df$RCB_full), ]
df_mv$stage <- droplevels(df_mv$stage)
df_mv$RCB_full <- droplevels(df_mv$RCB_full)

cox_mv <- coxph(Surv(time_months, event) ~ MAP3K8_grp + age + stage + RCB_full,
                data = df_mv)
s_mv <- summary(cox_mv)
mv_tab <- data.frame(
  term  = rownames(s_mv$conf.int),
  HR    = signif(s_mv$conf.int[, "exp(coef)"], 3),
  LCL95 = signif(s_mv$conf.int[, "lower .95"], 3),
  UCL95 = signif(s_mv$conf.int[, "upper .95"], 3),
  p     = signif(s_mv$coefficients[, "Pr(>|z|)"], 3))
print(mv_tab)
write.csv(mv_tab, "multivariable_Cox.csv", row.names = FALSE)

HR_adj  <- mv_tab$HR[mv_tab$term == "MAP3K8_grpMAP3K8-high"]
LCL_adj <- mv_tab$LCL95[mv_tab$term == "MAP3K8_grpMAP3K8-high"]
UCL_adj <- mv_tab$UCL95[mv_tab$term == "MAP3K8_grpMAP3K8-high"]
P_adj   <- mv_tab$p[mv_tab$term == "MAP3K8_grpMAP3K8-high"]

cat(sprintf("\n[B] Adjusted for age + stage + RCB: HR = %.2f (%.2f-%.2f), P = %.3g\n",
            HR_adj, LCL_adj, UCL_adj, P_adj))
cat(sprintf("    RCB-III vs RCB-0/I (sanity check): HR = %s, P = %s\n",
            mv_tab$HR[mv_tab$term == "RCB_fullRCB-III"],
            mv_tab$p[mv_tab$term == "RCB_fullRCB-III"]))

## ============================================================
##  (C) Pre-specified subgroup analyses
## ============================================================
run_subgroup <- function(label, dsub) {
  fit <- survfit(Surv(time_months, event) ~ MAP3K8_grp, data = dsub)
  sd  <- survdiff(Surv(time_months, event) ~ MAP3K8_grp, data = dsub)
  p_lr <- 1 - pchisq(sd$chisq, length(sd$n) - 1)
  cox <- coxph(Surv(time_months, event) ~ MAP3K8_grp, data = dsub)
  s   <- summary(cox)
  HR  <- s$conf.int[1, "exp(coef)"]
  LCL <- s$conf.int[1, "lower .95"]
  UCL <- s$conf.int[1, "upper .95"]
  cat(sprintf("\n[C] %s (n=%d, events=%d): HR = %.2f (%.2f-%.2f), log-rank P = %.3g\n",
              label, nrow(dsub), sum(dsub$event), HR, LCL, UCL, p_lr))

  list(fit = fit, dsub = dsub, HR = HR, LCL = LCL, UCL = UCL, p = p_lr,
       n = nrow(dsub), events = sum(dsub$event))
}
res_rcb0I  <- run_subgroup("RCB-0/I  (chemo-responders)",
                           df[df$RCB_cat == "Chemo-responders (RCB-0/I)" & !is.na(df$RCB_cat), ])
res_rcbHi  <- run_subgroup("RCB-II/III (residual disease)",
                           df[df$RCB_cat == "Residual disease (RCB-II/III)" & !is.na(df$RCB_cat), ])

## ============================================================
##  Plots
## ============================================================
make_km <- function(label, fit, dsub, HR, LCL, UCL, p_lr) {
  n_high <- sum(dsub$MAP3K8_grp == "MAP3K8-high")
  n_low  <- sum(dsub$MAP3K8_grp == "MAP3K8-low")
  events <- sum(dsub$event)
  annot  <- sprintf("HR = %.2f (%.2f-%.2f)\nLog-rank P = %.3g",
                    HR, LCL, UCL, p_lr)
  title  <- sprintf("%s\n(n = %d, %d DRFS events)", label, nrow(dsub), events)
  p <- ggsurvplot(
    fit, data = dsub,
    risk.table = TRUE, pval = FALSE, conf.int = FALSE,
    palette = c("#1F77B4", "#D7263D"),
    legend = "top",
    legend.title = "MAP3K8",
    legend.labs = c(paste0("low  (n=", n_low,  ")"),
                    paste0("high (n=", n_high, ")")),
    xlab = "Time (months)", ylab = "DRFS probability",
    title = title,
    risk.table.title = "Number at risk",
    risk.table.height = 0.27,
    ggtheme = theme_classic(base_size = 12),
    font.title = c(11, "bold"))
  p$plot <- p$plot +
    annotate("text",
             x = max(dsub$time_months, na.rm = TRUE) * 0.55, y = 0.92,
             label = annot, hjust = 0, size = 4)
  p
}

##  Full cohort plot
p_all <- make_km("GSE25066 TNBC — full chemo-treated cohort",
                 fit_all, df, HR_all, LCL_all, UCL_all, p_all_lr)
pdf("KM_MAP3K8_overall.pdf", width = 7.8, height = 7.5, useDingbats = FALSE)
print(p_all); dev.off()
png("KM_MAP3K8_overall.png", width = 7.8, height = 7.5, units = "in", res = 300)
print(p_all); dev.off()
tiff("KM_MAP3K8_overall.tiff", width = 7.8, height = 7.5, units = "in", res = 600)
print(p_all); dev.off()

##  RCB-0/I subgroup (chemo responders)
p_r0 <- make_km("RCB-0/I subgroup — chemo responders",
                res_rcb0I$fit, res_rcb0I$dsub,
                res_rcb0I$HR, res_rcb0I$LCL, res_rcb0I$UCL, res_rcb0I$p)
pdf("KM_MAP3K8_RCB0I.pdf", width = 7.5, height = 7.0, useDingbats = FALSE)
print(p_r0); dev.off()
png("KM_MAP3K8_RCB0I.png", width = 7.5, height = 7.0, units = "in", res = 300)
print(p_r0); dev.off()

##  RCB-II/III subgroup (residual disease)
p_rh <- make_km("RCB-II/III subgroup — residual disease",
                res_rcbHi$fit, res_rcbHi$dsub,
                res_rcbHi$HR, res_rcbHi$LCL, res_rcbHi$UCL, res_rcbHi$p)
pdf("KM_MAP3K8_RCBII_III.pdf", width = 7.5, height = 7.0, useDingbats = FALSE)
print(p_rh); dev.off()
png("KM_MAP3K8_RCBII_III.png", width = 7.5, height = 7.0, units = "in", res = 300)
print(p_rh); dev.off()

## ============================================================
##  (D) Forest plot of HRs across overall + subgroups
## ============================================================
sub_df <- data.frame(
  subgroup = c("All chemo-treated TNBC (univariable)",
               "All chemo-treated TNBC (adj. age + stage + RCB)",
               "RCB-0/I (chemo responders)",
               "RCB-II/III (residual disease)"),
  n = c(nrow(df), nrow(df_mv), res_rcb0I$n, res_rcbHi$n),
  events = c(sum(df$event), sum(df_mv$event),
             res_rcb0I$events, res_rcbHi$events),
  HR  = c(HR_all, HR_adj, res_rcb0I$HR, res_rcbHi$HR),
  LCL = c(LCL_all, LCL_adj, res_rcb0I$LCL, res_rcbHi$LCL),
  UCL = c(UCL_all, UCL_adj, res_rcb0I$UCL, res_rcbHi$UCL),
  p   = c(p_all_lr, P_adj, res_rcb0I$p, res_rcbHi$p))
sub_df$label <- sprintf("%s\n(n=%d, %d events)",
                        sub_df$subgroup, sub_df$n, sub_df$events)
sub_df$label <- factor(sub_df$label, levels = rev(sub_df$label))
print(sub_df)
write.csv(sub_df[, c("subgroup","n","events","HR","LCL","UCL","p")],
          "subgroup_HR_table.csv", row.names = FALSE)

fp <- ggplot(sub_df, aes(x = HR, y = label)) +
  geom_vline(xintercept = 1, linetype = 2, color = "grey50") +
  geom_errorbarh(aes(xmin = LCL, xmax = UCL), height = 0.18,
                 color = "#1F2A44", linewidth = 0.7) +
  geom_point(size = 3.6, color = "#D7263D") +
  geom_text(aes(label = sprintf("HR %.2f (%.2f-%.2f)  P=%.3g",
                                HR, LCL, UCL, p)),
            x = max(c(sub_df$UCL, 5)) * 1.1, hjust = 0, size = 3.5,
            color = "#1F2A44") +
  scale_x_log10(limits = c(min(sub_df$LCL) * 0.6,
                           max(c(sub_df$UCL, 5)) * 3.5)) +
  labs(x = "HR for MAP3K8-high vs MAP3K8-low (log scale)", y = NULL,
       title = "MAP3K8-high marks chemoresistance — DRFS HRs across subgroups",
       subtitle = "GSE25066 TNBC, neoadjuvant anthracycline-taxane (n=178)") +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(size = 12, face = "bold"),
        plot.subtitle = element_text(size = 10, color = "grey40"))
ggsave("forest_subgroup_HRs.pdf", fp, width = 9, height = 4.0)
ggsave("forest_subgroup_HRs.png", fp, width = 9, height = 4.0, dpi = 300)

## ============================================================
##  Summary
## ============================================================
summary_df <- data.frame(
  metric = c(
    "Cohort n",
    "Events (distant relapse)",
    "MAP3K8 cutoff (median log2)",
    "Full cohort: HR high vs low",
    "Full cohort: log-rank P",
    "Adjusted (age + stage + RCB): HR high vs low",
    "Adjusted: P",
    "RCB-0/I subgroup: HR high vs low",
    "RCB-0/I subgroup: log-rank P",
    "RCB-II/III subgroup: HR high vs low",
    "RCB-II/III subgroup: log-rank P",
    "Cumulative relapse @ 60mo (MAP3K8-high)",
    "Cumulative relapse @ 60mo (MAP3K8-low)"),
  value = c(
    nrow(df),
    sum(df$event),
    signif(cut, 4),
    sprintf("%.2f (%.2f-%.2f)", HR_all, LCL_all, UCL_all),
    signif(p_all_lr, 3),
    sprintf("%.2f (%.2f-%.2f)", HR_adj, LCL_adj, UCL_adj),
    signif(P_adj, 3),
    sprintf("%.2f (%.2f-%.2f)", res_rcb0I$HR, res_rcb0I$LCL, res_rcb0I$UCL),
    signif(res_rcb0I$p, 3),
    sprintf("%.2f (%.2f-%.2f)", res_rcbHi$HR, res_rcbHi$LCL, res_rcbHi$UCL),
    signif(res_rcbHi$p, 3),
    sprintf("%.1f%%", subset(cum_df, time_months == 60 &
                              grepl("high", group))$relapse_pct[1]),
    sprintf("%.1f%%", subset(cum_df, time_months == 60 &
                              grepl("low",  group))$relapse_pct[1])))
print(summary_df)
write.csv(summary_df, "SUMMARY.csv", row.names = FALSE)

cat("\nAnalysis 08 complete. Files in:\n  ", this_dir, "\n", sep = "")
