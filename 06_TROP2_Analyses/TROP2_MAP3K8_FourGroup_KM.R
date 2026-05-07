## ============================================================
##  Analysis 06b - TROP2 x MAP3K8 four-group interaction KM
##  Cohort : GSE25066 TNBC (n = 178), DRFS endpoint (months)
##  Design :
##    Median split TROP2 (high/low) x median split MAP3K8 (high/low)
##    -> 4 groups: TROP2_low / MAP3K8_low      (reference)
##                 TROP2_low / MAP3K8_high
##                 TROP2_high / MAP3K8_low
##                 TROP2_high / MAP3K8_high
##  Tests:
##    - Overall log-rank across 4 strata
##    - Pairwise log-rank (Benjamini-Hochberg adjusted)
##    - Cox model with interaction: DRFS ~ TROP2 * MAP3K8
##    - Group4-vs-reference Cox
##  Outputs (this folder):
##    KMplot_TROP2_MAP3K8_4groups.{pdf,png,tiff}
##    forest_TROP2_MAP3K8_4groups.{pdf,png}
##    interaction_TROP2_MAP3K8.csv
##    pairwise_TROP2_MAP3K8.csv
##    summary_TROP2_MAP3K8.csv
##    group4_vs_reference_TROP2_MAP3K8.csv
## ============================================================

need <- c("survival","survminer","ggplot2","dplyr")
for (pkg in need) if (!requireNamespace(pkg, quietly=TRUE))
  install.packages(pkg, repos="https://cloud.r-project.org")
suppressPackageStartupMessages({
  library(survival); library(survminer); library(ggplot2); library(dplyr)
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

## ---- Pick TROP2 probe ---------------------------------------------------
sym_col      <- "Gene Symbol"
trop2_rows   <- which(gpl_table[[sym_col]] == "TACSTD2")
trop2_probes <- as.character(gpl_table$ID[trop2_rows])
trop2_probe  <- names(which.max(rowMeans(expr_tnbc[trop2_probes, , drop = FALSE])))
cat("TROP2 probe used:", trop2_probe, "\n")
trop2 <- as.numeric(expr_tnbc[trop2_probe, ])

stopifnot(identical(rownames(clinical_tnbc), colnames(expr_tnbc)))

t_months <- as.numeric(clinical_tnbc[["drfs_even_time_years:ch1"]]) * 12
event    <- as.numeric(clinical_tnbc[["drfs_1_event_0_censored:ch1"]])

df <- data.frame(geo_accession = clinical_tnbc$geo_accession,
                 TROP2 = trop2, MAP3K8 = as.numeric(map3k8),
                 time_months = t_months, event = event)
df <- df[is.finite(df$TROP2) & is.finite(df$MAP3K8) &
         is.finite(df$time_months) & is.finite(df$event), ]

## ---- Median splits + 4-group factor ------------------------------------
t_cut <- median(df$TROP2)
m_cut <- median(df$MAP3K8)
df$TROP2_grp  <- factor(ifelse(df$TROP2  >= t_cut, "TROP2_high",  "TROP2_low"),
                        levels = c("TROP2_low","TROP2_high"))
df$MAP3K8_grp <- factor(ifelse(df$MAP3K8 >= m_cut, "MAP3K8_high", "MAP3K8_low"),
                        levels = c("MAP3K8_low","MAP3K8_high"))

df$Group4 <- factor(
  paste(df$TROP2_grp, df$MAP3K8_grp, sep = " / "),
  levels = c("TROP2_low / MAP3K8_low",
             "TROP2_low / MAP3K8_high",
             "TROP2_high / MAP3K8_low",
             "TROP2_high / MAP3K8_high"))
cat("4-group sizes:\n"); print(table(df$Group4))

## ---- Survival models ----------------------------------------------------
so <- Surv(df$time_months, df$event)

fit  <- survfit(so ~ Group4, data = df)
ovr  <- survdiff(so ~ Group4, data = df)
p_overall <- 1 - pchisq(ovr$chisq, length(ovr$n) - 1)

pw <- pairwise_survdiff(Surv(time_months, event) ~ Group4, data = df,
                        p.adjust.method = "BH")
cat("\nPairwise log-rank (BH):\n"); print(pw$p.value)

cox_int  <- coxph(so ~ TROP2_grp * MAP3K8_grp, data = df)
s_int    <- summary(cox_int)
cox_main <- coxph(so ~ TROP2_grp + MAP3K8_grp, data = df)
lrt      <- anova(cox_main, cox_int, test = "Chisq")
p_interaction <- lrt[2, "Pr(>|Chi|)"]

hr_tab <- data.frame(
  term  = rownames(s_int$conf.int),
  HR    = signif(s_int$conf.int[, "exp(coef)"], 3),
  LCL95 = signif(s_int$conf.int[, "lower .95"], 3),
  UCL95 = signif(s_int$conf.int[, "upper .95"], 3),
  p     = signif(s_int$coefficients[, "Pr(>|z|)"], 3))
print(hr_tab)
write.csv(hr_tab, "interaction_TROP2_MAP3K8.csv", row.names = FALSE)
write.csv(as.data.frame(pw$p.value), "pairwise_TROP2_MAP3K8.csv")

summary_stats <- data.frame(
  metric = c("n_total",
             "n_TROP2low_MAP3K8low",
             "n_TROP2low_MAP3K8high",
             "n_TROP2high_MAP3K8low",
             "n_TROP2high_MAP3K8high",
             "events_total",
             "TROP2_probe",
             "logrank_overall_p",
             "interaction_LRT_p"),
  value  = c(nrow(df),
             unname(table(df$Group4)["TROP2_low / MAP3K8_low"]),
             unname(table(df$Group4)["TROP2_low / MAP3K8_high"]),
             unname(table(df$Group4)["TROP2_high / MAP3K8_low"]),
             unname(table(df$Group4)["TROP2_high / MAP3K8_high"]),
             sum(df$event),
             trop2_probe,
             signif(p_overall, 3),
             signif(p_interaction, 3)))
write.csv(summary_stats, "summary_TROP2_MAP3K8.csv", row.names = FALSE)

## ---- Group4 vs reference HRs (cleaner readout) -------------------------
cox_g4 <- coxph(so ~ Group4, data = df)
sg4    <- summary(cox_g4)
g4_tab <- data.frame(
  contrast = sub("^Group4", "", rownames(sg4$conf.int)),
  HR    = signif(sg4$conf.int[, "exp(coef)"], 3),
  LCL95 = signif(sg4$conf.int[, "lower .95"], 3),
  UCL95 = signif(sg4$conf.int[, "upper .95"], 3),
  p     = signif(sg4$coefficients[, "Pr(>|z|)"], 3))
print(g4_tab)
write.csv(g4_tab, "group4_vs_reference_TROP2_MAP3K8.csv", row.names = FALSE)

##  Forest plot of group-vs-reference HRs
fp_df <- g4_tab
fp_df$contrast <- factor(fp_df$contrast, levels = rev(fp_df$contrast))
fp <- ggplot(fp_df, aes(x = HR, y = contrast)) +
  geom_vline(xintercept = 1, linetype = 2, color = "grey50") +
  geom_errorbarh(aes(xmin = LCL95, xmax = UCL95), height = 0.2) +
  geom_point(size = 3.2, color = "#D7263D") +
  scale_x_log10() +
  labs(x = "HR vs TROP2_low / MAP3K8_low (log scale)", y = NULL,
       title = "DRFS HR by TROP2 x MAP3K8 group (GSE25066 TNBC, n=178)") +
  theme_classic(base_size = 11)
ggsave("forest_TROP2_MAP3K8_4groups.pdf", fp, width = 7.2, height = 3.6)
ggsave("forest_TROP2_MAP3K8_4groups.png", fp, width = 7.2, height = 3.6, dpi = 300)

## ---- KM plot ------------------------------------------------------------
title_block <- paste0(
  "GSE25066 - TNBC (n = ", nrow(df), ")\n",
  "DRFS by TROP2 x MAP3K8 median-split groups")
annot_txt <- sprintf("Overall log-rank P = %.3g\nInteraction (LRT) P = %.3g",
                     p_overall, p_interaction)

short_labs <- c("TROP2↓  MAP3K8↓",
                "TROP2↓  MAP3K8↑",
                "TROP2↑  MAP3K8↓",
                "TROP2↑  MAP3K8↑")

p <- ggsurvplot(
  fit, data = df,
  risk.table = TRUE, pval = FALSE, conf.int = FALSE,
  palette = c("#1F77B4", "#2CA02C", "#FF7F0E", "#D7263D"),
  legend = "bottom",
  legend.title = "Group",
  legend.labs  = short_labs,
  xlab = "Time (months)",
  ylab = "DRFS probability",
  title = title_block,
  risk.table.title = "Number at risk",
  risk.table.height = 0.30,
  ggtheme = theme_classic(base_size = 11),
  font.title = c(11, "bold"),
  font.legend = c(10, "plain"))
p$plot <- p$plot +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE)) +
  annotate("text",
           x = max(df$time_months, na.rm = TRUE) * 0.45, y = 0.95,
           label = annot_txt, hjust = 0, size = 3.5)

pdf("KMplot_TROP2_MAP3K8_4groups.pdf", width = 8.5, height = 8.0, useDingbats = FALSE)
print(p); dev.off()
png("KMplot_TROP2_MAP3K8_4groups.png", width = 8.5, height = 8.0, units = "in", res = 300)
print(p); dev.off()
tiff("KMplot_TROP2_MAP3K8_4groups.tiff", width = 8.5, height = 8.0, units = "in", res = 600)
print(p); dev.off()
cat("TROP2 x MAP3K8 4-group KM done. Files in:", this_dir, "\n")
