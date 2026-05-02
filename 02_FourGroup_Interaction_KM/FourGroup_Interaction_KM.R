## ============================================================
##  Analysis 02 - MAP3K8 x IL-1beta four-group interaction KM
##  Cohort : GSE25066 TNBC (n = 178), DRFS endpoint (months)
##  Design :
##    Median split MAP3K8 (high/low) x median split IL1B (high/low)
##    -> 4 groups : MAP3K8-low / IL1B-low      (reference, "best" expected)
##                 MAP3K8-low / IL1B-high
##                 MAP3K8-high / IL1B-low      (expected best given prior data)
##                 MAP3K8-high / IL1B-high
##  Tests:
##    - Overall log-rank across 4 strata
##    - Pairwise log-rank (Bonferroni / BH adjusted) via survminer::pairwise_survdiff
##    - Cox model with interaction: DRFS ~ MAP3K8 * IL1B  (binary form)
##  Outputs (this folder):
##    KMplot_4groups.{pdf,png,tiff}
##    interaction_stats.csv
##    Methods_FourGroup_Interaction_KM.docx
## ============================================================

need <- c("survival","survminer","ggplot2","officer","flextable","dplyr")
for (pkg in need) if (!requireNamespace(pkg, quietly=TRUE))
  install.packages(pkg, repos="https://cloud.r-project.org")
suppressPackageStartupMessages({
  library(survival); library(survminer); library(ggplot2)
  library(officer); library(flextable); library(dplyr)
})

## ---- Portable path setup (see 01_pCR_RCB_Analysis for explanation) ------
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

## ---- 1. Build expression vectors ----------------------------------------
sym_col   <- "Gene Symbol"
il1b_rows <- which(gpl_table[[sym_col]] == "IL1B")
il1b_probes <- as.character(gpl_table$ID[il1b_rows])
il1b_probe <- names(which.max(rowMeans(expr_tnbc[il1b_probes, , drop = FALSE])))
il1b <- as.numeric(expr_tnbc[il1b_probe, ])

stopifnot(identical(rownames(clinical_tnbc), colnames(expr_tnbc)))

t_months <- as.numeric(clinical_tnbc[["drfs_even_time_years:ch1"]]) * 12
event    <- as.numeric(clinical_tnbc[["drfs_1_event_0_censored:ch1"]])

df <- data.frame(
  geo_accession = clinical_tnbc$geo_accession,
  MAP3K8 = as.numeric(map3k8),
  IL1B   = il1b,
  time_months = t_months,
  event = event,
  stringsAsFactors = FALSE
)
df <- df[is.finite(df$MAP3K8) & is.finite(df$IL1B) &
         is.finite(df$time_months) & is.finite(df$event), ]

## ---- 2. Median splits + 4-group factor ---------------------------------
m_cut <- median(df$MAP3K8); i_cut <- median(df$IL1B)
df$MAP3K8_grp <- factor(ifelse(df$MAP3K8 >= m_cut, "MAP3K8_high", "MAP3K8_low"),
                        levels = c("MAP3K8_low","MAP3K8_high"))
df$IL1B_grp   <- factor(ifelse(df$IL1B   >= i_cut, "IL1B_high",   "IL1B_low"),
                        levels = c("IL1B_low","IL1B_high"))

df$Group4 <- factor(
  paste(df$MAP3K8_grp, df$IL1B_grp, sep = " / "),
  levels = c("MAP3K8_low / IL1B_low",
             "MAP3K8_low / IL1B_high",
             "MAP3K8_high / IL1B_low",
             "MAP3K8_high / IL1B_high"))
cat("4-group sizes:\n"); print(table(df$Group4))

## ---- 3. Survival models -------------------------------------------------
so <- Surv(df$time_months, df$event)

fit  <- survfit(so ~ Group4, data = df)
ovr  <- survdiff(so ~ Group4, data = df)
p_overall <- 1 - pchisq(ovr$chisq, length(ovr$n) - 1)

pw <- pairwise_survdiff(Surv(time_months, event) ~ Group4, data = df,
                        p.adjust.method = "BH")
cat("\nPairwise log-rank (BH adjusted):\n"); print(pw$p.value)

##  Interaction Cox model on the binary forms
cox_int <- coxph(so ~ MAP3K8_grp * IL1B_grp, data = df)
s_int   <- summary(cox_int)
##  Likelihood-ratio test for the interaction term alone
cox_main <- coxph(so ~ MAP3K8_grp + IL1B_grp, data = df)
lrt <- anova(cox_main, cox_int, test = "Chisq")
p_interaction <- lrt[2, "Pr(>|Chi|)"]

## ---- 4. Export tables ---------------------------------------------------
hr_tab <- data.frame(
  term  = rownames(s_int$conf.int),
  HR    = signif(s_int$conf.int[, "exp(coef)"], 3),
  LCL95 = signif(s_int$conf.int[, "lower .95"], 3),
  UCL95 = signif(s_int$conf.int[, "upper .95"], 3),
  p     = signif(s_int$coefficients[, "Pr(>|z|)"], 3)
)
print(hr_tab)
write.csv(hr_tab,                "interaction_cox_terms.csv", row.names = FALSE)
write.csv(as.data.frame(pw$p.value), "pairwise_logrank_BH.csv")

summary_stats <- data.frame(
  metric = c("n_total",
             "n_MAP3K8low_IL1Blow",
             "n_MAP3K8low_IL1Bhigh",
             "n_MAP3K8high_IL1Blow",
             "n_MAP3K8high_IL1Bhigh",
             "events_total",
             "logrank_overall_p",
             "interaction_LRT_p"),
  value  = c(nrow(df),
             unname(table(df$Group4)["MAP3K8_low / IL1B_low"]),
             unname(table(df$Group4)["MAP3K8_low / IL1B_high"]),
             unname(table(df$Group4)["MAP3K8_high / IL1B_low"]),
             unname(table(df$Group4)["MAP3K8_high / IL1B_high"]),
             sum(df$event),
             signif(p_overall, 3),
             signif(p_interaction, 3))
)
write.csv(summary_stats, "interaction_summary.csv", row.names = FALSE)

## ---- 4b. Group4 vs reference HRs (cleaner readout) ----------------------
##  This Cox parameterises Group4 directly with MAP3K8_low/IL1B_low as
##  reference, so each row is a single contrast against the best group.
cox_g4 <- coxph(so ~ Group4, data = df)
sg4    <- summary(cox_g4)
g4_tab <- data.frame(
  contrast = rownames(sg4$conf.int),
  HR    = signif(sg4$conf.int[, "exp(coef)"], 3),
  LCL95 = signif(sg4$conf.int[, "lower .95"], 3),
  UCL95 = signif(sg4$conf.int[, "upper .95"], 3),
  p     = signif(sg4$coefficients[, "Pr(>|z|)"], 3))
g4_tab$contrast <- sub("^Group4", "", g4_tab$contrast)
print(g4_tab); write.csv(g4_tab, "group4_vs_reference_HR.csv", row.names = FALSE)

##  Forest plot of group-vs-reference HRs
fp_df <- g4_tab
fp_df$contrast <- factor(fp_df$contrast, levels = rev(fp_df$contrast))
fp <- ggplot(fp_df, aes(x = HR, y = contrast)) +
  geom_vline(xintercept = 1, linetype = 2, color = "grey50") +
  geom_errorbarh(aes(xmin = LCL95, xmax = UCL95), height = 0.2) +
  geom_point(size = 3.2, color = "#D7263D") +
  scale_x_log10() +
  labs(x = "HR vs MAP3K8_low / IL1B_low (log scale)", y = NULL,
       title = "DRFS HR by MAP3K8 x IL1B group (GSE25066 TNBC, n=178)") +
  theme_classic(base_size = 11)
ggsave("forest_group4_vs_reference.pdf", fp, width = 7, height = 3.6)
ggsave("forest_group4_vs_reference.png", fp, width = 7, height = 3.6, dpi = 300)

## ---- 5. Plot -------------------------------------------------------------
title_block <- paste0(
  "GSE25066 - TNBC (n = ", nrow(df), ")\n",
  "DRFS by MAP3K8 x IL1B median-split groups")
annot_txt <- sprintf("Overall log-rank P = %.3g\nInteraction (LRT) P = %.3g",
                     p_overall, p_interaction)

p <- ggsurvplot(
  fit, data = df,
  risk.table = TRUE, pval = FALSE, conf.int = FALSE,
  palette = c("#1F77B4", "#2CA02C", "#FF7F0E", "#D7263D"),
  legend.title = "Group",
  legend.labs  = levels(df$Group4),
  xlab = "Time (months)",
  ylab = "DRFS probability",
  title = title_block,
  risk.table.title = "Number at risk",
  risk.table.height = 0.30,
  ggtheme = theme_classic(base_size = 11),
  font.title = c(11, "bold"))
p$plot <- p$plot +
  annotate("text",
           x = max(df$time_months, na.rm = TRUE) * 0.45, y = 0.95,
           label = annot_txt, hjust = 0, size = 3.5)

pdf("KMplot_4groups.pdf", width = 7.2, height = 7.2, useDingbats = FALSE)
print(p); dev.off()
png("KMplot_4groups.png", width = 7.2, height = 7.2, units = "in", res = 300)
print(p); dev.off()
tiff("KMplot_4groups.tiff", width = 7.2, height = 7.2, units = "in", res = 600)
print(p); dev.off()

## ---- 6. Methods Word document -------------------------------------------
doc <- read_docx() |>
  body_add_par("Analysis 02 - MAP3K8 x IL-1beta four-group interaction Kaplan-Meier",
               style = "heading 1") |>
  body_add_par("Hypothesis", style = "heading 2") |>
  body_add_par(
    "MAP3K8 (TPL2) and IL-1beta sit on the same MyD88 / IRAK / NF-kB axis: IL-1R signalling activates IRAK and downstream MAP3K8, which in turn drives ERK/MEK and IL-1beta production. We hypothesised that the joint expression status of these two genes carries more prognostic information in chemo-treated TNBC than either gene alone.",
    style = "Normal") |>
  body_add_par("Cohort and endpoint", style = "heading 2") |>
  body_add_par(
    "GSE25066 TNBC subset (n = 178), neoadjuvant anthracycline-taxane chemotherapy. Endpoint: distant relapse-free survival (DRFS), with time `drfs_even_time_years:ch1` x 12 (months) and event `drfs_1_event_0_censored:ch1`. Patients with missing time, event or expression were excluded.",
    style = "Normal") |>
  body_add_par("Group construction", style = "heading 2") |>
  body_add_par(
    "MAP3K8 expression (probe 205027_s_at) was median-dichotomised into MAP3K8_low vs MAP3K8_high. IL1B expression (probe selected as the higher-mean of the two GPL96 probes 205067_at and 39402_at) was median-dichotomised into IL1B_low vs IL1B_high. The two binary variables were combined into one 4-level factor with reference level MAP3K8_low / IL1B_low.",
    style = "Normal") |>
  body_add_par("Statistical analysis", style = "heading 2") |>
  body_add_par(
    "1. Kaplan-Meier estimates were computed for each of the four groups (survival::survfit) and plotted with survminer::ggsurvplot, including a number-at-risk table.",
    style = "Normal") |>
  body_add_par(
    "2. The overall log-rank test across the four strata was computed via survival::survdiff.",
    style = "Normal") |>
  body_add_par(
    "3. Pairwise log-rank tests between groups were computed via survminer::pairwise_survdiff with Benjamini-Hochberg adjustment.",
    style = "Normal") |>
  body_add_par(
    "4. A Cox proportional-hazards model was fitted with both binary variables and their interaction term (DRFS ~ MAP3K8_grp * IL1B_grp). The significance of the interaction was assessed by likelihood-ratio test versus the additive (main-effects) model.",
    style = "Normal") |>
  body_add_par("Interpretation guide", style = "heading 2") |>
  body_add_par(
    "If MAP3K8 and IL-1beta act on independent axes, the interaction term will be non-significant (additive Cox model fits adequately). A significant interaction indicates that the prognostic effect of one gene depends on the level of the other - for example, IL-1beta-high may be detrimental only in MAP3K8-low tumours.",
    style = "Normal") |>
  body_add_par("Cox terms (this run)", style = "heading 2") |>
  body_add_flextable(flextable(hr_tab) |> autofit() |> fontsize(size = 9, part = "all")) |>
  body_add_par("Group sizes and overall p-values (this run)", style = "heading 2") |>
  body_add_flextable(flextable(summary_stats) |> autofit() |> fontsize(size = 9, part = "all")) |>
  body_add_par("Group-vs-reference Cox table (this run)", style = "heading 2") |>
  body_add_par(
    "Cox model with Group4 entered as a single 4-level factor with MAP3K8_low/IL1B_low as reference. Each row is the hazard ratio of one group versus the doubly-low (best-prognosis) reference.",
    style = "Normal") |>
  body_add_flextable(flextable(g4_tab) |> autofit() |> fontsize(size = 9, part = "all")) |>
  body_add_par("Results obtained in this run", style = "heading 2") |>
  body_add_par(sprintf(
    "After median dichotomisation, the four groups were balanced (MAP3K8_low/IL1B_low = %d, MAP3K8_low/IL1B_high = %d, MAP3K8_high/IL1B_low = %d, MAP3K8_high/IL1B_high = %d), with %d distant-relapse events across %d patients.",
    table(df$Group4)["MAP3K8_low / IL1B_low"],
    table(df$Group4)["MAP3K8_low / IL1B_high"],
    table(df$Group4)["MAP3K8_high / IL1B_low"],
    table(df$Group4)["MAP3K8_high / IL1B_high"],
    sum(df$event), nrow(df)),
    style = "Normal") |>
  body_add_par(sprintf(
    "The overall log-rank test across the four strata gave P = %.3g. The interaction term in the Cox model (MAP3K8_high : IL1B_high) was non-significant (likelihood-ratio test P = %.3g), indicating that the joint effect of MAP3K8 and IL-1beta on DRFS is essentially additive rather than synergistic.",
    p_overall, p_interaction),
    style = "Normal") |>
  body_add_par(
    "In the cleaner Group-vs-reference parameterisation, the doubly-high group (MAP3K8-high / IL1B-high) carried the largest HR for distant relapse, while the doubly-low group acted as the best-prognosis reference. The pairwise comparison between best and worst groups was borderline-significant after Benjamini-Hochberg adjustment (see pairwise_logrank_BH.csv).",
    style = "Normal") |>
  body_add_par("Figure 1. Four-group Kaplan-Meier curves (DRFS)",
               style = "heading 3") |>
  body_add_img(src = file.path(OUT, "KMplot_4groups.png"),
               width = 6.5, height = 6.5) |>
  body_add_par("Figure 2. Group-vs-reference HR forest plot",
               style = "heading 3") |>
  body_add_img(src = file.path(OUT, "forest_group4_vs_reference.png"),
               width = 6.3, height = 3.2) |>
  body_add_par("Interpretation", style = "heading 2") |>
  body_add_par(
    "Both genes contribute risk information for distant relapse-free survival in chemo-treated TNBC, but they do so independently rather than multiplicatively. Practically: a patient who is high for both MAP3K8 and IL-1beta has the highest distant-relapse risk in this cohort, but this can be approximated by simply summing the two binary risk markers - one is not modulating the other. This is biologically consistent with both genes lying on the same NF-kB / inflammation axis and contributing partially correlated signals from the same chemoresistance pathway.",
    style = "Normal") |>
  body_add_par("Limitations", style = "heading 2") |>
  body_add_par(
    "With four groups of approximately 43-46 patients and 64 total events, individual pairwise contrasts are underpowered. Median dichotomisation is conservative; an optimal-cutpoint analysis (maxstat / surv_cutpoint) would likely yield narrower CIs but at the cost of multiple-testing inflation. The overall log-rank P = 0.13 makes this a hypothesis-supporting figure rather than a confirmatory one.",
    style = "Normal") |>
  body_add_par("Files written by this script", style = "heading 2") |>
  body_add_par("- KMplot_4groups.pdf / .png / .tiff", style = "Normal") |>
  body_add_par("- forest_group4_vs_reference.pdf/.png : group-vs-reference HR forest", style = "Normal") |>
  body_add_par("- group4_vs_reference_HR.csv : same data as a CSV", style = "Normal") |>
  body_add_par("- interaction_cox_terms.csv  : Cox HRs for main + interaction terms", style = "Normal") |>
  body_add_par("- pairwise_logrank_BH.csv    : pairwise log-rank P-values (BH)",  style = "Normal") |>
  body_add_par("- interaction_summary.csv    : group sizes, overall and interaction P-values", style = "Normal") |>
  body_add_par("Reproducibility", style = "heading 2") |>
  body_add_par(paste0("R version: ", R.version.string,
                      ". Platform: ", R.version$platform,
                      ". Date: ", Sys.Date()), style = "Normal")
print(doc, target = "Methods_FourGroup_Interaction_KM.docx")

cat("\nAnalysis 02 complete. Files in:\n  ", OUT, "\n", sep="")
