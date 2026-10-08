# ============================================================================
# FOREST PLOT — Logistic Regression OR (computed from the fitted model)
# ============================================================================

library(ggplot2)
library(dplyr)
library(ggpubr)

# FIX: the original script had the OR/CI/p values typed in by hand as
# literal numbers, disconnected from `logit_model`
# (scripts/11_tcga_validation/risk_table.R). Any re-fit of the model (new
# data, different covariates) would leave this plot silently showing stale
# numbers. Pulled directly from the model instead.

or_ci <- exp(cbind(OR = coef(logit_model), confint(logit_model)))
pvals <- summary(logit_model)$coefficients[, "Pr(>|z|)"]

# Drop the intercept row — not a predictor to plot
or_ci <- or_ci[rownames(or_ci) != "(Intercept)", , drop = FALSE]
pvals <- pvals[names(pvals) != "(Intercept)"]

variable_labels <- c(
  Stromal_High = "High Stromal Score",
  Age          = "Age",
  Stage        = "AJCC Stage (non-I)"
)

format_p_sci <- function(p) {
  # e.g. 3.85e-09 -> "3.85×10⁻⁹"
  exp10    <- floor(log10(p))
  mantissa <- p / 10^exp10
  superscripts <- c("0"="\u2070","1"="\u00b9","2"="\u00b2","3"="\u00b3","4"="\u2074",
                    "5"="\u2075","6"="\u2076","7"="\u2077","8"="\u2078","9"="\u2079",
                    "-"="\u207b")
  exp_chars <- strsplit(as.character(exp10), "")[[1]]
  exp_str   <- paste(superscripts[exp_chars], collapse = "")
  sprintf("%.2f\u00d710%s", mantissa, exp_str)
}

forest_df <- data.frame(
  Variable = unname(variable_labels[rownames(or_ci)]),
  OR       = or_ci[, "OR"],
  Lower    = or_ci[, "2.5 %"],
  Upper    = or_ci[, "97.5 %"],
  P        = sapply(pvals, format_p_sci),
  row.names = NULL
)

forest_df$Variable <- factor(forest_df$Variable,
                             levels = c("AJCC Stage (non-I)", "Age",
                                        "High Stromal Score"))

forest_df$Significant <- ifelse(forest_df$Lower > 1 | forest_df$Upper < 1,
                                "Yes", "No")

p_forest <- ggplot(forest_df, aes(x = OR, y = Variable, color = Significant)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.8) +
  geom_errorbarh(aes(xmin = Lower, xmax = Upper),
                 height = 0.2, linewidth = 0.9) +
  geom_point(size = 4) +
  geom_text(aes(label = paste0("OR = ", round(OR, 2),
                               "\n95% CI: ", round(Lower, 2),
                               "\u2013", round(Upper, 2),
                               "\np = ", P)),
            hjust = -0.1, size = 3.2, color = "black") +
  scale_color_manual(values = c("Yes" = "#E74C3C", "No" = "#95A5A6")) +
  scale_x_log10(limits = c(0.8, 120),
                breaks = c(1, 2, 5, 10, 20, 50),
                labels = c("1", "2", "5", "10", "20", "50")) +
  labs(title = "Independent Predictors of Non-Low Clinical Risk",
       subtitle = paste0("Multivariable logistic regression — TCGA-THCA (n = ",
                         nrow(logit_model$model), ")"),
       x = "Odds Ratio (log scale)", y = NULL) +
  theme_pubr(base_size = 13) +
  theme(legend.position = "none",
        plot.title    = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, color = "grey50"),
        axis.text.y   = element_text(size = 12, face = "bold"))

ggsave("figures/Forest_Plot_Risk_Predictors.png", p_forest,
       width = 8, height = 4, dpi = 300)
cat("Forest plot saved.\n")

# ============================================================================
# EXPRESSION HEATMAP — Top Hub Genes по TME субгрупи (Discovery Cohort)
# ============================================================================

library(pheatmap)
library(tidyverse)

# Топ hub гени от steelblue модула
hub_genes <- c("CTHRC1", "INHBA", "COL3A1", "LAMP5", "LINC01614",
               "COL1A1", "COL1A2", "POSTN", "THBS2", "VCAN",
               "LOXL1", "LOXL2", "THBS1", "ASPN", "LUM",
               "SRPX2", "ANXA8L1")

# Uses `vsd` from the 19-sample discovery-cohort DESeq2 object (see
# scripts/06_differential_expression/pca_volcano.R or
# scripts/07_ecm_stromal/ecm.R — both create `vsd <- vst(dds, blind = FALSE)`
# from the full, non-WGCNA-outlier-excluded `dds`). Do NOT use the WGCNA
# `dds_wgcna`/`norm.counts` object here — it has NG74 excluded (18 samples)
# and would misalign against `colData_TME4` (19 samples) used below.
hub_present <- intersect(hub_genes, rownames(vsd))
cat("Hub гени налични:", length(hub_present), "от", length(hub_genes), "\n")

expr_hub <- assay(vsd)[hub_present, ]

expr_scaled <- t(scale(t(expr_hub)))
expr_scaled <- expr_scaled[complete.cases(expr_scaled), ]

# Metadata от colData, за анотация на колоните по TME група
col_annotation <- data.frame(
  TME_group = colData_TME4$TME_Class2,
  row.names = colnames(expr_hub)
)
table(col_annotation$TME_group)

ann_colors <- list(
  TME_group = c(
    "NG Immune-High"     = "forestgreen",
    "NG Immune-Low"      = "darkolivegreen3",
    "TC Immune-Enriched" = "firebrick",
    "TC Immune-Depleted" = "steelblue"
  )
)

# Нареди колоните по TME група
col_order <- order(col_annotation$TME_group)

# FIX: the original script had TWO pheatmap() calls to the same output
# file. The first one ran BEFORE `col_order` was defined (it was only
# created between the two calls) and would error immediately on a clean
# run. The second call also passed `annotation_colors = colors4` — a flat
# named color vector — where pheatmap requires a *list* keyed by the
# annotation column name (i.e. `ann_colors`, defined above); `colors4` on
# its own is the wrong shape for this argument. Both issues are fixed below
# in a single, working call.
pheatmap(
  expr_scaled[, col_order],
  annotation_col    = col_annotation,
  annotation_colors = ann_colors,
  cluster_cols      = FALSE,
  cluster_rows      = TRUE,
  show_colnames     = TRUE,
  show_rownames     = TRUE,
  fontsize_row      = 10,
  fontsize_col      = 8,
  color             = colorRampPalette(c("#2980B9", "white", "#E74C3C"))(100),
  breaks            = seq(-2, 2, length.out = 101),
  border_color      = NA,
  main              = "Hub Gene Expression — TME Subgroups",
  filename          = "figures/Hub_Genes_Heatmap_Discovery.png",
  width             = 8,
  height            = 7
)

cat("Heatmap saved.\n")
