# TCGA VALIDATION — STROMAL SCORE (tc_unique gene set)
# Computes the ssGSEA "Stromal_Score" for the TCGA-THCA
# cohort using the discovery-cohort tc_unique gene set, and validates
# it against the paper-reported Thyroid Differentiation Score (TDS).

library(tidyverse)
library(ggpubr)
library(rstatix)

# already saved  TCGA data - only pappilary thyroid cancer, after combat batch correction
vst_matrix_tcga <- readRDS("data/TCGA_PTC_vst_combat.rds")
metadata_ptc    <- readRDS("data/TCGA_PTC_metadata.rds")

# Gene sets 
# tc_unique = discovery-cohort stromal gene set
gs_list <- list(Stromal = tc_unique)

tds_genes <- c("DIO1", "DIO2", "DUOX1", "DUOX2", "FOXE1", "GLIS3", "NKX2-1",
               "PAX8", "SLC26A4", "SLC5A5", "SLC5A8", "TG", "THRA", "THRB", "TPO", "TSHR")

genes_in_tcga <- intersect(tc_unique, rownames(vst_matrix_tcga))
cat("Stromal genes в TCGA VST:", length(genes_in_tcga), "от", length(tc_unique), "\n")

# Functions - calculate_tds; run_ssGSEA_native
calculate_tds <- function(matrix, genes) {
  genes_present <- intersect(genes, rownames(matrix))
  gene_expr     <- matrix[genes_present, , drop = FALSE]
  gene_medians  <- apply(gene_expr, 1, median, na.rm = TRUE)
  deviations    <- sweep(gene_expr, 1, gene_medians, "-")
  return(colMeans(deviations, na.rm = TRUE))
}

run_ssGSEA_native <- function(X, gene_sets) {
  R        <- apply(X, 2, rank, ties.method = "average")
  num_genes <- nrow(X)
  es_results <- matrix(NA, nrow = length(gene_sets), ncol = ncol(X),
                       dimnames = list(names(gene_sets), colnames(X)))
  for (set_name in names(gene_sets)) {
    genes_in_set <- intersect(gene_sets[[set_name]], rownames(X))
    n_h <- length(genes_in_set)
    if (n_h > 0) {
      for (j in 1:ncol(X)) {
        sorted_indices <- order(R[, j], decreasing = TRUE)
        is_hit <- rownames(X)[sorted_indices] %in% genes_in_set
        walk   <- cumsum((is_hit / n_h) - ((!is_hit) / (num_genes - n_h)))
        es_results[set_name, j] <- sum(walk)
      }
    }
  }
  return(as.data.frame(t(es_results)))
}

# compute scores 
TDS_tcga    <- calculate_tds(vst_matrix_tcga, tds_genes)
ssgsea_tcga <- run_ssGSEA_native(as.matrix(vst_matrix_tcga), gs_list)

scores_tcga <- data.frame(
  SampleID      = rownames(ssgsea_tcga),
  TDS           = TDS_tcga[rownames(ssgsea_tcga)],
  Stromal_Score = ssgsea_tcga$Stromal
) %>%
  mutate(barcode = substr(gsub("\\.", "-", SampleID), 1, 15))

# metadata merge 
metadata_tcga <- metadata_ptc %>%
  mutate(barcode = substr(gsub("\\.", "-", barcode), 1, 15))

scores_tcga_annotated <- scores_tcga %>%
  left_join(metadata_tcga, by = "barcode")

match_count <- sum(!is.na(scores_tcga_annotated$patient))
cat("Matched samples:", match_count, "от", nrow(scores_tcga), "\n")

# TDS validation vs paper-reported values
cor_test_tds <- cor.test(scores_tcga_annotated$TDS,
                         scores_tcga_annotated$TDS_paper,
                         method = "pearson", use = "complete.obs")
print(cor_test_tds)

p_tds_val <- ggplot(scores_tcga_annotated, aes(x = TDS_paper, y = TDS)) +
  geom_point(alpha = 0.3, color = "midnightblue") +
  geom_smooth(method = "lm", color = "red") +
  stat_cor(method = "pearson") +
  labs(title = "TDS Validation: Calculated vs. Official Values",
       x = "Official TDS (Paper)", y = "Calculated TDS") +
  theme_pubr()
ggsave("figures/TDS_validation_441.png", p_tds_val, width = 6, height = 4, dpi = 300)

# phenotypic validation — BRAF vs RAS 
plot_df <- scores_tcga_annotated %>%
  mutate(Subtype = case_when(
    paper_BRAFV600E_RAS == "BRAF_V600E" ~ "BRAF-like",
    paper_BRAFV600E_RAS == "RAS"        ~ "RAS-like",
    TRUE                                ~ "Other/WT"
  )) %>%
  filter(!is.na(TDS), !is.na(Stromal_Score))

plot_df$Subtype <- factor(plot_df$Subtype, levels = c("RAS-like", "Other/WT", "BRAF-like"))

gradient_colors <- c("RAS-like" = "#5DADE2", "Other/WT" = "#ABB2B9", "BRAF-like" = "#E74C3C")

create_validation_gradient <- function(data, variable, y_label, title, colors) {
  n_labels <- data %>%
    group_by(Subtype) %>% tally() %>%
    mutate(label = paste0(Subtype, "\n(n=", n, ")"))
  label_map <- setNames(n_labels$label, n_labels$Subtype)

  ggplot(data, aes(x = Subtype, y = .data[[variable]], fill = Subtype)) +
    geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.8, color = "#2c3e50") +
    geom_jitter(width = 0.2, alpha = 0.2, size = 1.2, color = "black") +
    stat_summary(geom = "text", fun = median,
                 aes(label = sprintf("%.0f", after_stat(y))),
                 hjust = -0.15, vjust = -0.5,
                 fontface = "bold", size = 3.5) +
    stat_compare_means(label.y = max(data[[variable]], na.rm = TRUE) * 1.1) +
    stat_compare_means(comparisons = list(c("RAS-like", "BRAF-like"),
                                          c("RAS-like", "Other/WT"),
                                          c("Other/WT", "BRAF-like")),
                       label = "p.signif") +
    scale_fill_manual(values = colors) +
    scale_x_discrete(labels = label_map) +
    labs(title = title, y = y_label, x = NULL) +
    theme_pubr(base_size = 14) +
    theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
}

p_stromal <- create_validation_gradient(
  plot_df, "Stromal_Score", "Stromal ssGSEA Score (441 genes)",
  "Stromal Infiltration Gradient — TCGA Validation", gradient_colors
)
ggsave("figures/TCGA_Stromal_Gradient_441genes.png", p_stromal, width = 7, height = 6, dpi = 300)

# export results 
write.csv(scores_tcga_annotated,
          "results/TCGA_Validation_Stromal441_annotated.csv", row.names = FALSE)

# TDS ~ STROMAL SCORE CORRELATION — TCGA VALIDATION

cor_test_stromal <- cor.test(scores_tcga_annotated$TDS,
                             scores_tcga_annotated$Stromal_Score,
                             method = "pearson", use = "complete.obs")
print(cor_test_stromal)

p_corr <- ggplot(scores_tcga_annotated, aes(x = TDS, y = Stromal_Score)) +
  geom_point(alpha = 0.35, size = 1.8, color = "#2C3E50") +
  geom_smooth(method = "lm", color = "#E74C3C", linewidth = 1) +
  stat_cor(method = "pearson",
           label.x = 0.5, label.y = 6800,
           size = 4.5, fontface = "italic") +
  labs(
    title = "Thyroid Differentiation vs. Stromal Infiltration",
    subtitle = "TCGA-THCA cohort (n = 456)",
    x = "Thyroid Differentiation Score (TDS)",
    y = "Stromal ssGSEA Score (441 genes)"
  ) +
  theme_pubr(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, color = "grey50"))

ggsave("figures/TCGA_TDS_vs_Stromal_correlation.png", p_corr, width = 6, height = 5, dpi = 300)


