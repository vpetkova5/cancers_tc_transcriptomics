# CLINICAL & MOLECULAR VALIDATION PANELS
library(ggpubr)
library(tidyverse)
library(rstatix)

# PANEL 1: CLINICAL

# 1A: Histological type 
hist_df <- scores_tcga_annotated %>%
  mutate(Histology = case_when(
    paper_histological_type == "Classical/usual"                        ~ "Classical",
    paper_histological_type == "Follicular (>= 99% follicular patterned)" ~ "Follicular",
    paper_histological_type == "Tall Cell (>= 50% tall cell features)"  ~ "Tall Cell",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(Histology), !is.na(Stromal_Score))

hist_df$Histology <- factor(hist_df$Histology, levels = c("Follicular", "Classical", "Tall Cell"))

hist_colors <- c("Follicular" = "#5DADE2", "Classical" = "#ABB2B9", "Tall Cell" = "#E74C3C")

p1a <- ggplot(hist_df, aes(x = Histology, y = Stromal_Score, fill = Histology)) +
  geom_boxplot(outlier.shape = NA, width = 0.55, alpha = 0.85, color = "#2c3e50") +
  geom_jitter(width = 0.18, alpha = 0.25, size = 1.2, color = "black") +
  stat_summary(geom = "text", fun = median,
               aes(label = sprintf("%.0f", after_stat(y))),
               hjust = -0.15, vjust = -0.5, fontface = "bold", size = 3.2) +
  stat_compare_means(method = "kruskal.test",
                     label.y = max(hist_df$Stromal_Score, na.rm = TRUE) * 1.35) +
  stat_compare_means(comparisons = list(c("Follicular", "Classical"),
                                        c("Follicular", "Tall Cell"),
                                        c("Classical",  "Tall Cell")),
                     label = "p.signif", step.increase = 0.08) +
  scale_fill_manual(values = hist_colors) +
  scale_x_discrete(labels = function(x) {
    n <- table(hist_df$Histology)
    paste0(x, "\n(n=", n[x], ")")
  }) +
  labs(title = "Histological Subtype", y = "Stromal ssGSEA Score (441 genes)", x = NULL) +
  theme_pubr(base_size = 12) +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

# 1B: AJCC Stage 
stage_df <- scores_tcga_annotated %>%
  mutate(Stage = case_when(
    ajcc_pathologic_stage == "Stage I"   ~ "I",
    ajcc_pathologic_stage == "Stage II"  ~ "II",
    ajcc_pathologic_stage == "Stage III" ~ "III",
    ajcc_pathologic_stage %in% c("Stage IV", "Stage IVA", "Stage IVC") ~ "IV",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(Stage), !is.na(Stromal_Score))

stage_df$Stage <- factor(stage_df$Stage, levels = c("I", "II", "III", "IV"))

stage_colors <- c("I" = "#AED6F1", "II" = "#5DADE2", "III" = "#E59866", "IV" = "#E74C3C")

p1b <- ggplot(stage_df, aes(x = Stage, y = Stromal_Score, fill = Stage)) +
  geom_boxplot(outlier.shape = NA, width = 0.55, alpha = 0.85, color = "#2c3e50") +
  geom_jitter(width = 0.18, alpha = 0.25, size = 1.2, color = "black") +
  stat_summary(geom = "text", fun = median,
               aes(label = sprintf("%.0f", after_stat(y))),
               hjust = -0.15, vjust = -0.5, fontface = "bold", size = 3.2) +
  stat_compare_means(method = "kruskal.test",
                     label.y = max(stage_df$Stromal_Score, na.rm = TRUE) * 1.35) +
  stat_compare_means(comparisons = list(c("I", "III"), c("I", "IV"), c("II", "IV")),
                     label = "p.signif", step.increase = 0.08) +
  scale_fill_manual(values = stage_colors) +
  scale_x_discrete(labels = function(x) {
    n <- table(stage_df$Stage)
    paste0("Stage ", x, "\n(n=", n[x], ")")
  }) +
  labs(title = "AJCC Pathologic Stage", y = NULL, x = NULL) +
  theme_pubr(base_size = 12) +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

# 1C: Extrathyroidal extension 
ext_df <- scores_tcga_annotated %>%
  mutate(Extension = case_when(
    paper_Extrathyroidal_extension == "None"                    ~ "None",
    paper_Extrathyroidal_extension == "Minimal (T3)"            ~ "Minimal",
    paper_Extrathyroidal_extension == "Moderate/Advanced (T4a)" ~ "Moderate/\nAdvanced",
    paper_Extrathyroidal_extension == "Very Advanced (T4b)"     ~ "Moderate/\nAdvanced",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(Extension), !is.na(Stromal_Score))

ext_df$Extension <- factor(ext_df$Extension, levels = c("None", "Minimal", "Moderate/\nAdvanced"))

ext_colors <- c("None" = "#AED6F1", "Minimal" = "#E59866", "Moderate/\nAdvanced" = "#E74C3C")

p1c <- ggplot(ext_df, aes(x = Extension, y = Stromal_Score, fill = Extension)) +
  geom_boxplot(outlier.shape = NA, width = 0.55, alpha = 0.85, color = "#2c3e50") +
  geom_jitter(width = 0.18, alpha = 0.25, size = 1.2, color = "black") +
  stat_summary(geom = "text", fun = median,
               aes(label = sprintf("%.0f", after_stat(y))),
               hjust = -0.15, vjust = -0.5, fontface = "bold", size = 3.2) +
  stat_compare_means(method = "kruskal.test",
                     label.y = max(ext_df$Stromal_Score, na.rm = TRUE) * 1.35) +
  stat_compare_means(comparisons = list(c("None", "Minimal"),
                                        c("None", "Moderate/\nAdvanced"),
                                        c("Minimal", "Moderate/\nAdvanced")),
                     label = "p.signif", step.increase = 0.08) +
  scale_fill_manual(values = ext_colors) +
  scale_x_discrete(labels = function(x) {
    n <- table(ext_df$Extension)
    paste0(x, "\n(n=", n[x], ")")
  }) +
  labs(title = "Extrathyroidal Extension", y = NULL, x = NULL) +
  theme_pubr(base_size = 12) +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

# combine Panel 1 
panel1 <- ggarrange(p1a, p1b, p1c, ncol = 3, align = "hv")
panel1 <- annotate_figure(panel1,
                          top = text_grob("Stromal Signature — Clinical Associations (TCGA-THCA)",
                                          face = "bold", size = 14))

ggsave("figures/TCGA_Clinical_Panel_3up_preview.png", panel1, width = 16, height = 6, dpi = 300)

# PANEL 2: MOLECULAR assosiations

# 2A: mRNA Cluster 
cluster_df <- scores_tcga_annotated %>%
  filter(!is.na(paper_mRNA_Cluster_number), !is.na(Stromal_Score)) %>%
  mutate(Cluster = factor(paper_mRNA_Cluster_number))

# TCGA mRNA cluster labels (от оригиналната TCGA THCA paper)
cluster_labels <- c("1" = "C1\nBRAF-like\nn=124",
                    "2" = "C2\nBRAF-like\nn=61",
                    "3" = "C3\nRAS-like\nn=75",
                    "4" = "C4\nRAS-like\nn=98",
                    "5" = "C5\nMixed\nn=87")

cluster_colors <- c("1" = "#C0392B", "2" = "#E67E22",
                    "3" = "#3498DB", "4" = "#2980B9", "5" = "#ABB2B9")

p2a <- ggplot(cluster_df, aes(x = Cluster, y = Stromal_Score, fill = Cluster)) +
  geom_boxplot(outlier.shape = NA, width = 0.55, alpha = 0.85, color = "#2c3e50") +
  geom_jitter(width = 0.18, alpha = 0.25, size = 1.2, color = "black") +
  stat_summary(geom = "text", fun = median,
               aes(label = sprintf("%.0f", after_stat(y))),
               hjust = -0.15, vjust = -0.5, fontface = "bold", size = 3.2) +
  stat_compare_means(method = "kruskal.test",
                     label.y = max(cluster_df$Stromal_Score, na.rm = TRUE) * 1.08) +
  scale_fill_manual(values = cluster_colors) +
  scale_x_discrete(labels = cluster_labels) +
  labs(title = "mRNA Expression Cluster (TCGA)",
       y = "Stromal ssGSEA Score (441 genes)", x = NULL) +
  theme_pubr(base_size = 12) +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

ggsave("figures/TCGA_Molecular_Panel2_cluster_only_preview.png", p2a, width = 8, height = 6, dpi = 300)

panel2 <- ggarrange(p2a, p_corr, ncol = 2, widths = c(1.3, 1), align = "hv")
panel2 <- annotate_figure(panel2,
                          top = text_grob("Stromal Signature — Molecular Associations (TCGA-THCA)",
                                          face = "bold", size = 14))

ggsave("figures/TCGA_Molecular_Panel2_cluster_and_corr_preview.png", panel2, width = 15, height = 6, dpi = 300)

# TDS categories
tds_df <- scores_tcga_annotated %>%
  filter(!is.na(TDS), !is.na(Stromal_Score)) %>%
  mutate(TDS_group = case_when(
    TDS >= quantile(TDS, 0.67, na.rm = TRUE) ~ "High TDS",
    TDS <= quantile(TDS, 0.33, na.rm = TRUE) ~ "Low TDS",
    TRUE ~ "Mid TDS"
  )) %>%
  mutate(TDS_group = factor(TDS_group, levels = c("Low TDS", "Mid TDS", "High TDS")))

table(tds_df$TDS_group)

tds_colors <- c("Low TDS" = "#E74C3C", "Mid TDS" = "#ABB2B9", "High TDS" = "#5DADE2")

p_tds_box <- ggplot(tds_df, aes(x = TDS_group, y = Stromal_Score, fill = TDS_group)) +
  geom_boxplot(outlier.shape = NA, width = 0.55, alpha = 0.85, color = "#2c3e50") +
  geom_jitter(width = 0.18, alpha = 0.25, size = 1.2, color = "black") +
  stat_summary(geom = "text", fun = median,
               aes(label = sprintf("%.0f", after_stat(y))),
               hjust = -0.15, vjust = -0.5, fontface = "bold", size = 3.5) +
  stat_compare_means(method = "kruskal.test",
                     label.y = max(tds_df$Stromal_Score, na.rm = TRUE) * 1.35) +
  stat_compare_means(comparisons = list(c("Low TDS", "Mid TDS"),
                                        c("Mid TDS", "High TDS"),
                                        c("Low TDS", "High TDS")),
                     label = "p.signif", step.increase = 0.08) +
  scale_fill_manual(values = tds_colors) +
  scale_x_discrete(labels = function(x) {
    n <- table(tds_df$TDS_group)
    paste0(x, "\n(n=", n[x], ")")
  }) +
  labs(title = "Stromal Infiltration by Differentiation Status",
       subtitle = "TCGA-THCA cohort — TDS tertiles",
       x = NULL, y = "Stromal ssGSEA Score (441 genes)") +
  theme_pubr(base_size = 13) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, color = "grey50"))

ggsave("figures/TCGA_Stromal_by_TDS_tertiles.png", p_tds_box, width = 6, height = 5, dpi = 300)

