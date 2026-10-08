library(ggplot2)
library(ggpubr)
library(rstatix)
library(tidyr)
library(dplyr)

# tme_df from ESTIMATE + quanTIseq data
est_df <- as.data.frame(res_est) %>%
  filter(cell_type %in% c("immune score", "stroma score", "estimate score", "tumor purity")) %>%
  pivot_longer(-cell_type, names_to = "Sample", values_to = "Score") %>%
  pivot_wider(names_from = cell_type, values_from = Score)

quanti_df <- as.data.frame(res_quanti) %>%
  filter(cell_type %in% c("T cell CD8+", "T cell regulatory (Tregs)",
                          "Macrophage M1", "Macrophage M2")) %>%
  pivot_longer(-cell_type, names_to = "Sample", values_to = "Fraction") %>%
  pivot_wider(names_from = cell_type, values_from = Fraction)

tme_df <- left_join(est_df, quanti_df, by = "Sample") %>%
  mutate(
    Sample_Type = case_when(
      startsWith(Sample, "TC") ~ "TC",
      startsWith(Sample, "NG") ~ "NG",
      TRUE ~ "Unknown"
    )
  )

# TME classification - 4 groups
# Thresholds were established separately for the malignant and non-malignant groups.

tc_median <- tme_df %>% filter(Sample_Type == "TC") %>% pull(`immune score`) %>% median()
ng_median <- tme_df %>% filter(Sample_Type == "NG") %>% pull(`immune score`) %>% median()

tme_df <- tme_df %>%
  mutate(
    TME_Class2 = case_when(
      Sample_Type == "NG" & `immune score` > ng_median ~ "NG Immune-High",
      Sample_Type == "NG" & `immune score` <= ng_median ~ "NG Immune-Low",
      Sample_Type == "TC" & `immune score` > tc_median ~ "TC Immune-Enriched",
      Sample_Type == "TC" & `immune score` <= tc_median ~ "TC Immune-Depleted",
      TRUE ~ "Other"
    )
  )

cat("TME groups:\n")
print(table(tme_df$TME_Class2))
print(tme_df %>% select(Sample, TME_Class2, `immune score`, `stroma score`))

write.csv(as.data.frame(tme_df), "results/tme_df_TME4.csv", row.names = FALSE)

colors4 <- c(
  "NG Immune-High"     = "forestgreen",
  "NG Immune-Low"      = "darkolivegreen3",
  "TC Immune-Enriched" = "firebrick",
  "TC Immune-Depleted" = "steelblue"
)

# scatter plots 
ggplot(tme_df, aes(x = `stroma score`, y = `immune score`,
                   color = TME_Class2, label = Sample)) +
  geom_point(size = 4) +
  ggrepel::geom_text_repel(size = 3) +
  scale_color_manual(values = colors4) +
  theme_bw() +
  labs(title = "TME Classification — Thyroid Carcinoma",
       x = "Stromal Score (ESTIMATE)", y = "Immune Score (ESTIMATE)",
       color = "TME Class")
ggsave("figures/TME_scatter_4groups_v2.png", width = 8, height = 6, dpi = 300)

ggplot(tme_df, aes(x = `tumor purity`, y = `immune score`,
                   color = TME_Class2, label = Sample)) +
  geom_point(size = 4) +
  ggrepel::geom_text_repel(size = 3) +
  scale_color_manual(values = colors4) +
  theme_bw() +
  labs(title = "Tumor Purity vs Immune Score — Thyroid Carcinoma",
       x = "Tumor Purity (ESTIMATE)", y = "Immune Score (ESTIMATE)",
       color = "TME Class")
ggsave("figures/purity_vs_immune_4groups_v2.png", width = 8, height = 6, dpi = 300)

# quanTIseq long format
quanti_long <- as.data.frame(res_quanti) %>%
  pivot_longer(-cell_type, names_to = "Sample", values_to = "Fraction") %>%
  left_join(tme_df %>% select(Sample, TME_Class2), by = "Sample")

cat("Distribution per group:\n")
print(table(quanti_long$TME_Class2))

# Kruskal-Wallis 
stats_results <- quanti_long %>%
  group_by(cell_type) %>%
  kruskal_test(Fraction ~ TME_Class2) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()

print(stats_results %>% select(cell_type, p, p.adj, p.adj.signif))
write.csv(stats_results, "results/kruskal_wallis_quantiseq.csv", row.names = FALSE)

# Dunn post-hoc 
dunn_results <- quanti_long %>%
  group_by(cell_type) %>%
  dunn_test(Fraction ~ TME_Class2, p.adjust.method = "BH") %>%
  add_significance()
write.csv(dunn_results, "results/dunn_test_quantiseq.csv", row.names = FALSE)

# significant cell types (p.adj)
sig_cells <- stats_results %>% filter(p.adj < 0.05) %>% pull(cell_type)
cat("Significant cell types:\n")
print(sig_cells)

y_max_df <- quanti_long %>%
  filter(cell_type %in% sig_cells) %>%
  group_by(cell_type) %>%
  summarise(y_max = max(Fraction, na.rm = TRUE))

dunn_sig <- dunn_results %>%
  filter(cell_type %in% sig_cells, p.adj < 0.05) %>%
  left_join(y_max_df, by = "cell_type") %>%
  group_by(cell_type) %>%
  mutate(y.position = y_max * (1.1 + (row_number() - 1) * 0.15)) %>%
  ungroup()

quanti_long %>%
  filter(cell_type %in% sig_cells) %>%
  ggplot(aes(x = TME_Class2, y = Fraction, fill = TME_Class2)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.7) +
  facet_wrap(~cell_type, scales = "free_y") +
  scale_fill_manual(values = colors4) +
  stat_pvalue_manual(
    dunn_sig,
    label = "p.adj.signif",
    tip.length = 0.01
  ) +
  theme_bw() +
  theme(axis.text.x = element_blank()) +
  labs(title = "Immune Cell Fractions by TME Group (quanTIseq)",
       x = "", y = "Cell Fraction", fill = "TME Class")

ggsave("figures/boxplot_quantiseq_TME_groups_padj.png", width = 12, height = 8, dpi = 300)
