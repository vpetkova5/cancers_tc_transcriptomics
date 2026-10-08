# CAF analysis

# descriptive_stats
caf_stats <- epic_long %>%
  filter(cell_type == "Cancer associated fibroblast") %>%
  group_by(TME_Class2) %>%
  summarise(
    n = n(),
    mean_pct = mean(Fraction) * 100,
    sd_pct = sd(Fraction) * 100,
    median_pct = median(Fraction) * 100,
    .groups = "drop"
  )
print(caf_stats)

# Kruskal-Wallis + Dunn 
kw_caf <- epic_long %>%
  filter(cell_type == "Cancer associated fibroblast") %>%
  kruskal_test(Fraction ~ TME_Class2)
print(kw_caf)

dunn_caf <- epic_long %>%
  filter(cell_type == "Cancer associated fibroblast") %>%
  dunn_test(Fraction ~ TME_Class2, p.adjust.method = "BH") %>%
  add_significance()
print(dunn_caf)

# correlation CAF vs ESTIMATE scores 
caf_cor_df <- epic_long %>%
  filter(cell_type == "Cancer associated fibroblast") %>%
  left_join(tme_df %>% select(Sample, `immune score`,
                              `stroma score`, `tumor purity`),
            by = "Sample")

cor_caf_immune  <- cor.test(caf_cor_df$Fraction, caf_cor_df$`immune score`,
                            method = "spearman")
cor_caf_stromal <- cor.test(caf_cor_df$Fraction, caf_cor_df$`stroma score`,
                            method = "spearman")
cor_caf_purity  <- cor.test(caf_cor_df$Fraction, caf_cor_df$`tumor purity`,
                            method = "spearman")

cat("CAF vs Immune Score:  rho =", round(cor_caf_immune$estimate, 3),
    " p =", round(cor_caf_immune$p.value, 4), "\n")
cat("CAF vs Stromal Score: rho =", round(cor_caf_stromal$estimate, 3),
    " p =", round(cor_caf_stromal$p.value, 4), "\n")
cat("CAF vs Tumor Purity:  rho =", round(cor_caf_purity$estimate, 3),
    " p =", round(cor_caf_purity$p.value, 4), "\n")

# CAF vs MCP-counter Fibroblast (validation)
mcp_long <- as.data.frame(res_mcp) %>%
  pivot_longer(-cell_type, names_to = "Sample", values_to = "Score") %>%
  left_join(tme_df %>% select(Sample, TME_Class2), by = "Sample")

caf_mcp <- mcp_long %>%
  filter(grepl("ibroblast|troma", cell_type, ignore.case = TRUE))

if (nrow(caf_mcp) > 0) {
  caf_validation <- epic_long %>%
    filter(cell_type == "Cancer associated fibroblast") %>%
    select(Sample, CAF_EPIC = Fraction) %>%
    left_join(
      caf_mcp %>% select(Sample, MCP_score = Score),
      by = "Sample"
    )

  cor_validation <- cor.test(caf_validation$CAF_EPIC,
                             caf_validation$MCP_score,
                             method = "spearman")
  cat("CAF EPIC vs MCP-counter: rho =", round(cor_validation$estimate, 3),
      " p =", round(cor_validation$p.value, 4), "\n")
}

# CAF vs Endothelial (stroma)
epic_stromal <- epic_long %>%
  filter(cell_type %in% c("Cancer associated fibroblast",
                          "Endothelial cell")) %>%
  group_by(TME_Class2, cell_type) %>%
  summarise(mean_pct = mean(Fraction) * 100,
            sd_pct = sd(Fraction) * 100, .groups = "drop")
print(epic_stromal)

#  Plots 

# Boxplot CAF
y_max_caf <- epic_long %>%
  filter(cell_type == "Cancer associated fibroblast") %>%
  summarise(y = quantile(Fraction, 0.95)) %>% pull(y)

# p.adj.signif
dunn_caf_plot <- dunn_caf %>%
  filter(p.adj < 0.05) %>%
  mutate(y.position = y_max_caf * (1.1 + (row_number() - 1) * 0.18) * 100)

pA_caf <- epic_long %>%
  filter(cell_type == "Cancer associated fibroblast") %>%
  ggplot(aes(x = TME_Class2, y = Fraction * 100, fill = TME_Class2)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 2.5, alpha = 0.7) +
  scale_fill_manual(values = colors4) +
  stat_pvalue_manual(dunn_caf_plot, label = "p.adj.signif",
                     tip.length = 0.01, size = 4) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_blank()) +
  labs(title = "A  CAF Infiltration by TME Subgroup",
       x = "", y = "CAF Fraction (%)", fill = "TME Subgroup")

# Scatter CAF vs Stromal Score
pB_caf <- ggplot(caf_cor_df, aes(x = `stroma score`, y = Fraction * 100,
                                 color = TME_Class2, label = Sample)) +
  geom_point(size = 4) +
  geom_smooth(method = "lm", se = TRUE, color = "gray40", linetype = "dashed") +
  ggrepel::geom_text_repel(size = 3) +
  scale_color_manual(values = colors4) +
  theme_bw(base_size = 12) +
  labs(title = "B  CAF vs Stromal Score",
       x = "Stromal Score (ESTIMATE)", y = "CAF Fraction (%)",
       color = "TME Subgroup",
       subtitle = paste0("Spearman ρ = ", round(cor_caf_stromal$estimate, 3),
                         ", p = ", format(cor_caf_stromal$p.value, digits = 3)))

# Scatter CAF vs Immune Score
pC_caf <- ggplot(caf_cor_df, aes(x = `immune score`, y = Fraction * 100,
                                 color = TME_Class2, label = Sample)) +
  geom_point(size = 4) +
  geom_smooth(method = "lm", se = TRUE, color = "gray40", linetype = "dashed") +
  ggrepel::geom_text_repel(size = 3) +
  scale_color_manual(values = colors4) +
  theme_bw(base_size = 12) +
  labs(title = "C  CAF vs Immune Score",
       x = "Immune Score (ESTIMATE)", y = "CAF Fraction (%)",
       color = "TME Subgroup",
       subtitle = paste0("Spearman ρ = ", round(cor_caf_immune$estimate, 3),
                         ", p = ", format(cor_caf_immune$p.value, digits = 3)))

# Stromal composition — CAF vs Endothelial stacked
pD_caf <- epic_long %>%
  filter(cell_type %in% c("Cancer associated fibroblast",
                          "Endothelial cell")) %>%
  mutate(cell_type = recode(cell_type,
                            "Cancer associated fibroblast" = "CAF",
                            "Endothelial cell" = "Endothelial")) %>%
  ggplot(aes(x = TME_Class2, y = Fraction * 100, fill = cell_type)) +
  geom_boxplot(position = position_dodge(0.8), outlier.shape = NA) +
  geom_jitter(position = position_jitterdodge(0.8), size = 2, alpha = 0.6) +
  scale_fill_manual(values = c("CAF" = "darkorange",
                               "Endothelial" = "mediumpurple")) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_blank()) +
  facet_wrap(~TME_Class2, scales = "free_x", nrow = 1) +
  labs(title = "D  Stromal Composition by TME Subgroup",
       x = "", y = "Cell Fraction (%)", fill = "Cell Type")

# Combine Figure 3
figure3 <- (pA_caf | pB_caf | pC_caf) / pD_caf +
  plot_annotation(
    title = "Figure 3. Cancer-Associated Fibroblast Infiltration in Thyroid Carcinoma TME",
    theme = theme(plot.title = element_text(size = 13, face = "bold"))
  ) +
  plot_layout(guides = "collect")

ggsave("figures/Figure3_CAF_analysis.png", figure3, width = 16, height = 10, dpi = 300)
ggsave("figures/Figure3_CAF_analysis.pdf", figure3, width = 16, height = 10)
ggsave("figures/pA_caf.png", pA_caf, width = 8, height = 6, dpi = 300)
ggsave("figures/pB_caf.png", pB_caf, width = 8, height = 6, dpi = 300)
ggsave("figures/pC_caf.png", pC_caf, width = 8, height = 6, dpi = 300)
ggsave("figures/pD_caf.png", pD_caf, width = 6, height = 4, dpi = 300)

library(Hmisc)

