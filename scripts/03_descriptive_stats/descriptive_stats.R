# descriptive_stats for the 4 TME groups
tme_stats <- tme_df %>%
  group_by(TME_Class2) %>%
  summarise(
    n = n(),
    immune_mean = mean(`immune score`),
    immune_sd = sd(`immune score`),
    immune_median = median(`immune score`),
    stromal_mean = mean(`stroma score`),
    stromal_sd = sd(`stroma score`),
    purity_mean = mean(`tumor purity`),
    purity_sd = sd(`tumor purity`),
    purity_mean_pct = mean(`tumor purity`) * 100,
    purity_sd_pct = sd(`tumor purity`) * 100
  )
print(tme_stats)
write.csv(tme_stats, "results/tme_descriptive_stats.csv", row.names = FALSE)

# Kruskal-Wallis for immune score 
kw_immune <- kruskal_test(tme_df, `immune score` ~ TME_Class2)
print(kw_immune)

# Dunn post-hoc for immune score 
dunn_immune <- dunn_test(tme_df, `immune score` ~ TME_Class2,
                         p.adjust.method = "BH")
print(dunn_immune)

# Kruskal-Wallis for tumor purity 
kw_purity <- kruskal_test(tme_df, `tumor purity` ~ TME_Class2)
print(kw_purity)

dunn_purity <- dunn_test(tme_df, `tumor purity` ~ TME_Class2,
                         p.adjust.method = "BH")
print(dunn_purity)

# correlation immune score vs tumor purity 
cor_test <- cor.test(tme_df$`immune score`, tme_df$`tumor purity`,
                     method = "spearman")
print(cor_test)

#  EPIC data — B cell / CAF / Macrophage per groups
epic_long <- as.data.frame(res_epic) %>%
  pivot_longer(-cell_type, names_to = "Sample", values_to = "Fraction") %>%
  left_join(tme_df %>% select(Sample, TME_Class2), by = "Sample")

epic_stats <- epic_long %>%
  group_by(TME_Class2, cell_type) %>%
  summarise(mean = mean(Fraction), sd = sd(Fraction), .groups = "drop")
print(epic_stats)

kw_epic <- epic_long %>%
  group_by(cell_type) %>%
  kruskal_test(Fraction ~ TME_Class2) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()
print(kw_epic)
print(table(tme_df$TME_Class2))
print(tme_df %>% select(Sample, TME_Class2, `immune score`) %>% arrange(TME_Class2))

# Dunn post-hoc for EPIC
dunn_epic <- epic_long %>%
  group_by(cell_type) %>%
  dunn_test(Fraction ~ TME_Class2, p.adjust.method = "BH") %>%
  add_significance()

# sig cell types adjust_pvalue 
sig_cells_epic <- kw_epic %>% filter(p.adj < 0.05) %>% pull(cell_type)
print(sig_cells_epic)

# descriptive_stats
epic_long %>%
  filter(cell_type %in% sig_cells_epic) %>%
  group_by(TME_Class2, cell_type) %>%
  summarise(
    mean_pct = mean(Fraction) * 100,
    sd_pct = sd(Fraction) * 100,
    .groups = "drop"
  ) %>%
  print()
