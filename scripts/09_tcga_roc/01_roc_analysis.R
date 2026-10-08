#ROC ANALYSIS — Stromal Score & TDS discriminative power (TCGA-THCA)

library(pROC)
library(dplyr)
library(ggplot2)

# Stromal Score: BRAF-like vs RAS-like 
roc_df_braf <- scores_tcga_annotated %>%
  mutate(Subtype = case_when(
    paper_BRAFV600E_RAS == "BRAF_V600E" ~ 1,
    paper_BRAFV600E_RAS == "RAS"        ~ 0,
    TRUE ~ NA_real_
  )) %>%
  filter(!is.na(Subtype), !is.na(Stromal_Score))

roc_braf <- roc(roc_df_braf$Subtype, roc_df_braf$Stromal_Score ci = TRUE)  
print(roc_braf)
cat("AUC BRAF vs RAS:", round(auc(roc_braf), 3),
    "95% CI:", round(ci(roc_braf)[1], 3), "-", round(ci(roc_braf)[3], 3), "\n")

# Stromal Score: Tall Cell vs Classical 
roc_df_hist <- scores_tcga_annotated %>%
  mutate(Histo = case_when(
    paper_histological_type == "Tall Cell (>= 50% tall cell features)" ~ 1,
    paper_histological_type == "Classical/usual"                        ~ 0,
    TRUE ~ NA_real_
  )) %>%
  filter(!is.na(Histo), !is.na(Stromal_Score))

roc_hist <- roc(roc_df_hist$Histo, roc_df_hist$Stromal_Score, ci = TRUE)  
print(roc_hist)
cat("AUC Tall Cell vs Classical:", round(auc(roc_hist), 3),
    "95% CI:", round(ci(roc_hist)[1], 3), "-", round(ci(roc_hist)[3], 3), "\n")

# ── 3. TDS: BRAF-like vs RAS-like ────────────────────────────────────────────
roc_df_tds_braf <- scores_tcga_annotated %>%
  mutate(Subtype = case_when(
    paper_BRAFV600E_RAS == "BRAF_V600E" ~ 1,
    paper_BRAFV600E_RAS == "RAS"        ~ 0,
    TRUE ~ NA_real_
  )) %>%
  filter(!is.na(Subtype), !is.na(TDS))

roc_tds_braf <- roc(roc_df_tds_braf$Subtype, roc_df_tds_braf$TDS, ci = TRUE)  
cat("AUC TDS BRAF vs RAS:", round(auc(roc_tds_braf), 3),
    "95% CI:", round(ci(roc_tds_braf)[1], 3), "-", round(ci(roc_tds_braf)[3], 3), "\n")

# Stromal Score: Follicular vs Classical 
roc_df_foll <- scores_tcga_annotated %>%
  mutate(Histo2 = case_when(
    paper_histological_type == "Follicular (>= 99% follicular patterned)" ~ 1,
    paper_histological_type == "Classical/usual"                           ~ 0,
    TRUE ~ NA_real_
  )) %>%
  filter(!is.na(Histo2), !is.na(Stromal_Score))

roc_foll <- roc(roc_df_foll$Histo2, roc_df_foll$Stromal_Score, ci = TRUE) 
cat("AUC Stromal Follicular vs Classical:", round(auc(roc_foll), 3),
    "95% CI:", round(ci(roc_foll)[1], 3), "-", round(ci(roc_foll)[3], 3), "\n")

# Combined plot — всички 4 ROC криви 
make_roc_df <- function(roc_obj, label) {
  data.frame(
    specificity = roc_obj$specificities,
    sensitivity = roc_obj$sensitivities,
    Comparison  = paste0(label, "\nAUC = ", round(auc(roc_obj), 3),
                         " (", round(ci(roc_obj)[1], 3),
                         "–", round(ci(roc_obj)[3], 3), ")")
  )
}

roc_combined <- bind_rows(
  make_roc_df(roc_braf,     "Stromal: BRAF-like vs RAS-like"),
  make_roc_df(roc_hist,     "Stromal: Tall Cell vs Classical"),
  make_roc_df(roc_tds_braf, "TDS: BRAF-like vs RAS-like"),
  make_roc_df(roc_foll,     "Stromal: Follicular vs Classical")
)

roc_colors <- setNames(
  c("#E74C3C", "#8E44AD", "#E67E22", "#27AE60"),
  unique(roc_combined$Comparison)
)

p_roc_combined <- ggplot(roc_combined, aes(x = 1 - specificity, y = sensitivity,
                                           color = Comparison)) +
  geom_line(linewidth = 1.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "grey60", linewidth = 0.8) +
  scale_color_manual(values = roc_colors) +
  scale_x_continuous(limits = c(0, 1), expand = c(0.01, 0.01)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0.01, 0.01)) +
  labs(title = "ROC Analysis — Stromal Signature & TDS Discriminative Power",
       subtitle = "TCGA-THCA cohort",
       x = "1 − Specificity (False Positive Rate)",
       y = "Sensitivity (True Positive Rate)",
       color = NULL) +
  theme_pubr(base_size = 13) +
  theme(plot.title      = element_text(hjust = 0.5),
        plot.subtitle   = element_text(hjust = 0.5, color = "grey50"),
        legend.position = c(0.60, 0.25),
        legend.background = element_rect(fill = "white", color = "grey80"),
        legend.text     = element_text(size = 9))

ggsave("figures/TCGA_ROC_Combined.png", p_roc_combined, width = 7, height = 6, dpi = 300)
