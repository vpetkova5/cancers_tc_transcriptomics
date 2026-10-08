library(ggplot2)
library(dplyr)
library(patchwork)
library(enrichR)

# TC-unique cancer-intrinsic genes after intersection analysis
genes_tc_unique <- data.frame(V1 = tc_unique)

# Shared inflammatory гени
genes_shared <- data.frame(V1 = shared_inflammatory)

cat("TC-unique genes:", length(tc_unique), "\n")
cat("Shared inflammatory:", length(shared_inflammatory), "\n")


# ENRICHR
dbs <- c("GO_Biological_Process_2025",
         "GO_Molecular_Function_2025",
         "GO_Cellular_Component_2025",
         "Reactome_2022")

tc_unique_enrichment <- enrichr(tc_unique, dbs)
shared_enrichment    <- enrichr(shared_inflammatory, dbs)

# ENRICHR function

prepare_dotplot_data <- function(enrichr_df, category_name, gene_set_name, top_n = 10) {
  if (is.null(enrichr_df) || nrow(enrichr_df) == 0) {
    warning(paste("No data for", category_name))
    return(NULL)
  }
  enrichr_df %>%
    arrange(Adjusted.P.value) %>%
    head(top_n) %>%
    mutate(
      Category     = category_name,
      Gene_Set     = gene_set_name,
      GeneRatio    = sapply(strsplit(as.character(Overlap), "/"), function(x)
        if (length(x) == 2) as.numeric(x[1]) / as.numeric(x[2]) else NA),
      Count        = as.numeric(sapply(strsplit(as.character(Overlap), "/"), function(x)
        if (length(x) >= 1) x[1] else NA)),
      neg_log10_p  = -log10(Adjusted.P.value),
      Term_clean   = gsub(" \\(GO:.*\\)", "", Term),
      Term_clean   = gsub(" R-HSA-.*", "", Term_clean),
      Term_clean   = ifelse(nchar(Term_clean) > 50,
                            paste0(substr(Term_clean, 1, 47), "..."),
                            Term_clean)
    ) %>%
    arrange(desc(neg_log10_p))
}

# TC-unique stromal genes  

stromal_cc <- prepare_dotplot_data(
  tc_unique_enrichment$GO_Cellular_Component_2025,
  "GO Cellular Component",
  paste0("TC Immune-Enriched unique (", length(tc_unique), " genes)"),
  top_n = 10)

stromal_reactome <- prepare_dotplot_data(
  tc_unique_enrichment$Reactome_2022,
  "Reactome Pathways",
  paste0("TC Immune-Enriched unique (", length(tc_unique), " genes)"),
  top_n = 10)

# shared inflamatory 

inflammatory_cc <- prepare_dotplot_data(
  shared_enrichment$GO_Cellular_Component_2025,
  "GO Cellular Component",
  paste0("Shared inflammatory (", length(shared_inflammatory), " genes)"),
  top_n = 10)

inflammatory_reactome <- prepare_dotplot_data(
  shared_enrichment$Reactome_2022,
  "Reactome Pathways",
  paste0("Shared inflammatory (", length(shared_inflammatory), " genes)"),
  top_n = 10)

# plots

p_stromal_cc <- ggplot(stromal_cc,
                       aes(x = GeneRatio, y = reorder(Term_clean, GeneRatio))) +
  geom_point(aes(size = Count, color = neg_log10_p)) +
  scale_color_gradientn(
    colors = c("#FEE5D9", "#FC9272", "#DE2D26", "#A50F15"),
    name   = "-log10\n(adj.p)") +
  scale_size_continuous(name = "Gene\nCount", range = c(3, 10)) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.y        = element_text(size = 9),
    panel.grid.major.y = element_line(color = "grey90"),
    panel.grid.minor   = element_blank(),
    plot.title         = element_text(face = "bold", size = 11, color = "#A50F15"),
    plot.subtitle      = element_text(size = 9)) +
  labs(
    title    = "B. TC Immune-Enriched unique (Stromal)",
    subtitle = "GO Cellular Component",
    x = "Gene Ratio", y = "")

p_stromal_reactome <- ggplot(stromal_reactome,
                             aes(x = GeneRatio, y = reorder(Term_clean, GeneRatio))) +
  geom_point(aes(size = Count, color = neg_log10_p)) +
  scale_color_gradientn(
    colors = c("#FEE5D9", "#FC9272", "#DE2D26", "#A50F15"),
    name   = "-log10\n(adj.p)") +
  scale_size_continuous(name = "Gene\nCount", range = c(3, 10)) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.y        = element_text(size = 9),
    panel.grid.major.y = element_line(color = "grey90"),
    panel.grid.minor   = element_blank(),
    plot.title         = element_text(face = "bold", size = 11, color = "#A50F15"),
    plot.subtitle      = element_text(size = 9)) +
  labs(
    title    = "D. TC Immune-Enriched unique (Stromal)",
    subtitle = "Reactome Pathways",
    x = "Gene Ratio", y = "")

# plots

p_inflammatory_cc <- ggplot(inflammatory_cc,
                            aes(x = GeneRatio, y = reorder(Term_clean, GeneRatio))) +
  geom_point(aes(size = Count, color = neg_log10_p)) +
  scale_color_gradientn(
    colors = c("#DEEBF7", "#9ECAE1", "#3182BD", "#08519C"),
    name   = "-log10\n(adj.p)") +
  scale_size_continuous(name = "Gene\nCount", range = c(3, 10)) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.y        = element_text(size = 9),
    panel.grid.major.y = element_line(color = "grey90"),
    panel.grid.minor   = element_blank(),
    plot.title         = element_text(face = "bold", size = 11, color = "#08519C"),
    plot.subtitle      = element_text(size = 9)) +
  labs(
    title    = "C. Shared Inflammatory",
    subtitle = "GO Cellular Component",
    x = "Gene Ratio", y = "")

p_inflammatory_reactome <- ggplot(inflammatory_reactome,
                                  aes(x = GeneRatio, y = reorder(Term_clean, GeneRatio))) +
  geom_point(aes(size = Count, color = neg_log10_p)) +
  scale_color_gradientn(
    colors = c("#DEEBF7", "#9ECAE1", "#3182BD", "#08519C"),
    name   = "-log10\n(adj.p)") +
  scale_size_continuous(name = "Gene\nCount", range = c(3, 10)) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.y        = element_text(size = 9),
    panel.grid.major.y = element_line(color = "grey90"),
    panel.grid.minor   = element_blank(),
    plot.title         = element_text(face = "bold", size = 11, color = "#08519C"),
    plot.subtitle      = element_text(size = 9)) +
  labs(
    title    = "E. Shared Inflammatory",
    subtitle = "Reactome Pathways",
    x = "Gene Ratio", y = "")

# figure (CC + REACTOME)

compact_comparison <- (p_stromal_cc | p_inflammatory_cc) /
  (p_stromal_reactome | p_inflammatory_reactome) +
  plot_annotation(
    title = "Molecular basis of immune exclusion in TC Immune-Enriched",
    subtitle = paste0(
      "Left: Cancer-intrinsic stromal program (", length(tc_unique), " genes) | ",
      "Right: Shared inflammatory program (", length(shared_inflammatory), " genes)"),
    caption  = "GO Cellular Component and Reactome pathways | Top 10 terms",
    theme    = theme(
      plot.title    = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5)
    )
  )

ggsave("figures/Figure_Compact_Stromal_vs_Inflammatory.png",
       compact_comparison, width = 14, height = 12, dpi = 300)

ggsave("figures/Figure_Compact_Stromal_vs_Inflammatory.pdf",
       compact_comparison, width = 14, height = 12)

# export tables

write.csv(stromal_cc,           "results/Table_TCunique_CellularComponent.csv",  row.names = FALSE)
write.csv(inflammatory_cc,      "results/Table_Shared_CellularComponent.csv",    row.names = FALSE)
write.csv(stromal_reactome,     "results/Table_TCunique_Reactome.csv",           row.names = FALSE)
write.csv(inflammatory_reactome,"results/Table_Shared_Reactome.csv",             row.names = FALSE)

