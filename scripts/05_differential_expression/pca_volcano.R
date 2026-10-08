library(ggplot2)
library(ggrepel)
library(DESeq2)
library(dplyr)
library(patchwork)

# PCA 
vsd <- vst(dds, blind = FALSE)

pca_data <- plotPCA(vsd,
                    intgroup = "TME_Class2",
                    returnData = TRUE)

percent_var <- round(100 * attr(pca_data, "percentVar"))

# sample names
pca_data$sample <- rownames(pca_data)

colors4 <- c(
  "NG Immune-High"     = "forestgreen",
  "NG Immune-Low"      = "darkolivegreen3",
  "TC Immune-Enriched" = "firebrick",
  "TC Immune-Depleted" = "steelblue"
)

pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2,
                     color = TME_Class2,
                     label = sample)) +
  geom_point(size = 3) +
  geom_text_repel(size = 3,
                  max.overlaps = 20,
                  box.padding  = 0.4) +
  scale_color_manual(values = colors4, name = "group") +
  labs(
    x = paste0("PC1: ", percent_var[1], "% variance"),
    y = paste0("PC2: ", percent_var[2], "% variance")
  ) +
  theme_classic(base_size = 12)

ggsave("figures/pca_plot.png", pca_plot, width = 12, height = 8, dpi = 300)

# function for volcano plot
make_volcano <- function(res, title) {
  df <- as.data.frame(res) |>
    filter(!is.na(padj)) |>
    mutate(
      category = case_when(
        padj < 0.05 & log2FoldChange >= 1  ~ "up",
        padj < 0.05 & log2FoldChange <= -1 ~ "down",
        TRUE                                ~ "ns"
      )
    )

  counts <- table(df$category)
  up_n   <- ifelse("up"   %in% names(counts), counts["up"],   0)
  down_n <- ifelse("down" %in% names(counts), counts["down"], 0)

  ggplot(df, aes(x = log2FoldChange, y = -log10(padj), color = category)) +
    geom_point(alpha = 0.5, size = 1) +
    scale_color_manual(
      values = c("up" = "#E24B4A", "down" = "#378ADD", "ns" = "grey70"),
      labels = c("up"   = paste0("Up (", up_n, ")"),
                 "down" = paste0("Down (", down_n, ")"),
                 "ns"   = "Not sig"),
      name = NULL
    ) +
    geom_vline(xintercept = c(-1, 1),
               linetype = "dashed", color = "grey50", linewidth = 0.4) +
    geom_hline(yintercept = -log10(0.05),
               linetype = "dashed", color = "grey50", linewidth = 0.4) +
    labs(
      title    = title,
      subtitle = paste0("DEGs: ", up_n + down_n,
                        "  (|log2FC| ≥ 1, padj < 0.05)"),
      x        = "Log2 fold change",
      y        = "-log10(padj)"
    ) +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "top",
      plot.title      = element_text(size = 11, face = "bold"),
      plot.subtitle   = element_text(size = 10, color = "grey40")
    )
}

# all 6 volcano plots
p1 <- make_volcano(res_TC_Enr_vs_Dep,    "TC Immune-Enriched vs TC Immune-Depleted")
p2 <- make_volcano(res_TC_Enr_vs_NG_High,"TC Immune-Enriched vs NG Immune-High")
p3 <- make_volcano(res_TC_Enr_vs_NG_Low, "TC Immune-Enriched vs NG Immune-Low")
p4 <- make_volcano(res_TC_Dep_vs_NG_High,"TC Immune-Depleted vs NG Immune-High")
p5 <- make_volcano(res_TC_Dep_vs_NG_Low, "TC Immune-Depleted vs NG Immune-Low")
p6 <- make_volcano(res_NG_High_vs_Low,   "NG Immune-High vs NG Immune-Low")

# 2x3 grid
combined <- (p1 | p2 | p3) / (p4 | p5 | p6)

ggsave("figures/volcano_6contrasts.png",
       combined, width = 16, height = 10, dpi = 300)

ggsave("figures/p1.png", p1, width = 8, height = 6, dpi = 300)
ggsave("figures/p2.png", p2, width = 8, height = 6, dpi = 300)
ggsave("figures/p3.png", p3, width = 8, height = 6, dpi = 300)
ggsave("figures/p4.png", p4, width = 8, height = 6, dpi = 300)
ggsave("figures/p5.png", p5, width = 8, height = 6, dpi = 300)
ggsave("figures/p6.png", p6, width = 8, height = 6, dpi = 300)
