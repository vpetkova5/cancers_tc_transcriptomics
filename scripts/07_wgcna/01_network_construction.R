library(WGCNA)
library(DESeq2)
library(tidyverse)

# multi-threading
allowWGCNAThreads()

#  `countData`/`colData_TME4` from scripts/06_differential_expression/dseq2.R

# normalisation for WGCNA — VST without design
dds_wgcna <- DESeqDataSetFromMatrix(
  countData = countData,
  colData   = colData_TME4,
  design    = ~ 1)  

# filter low expressed genes as DESeq2
keep_wgcna <- rowSums(counts(dds_wgcna) >= 10) >= 2
dds_wgcna  <- dds_wgcna[keep_wgcna, ]
cat("Genes after the filter:", nrow(dds_wgcna), "\n")

# VST norm
dds_norm <- vst(dds_wgcna, blind = TRUE)

# t (samples x genes)
norm.counts <- assay(dds_norm) |> t()

cat("Matrix dim:", nrow(norm.counts), "samples x", ncol(norm.counts), "genes\n")
cat("Sample names:\n")
print(rownames(norm.counts))

# QC и outlier detection

gsg <- goodSamplesGenes(norm.counts)
summary(gsg)
gsg$allOK

table(gsg$goodGenes)
table(gsg$goodSamples)

if (!gsg$allOK) {
  norm.counts <- norm.counts[gsg$goodSamples, gsg$goodGenes]
  cat("After QC filter:", nrow(norm.counts), "samples x", ncol(norm.counts), "genes\n")
} else {
  cat("All genes and samples QC\n")
}

# hierarchical clustering for outlier detection
htree <- hclust(dist(norm.counts), method = "average")

png("figures/WGCNA_sample_clustering.png", width = 14, height = 8, units = "in", res = 300)
par(mar = c(0, 5, 4, 0))
plot(htree, main = "Sample clustering — outlier detection", sub = "", xlab = "", cex = 0.8)
dev.off()

# samples.to.exclude - NG74 (outlier after clustering)
samples.to.exclude <- c("NG74")
norm.counts <- norm.counts[!rownames(norm.counts) %in% samples.to.exclude, ]
cat("Samples after filter:", nrow(norm.counts), "\n")
print(rownames(norm.counts))

colData_wgcna <- colData_TME4[!rownames(colData_TME4) %in% samples.to.exclude, ]
cat("Samples in colData:", nrow(colData_wgcna), "\n")

# clustering 
htree2 <- hclust(dist(norm.counts), method = "average")

png("figures/WGCNA_sample_clustering_noNG74.png", width = 14, height = 8, units = "in", res = 300)
par(mar = c(0, 5, 4, 0))
dev.off()

# Soft threshold
power <- c(c(1:10), seq(from = 12, to = 30, by = 2))

sft <- pickSoftThreshold(norm.counts,
                         powerVector  = power,
                         networkType  = "signed",
                         verbose      = 5)

sft.data <- sft$fitIndices

png("figures/WGCNA_soft_threshold.png", width = 10, height = 8, units = "in", res = 300)
par(mfrow = c(1, 2))
plot(sft.data$Power, sft.data$SFT.R.sq,
     xlab = "Soft threshold (power)", ylab = "Scale free topology R²",
     main = "Scale independence", type = "n")
text(sft.data$Power, sft.data$SFT.R.sq, labels = sft.data$Power, cex = 0.9, col = "red")
abline(h = 0.8, col = "red", lty = 2)
plot(sft.data$Power, sft.data$mean.k.,
     xlab = "Soft threshold (power)", ylab = "Mean connectivity",
     main = "Mean connectivity", type = "n")
text(sft.data$Power, sft.data$mean.k., labels = sft.data$Power, cex = 0.9, col = "red")
dev.off()

# NOTE: soft_power below was chosen by inspecting the plot above (first
# power reaching R² ≈ 0.8 with reasonable mean connectivity) — a manual,
# visual decision. 
soft_power <- 18

# Етап 4 — Network construction

# WGCNA::cor overrides base cor while blockwiseModules runs, to avoid a namespace clash
temp_cor <- cor
cor <- WGCNA::cor

net <- blockwiseModules(norm.counts,
                        power             = soft_power,
                        networkType       = "signed",
                        TOMType           = "signed",
                        maxBlockSize      = 18000,
                        minModuleSize     = 40,
                        reassignThreshold = 0,
                        mergeCutHeight    = 0.25,
                        numericLabels     = TRUE,
                        pamRespectsDendro = FALSE,
                        saveTOMs          = FALSE,
                        verbose           = 3)

cor <- temp_cor

# Module summary
moduleLabels <- net$colors
moduleColors <- labels2colors(moduleLabels)


n_modules <- length(unique(moduleLabels)) - 1  # -1 excludes the grey/unassigned label (0)
n_samples <- nrow(norm.counts)

cat("\n Module Summary \n")
cat("Genes:", length(moduleLabels), "\n")
cat("Modules:", n_modules, "\n")
cat("Unassigned (grey):", sum(moduleLabels == 0), "\n\n")

cat("Module sizes:\n")
module_sizes <- sort(table(moduleColors), decreasing = TRUE)
print(module_sizes)

# Dendrogram
png("figures/WGCNA_dendrogram.png", width = 16, height = 9, units = "in", res = 300)
plotDendroAndColors(net$dendrograms[[1]],
                    moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE,
                    hang         = 0.03,
                    addGuide     = TRUE,
                    guideHang    = 0.05,
                    main         = paste0("Gene co-expression network dendrogram\n",
                                          "(Power=", soft_power, ", n=", n_samples,
                                          " samples, ", n_modules, " modules)"))
dev.off()


# Етап 5 — Module-Trait correlations

MEs0 <- moduleEigengenes(norm.counts, moduleColors)$eigengenes
MEs  <- orderMEs(MEs0)
cat("Module eigengenes изчислени за", ncol(MEs), "модула\n")

datTraits <- data.frame(
  TC_Immune_Enriched = as.numeric(colData_wgcna$TME_Class2 == "TC Immune-Enriched"),
  TC_Immune_Depleted = as.numeric(colData_wgcna$TME_Class2 == "TC Immune-Depleted"),
  NG_Immune_High     = as.numeric(colData_wgcna$TME_Class2 == "NG Immune-High"),
  NG_Immune_Low      = as.numeric(colData_wgcna$TME_Class2 == "NG Immune-Low"),
  Malignant          = as.numeric(colData_wgcna$TME_Class2 %in%
                                    c("TC Immune-Enriched", "TC Immune-Depleted")),
  row.names          = rownames(colData_wgcna)
)

cat("Trait summary:\n")
print(colSums(datTraits))

moduleTraitCor    <- cor(MEs, datTraits, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(norm.counts))

cat("\nStat sig modules (p < 0.05) \n")
for (trait in colnames(datTraits)) {
  sig <- sum(moduleTraitPvalue[, trait] < 0.05)
  cat(trait, ":", sig, "Stat. sig \n")
}

# --- Heatmap 1: correlation + significance stars ---
textMatrix_stars <- matrix("", nrow = nrow(moduleTraitCor), ncol = ncol(moduleTraitCor))
for (i in 1:nrow(moduleTraitPvalue)) {
  for (j in 1:ncol(moduleTraitPvalue)) {
    p    <- moduleTraitPvalue[i, j]
    corv <- signif(moduleTraitCor[i, j], 2)
    stars <- ifelse(p < 0.001, "***",
                    ifelse(p < 0.01,  "**",
                           ifelse(p < 0.05,  "*", "")))
    textMatrix_stars[i, j] <- paste0(corv, "\n", stars)
  }
}
dim(textMatrix_stars) <- dim(moduleTraitCor)

heatmap_title <- paste0("Module-Trait Relationships\n(n=", n_samples,
                        " samples, ", n_modules, " modules)")

png("figures/WGCNA_module_trait_heatmap_stars.png",
    width = 10, height = 16, units = "in", res = 300)
par(mar = c(6, 10, 3, 3))
labeledHeatmap(Matrix        = moduleTraitCor,
               xLabels       = colnames(datTraits),
               yLabels       = names(MEs),
               ySymbols      = names(MEs),
               colorLabels   = FALSE,
               colors        = blueWhiteRed(50),
               textMatrix    = textMatrix_stars,
               setStdMargins = FALSE,
               cex.text      = 0.55,
               zlim          = c(-1, 1),
               main          = heatmap_title)
dev.off()

# Heatmap 2: correlation + p-value 
textMatrix_pval <- paste(signif(moduleTraitCor, 2), "\n(",
                         signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix_pval) <- dim(moduleTraitCor)

png("figures/WGCNA_module_trait_heatmap_pvalues.png",
    width = 10, height = 16, units = "in", res = 300)
par(mar = c(6, 10, 3, 3))
labeledHeatmap(Matrix        = moduleTraitCor,
               xLabels       = colnames(datTraits),
               yLabels       = names(MEs),
               ySymbols      = names(MEs),
               colorLabels   = FALSE,
               colors        = blueWhiteRed(50),
               textMatrix    = textMatrix_pval,
               setStdMargins = FALSE,
               cex.text      = 0.55,
               zlim          = c(-1, 1),
               main          = heatmap_title)
dev.off()

# Stat sig moduls per trait + gene lists 

sig_modules_list <- list()
for (trait in colnames(datTraits)) {
  sig_idx <- which(moduleTraitPvalue[, trait] < 0.05)
  if (length(sig_idx) > 0) {
    sig_modules_list[[trait]] <- data.frame(
      Module      = gsub("ME", "", rownames(moduleTraitCor)[sig_idx]),
      Correlation = round(moduleTraitCor[sig_idx, trait], 3),
      Pvalue      = round(moduleTraitPvalue[sig_idx, trait], 4),
      Direction   = ifelse(moduleTraitCor[sig_idx, trait] > 0, "positive", "negative"),
      Trait       = trait
    ) |> arrange(Pvalue)
  }
}

for (trait in names(sig_modules_list)) {
  cat("\n===", trait, "===\n")
  print(sig_modules_list[[trait]])
}

all_sig <- bind_rows(sig_modules_list)
write.csv(all_sig, "results/WGCNA_significant_modules.csv", row.names = FALSE)

modules_significant <- unique(all_sig$Module)
print(modules_significant)

# genes for each module
module_gene_lists_sig <- list()
for (module in modules_significant) {
  genes <- names(moduleColors)[moduleColors == module]
  module_gene_lists_sig[[module]] <- genes
  cor_info <- all_sig[all_sig$Module == module, c("Trait", "Correlation", "Pvalue", "Direction")]
  cat("\n---", module, "---\n")
  cat("Genes:", length(genes), "\n")
  print(cor_info)
}

all_sig_genes <- bind_rows(lapply(modules_significant, function(m) {
  genes <- module_gene_lists_sig[[m]]
  cor_info <- all_sig[all_sig$Module == m, ]
  best_trait <- cor_info[which.min(cor_info$Pvalue), ]
  data.frame(
    Gene        = genes,
    Module      = m,
    Best_Trait  = best_trait$Trait,
    Best_Cor    = best_trait$Correlation,
    Best_Pvalue = best_trait$Pvalue,
    Direction   = best_trait$Direction,
    stringsAsFactors = FALSE
  )
}))

write.csv(all_sig_genes, "results/WGCNA_all_significant_modules_genes.csv", row.names = FALSE)

for (module in modules_significant) {
  write.table(module_gene_lists_sig[[module]],
              file  = paste0("results/genes_sig_", module, "_module.txt"),
              quote = FALSE, row.names = FALSE, col.names = FALSE)
}
