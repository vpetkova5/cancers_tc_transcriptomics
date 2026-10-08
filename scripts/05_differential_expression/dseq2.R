library(DESeq2)
library(ggplot2)
library(dplyr)
library(tidyr)
library(EnhancedVolcano)

colData_TME4 <- read.csv("data/colData_TME4.csv", row.names = 1)
count_table_19_Samples_raw_count_data <- read.delim(
  "data/count_table_19_Samples_raw_count_data.txt", row.names = 1)

colData   <- colData_TME4
countData <- count_table_19_Samples_raw_count_data

# cheking colData and countData sample order
stopifnot(all(rownames(colData_TME4) == colnames(countData)))

print(unique(colData$TME_Class2))

# reference level (
colData_TME4$TME_Class2 <- factor(colData_TME4$TME_Class2)
colData_TME4$TME_Class2 <- relevel(colData_TME4$TME_Class2,
                                   ref = "TC Immune-Depleted")

# DESeq2 object
dds <- DESeqDataSetFromMatrix(
  countData = countData,
  colData   = colData_TME4,
  design    = ~ TME_Class2
)

# pre-filter: keep genes with >= 10 counts in at least 2 samples
keep <- rowSums(counts(dds) >= 10) >= 2
dds  <- dds[keep, ]

# DESeq2
dds <- DESeq(dds)

# all 6 contrasts - full results
res_TC_Enr_vs_Dep <- results(dds,
                             contrast      = c("TME_Class2", "TC Immune-Enriched", "TC Immune-Depleted"),
                             alpha         = 0.05,
                             pAdjustMethod = "BH")

# LFC > 0 (up) = 1476 genes → upregulated в TC Immune-Enriched vs TC Immune-Depleted
# LFC < 0 (down) = 710 genes → downregulated в TC Immune-Enriched

res_TC_Enr_vs_NG_High <- results(dds,
                                 contrast      = c("TME_Class2", "TC Immune-Enriched", "NG Immune-High"),
                                 alpha         = 0.05,
                                 pAdjustMethod = "BH")

res_TC_Enr_vs_NG_Low <- results(dds,
                                contrast      = c("TME_Class2", "TC Immune-Enriched", "NG Immune-Low"),
                                alpha         = 0.05,
                                pAdjustMethod = "BH")

res_TC_Dep_vs_NG_High <- results(dds,
                                 contrast      = c("TME_Class2", "TC Immune-Depleted", "NG Immune-High"),
                                 alpha         = 0.05,
                                 pAdjustMethod = "BH")

res_TC_Dep_vs_NG_Low <- results(dds,
                                contrast      = c("TME_Class2", "TC Immune-Depleted", "NG Immune-Low"),
                                alpha         = 0.05,
                                pAdjustMethod = "BH")

res_NG_High_vs_Low <- results(dds,
                              contrast      = c("TME_Class2", "NG Immune-High", "NG Immune-Low"),
                              alpha         = 0.05,
                              pAdjustMethod = "BH")

# filter padj < 0.05 & abs(log2FoldChange) >= 1
res_TC_Enr_vs_Dep_filt <- subset(as.data.frame(res_TC_Enr_vs_Dep),
                                 !is.na(padj) & padj < 0.05 & abs(log2FoldChange) >= 1)

res_TC_Enr_vs_NG_High_filt <- subset(as.data.frame(res_TC_Enr_vs_NG_High),
                                     !is.na(padj) & padj < 0.05 & abs(log2FoldChange) >= 1)

res_TC_Enr_vs_NG_Low_filt <- subset(as.data.frame(res_TC_Enr_vs_NG_Low),
                                    !is.na(padj) & padj < 0.05 & abs(log2FoldChange) >= 1)

res_TC_Dep_vs_NG_High_filt <- subset(as.data.frame(res_TC_Dep_vs_NG_High),
                                     !is.na(padj) & padj < 0.05 & abs(log2FoldChange) >= 1)

res_TC_Dep_vs_NG_Low_filt <- subset(as.data.frame(res_TC_Dep_vs_NG_Low),
                                    !is.na(padj) & padj < 0.05 & abs(log2FoldChange) >= 1)

res_NG_High_vs_Low_filt <- subset(as.data.frame(res_NG_High_vs_Low),
                                  !is.na(padj) & padj < 0.05 & abs(log2FoldChange) >= 1)

# export 

write.csv(res_TC_Enr_vs_Dep_filt,
          "results/res_filt_TC_Enriched_vs_Depleted.csv")
write.csv(res_TC_Enr_vs_NG_High_filt,
          "results/res_filt_TC_Enriched_vs_NG_High.csv")
write.csv(res_TC_Enr_vs_NG_Low_filt,
          "results/res_filt_TC_Enriched_vs_NG_Low.csv")
write.csv(res_TC_Dep_vs_NG_High_filt,
          "results/res_filt_TC_Depleted_vs_NG_High.csv")
write.csv(res_TC_Dep_vs_NG_Low_filt,
          "results/res_filt_TC_Depleted_vs_NG_Low.csv")
write.csv(res_NG_High_vs_Low_filt,
          "results/res_filt_NG_High_vs_Low.csv")
