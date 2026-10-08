# Generic GO/Reactome enrichment + hub genes for ALL significant WGCNA modules 

library(enrichR)
library(dplyr)


# GO enrichment for all modules

dbs_go <- c("GO_Biological_Process_2025",
            "GO_Cellular_Component_2025",
            "GO_Molecular_Function_2025",
            "Reactome_2022")

all_enrichment_results <- list()

for (module in modules_significant) {
  cat("\nModule:", module, "(", length(module_gene_lists_sig[[module]]), "genes )\n")

  all_enrichment_results[[module]] <- enrichr(
    module_gene_lists_sig[[module]],
    dbs_go)

  top3 <- all_enrichment_results[[module]]$GO_Biological_Process_2025 |>
    arrange(Adjusted.P.value) |>
    head(3) |>
    select(Term, Overlap, Adjusted.P.value)

  cat("Top GO BP:\n")
  print(top3)
}

for (module in modules_significant) {
  write.csv(all_enrichment_results[[module]]$GO_Biological_Process_2025,
            paste0("results/GO_BP_", module, "_module.csv"), row.names = FALSE)
  write.csv(all_enrichment_results[[module]]$GO_Cellular_Component_2025,
            paste0("results/GO_CC_", module, "_module.csv"), row.names = FALSE)
  write.csv(all_enrichment_results[[module]]$Reactome_2022,
            paste0("results/Reactome_", module, "_module.csv"), row.names = FALSE)
}


# Hub genes 
hub_genes_all <- list()

for (module in modules_significant) {
  me_name <- paste0("ME", module)
  MM <- as.numeric(cor(norm.counts, MEs[, me_name], use = "p"))
  names(MM) <- colnames(norm.counts)

  genes_in_module <- names(moduleColors)[moduleColors == module]

  trait_rows <- all_sig[all_sig$Module == module, ]
  trait_rows <- trait_rows[order(trait_rows$Pvalue), ]
  best_trait_name <- as.character(trait_rows$Trait[1])
  best_trait_cor  <- trait_rows$Correlation[1]

  trait_vec <- datTraits[, best_trait_name]
  GS <- as.numeric(cor(norm.counts, trait_vec, use = "p"))
  names(GS) <- colnames(norm.counts)

  df <- data.frame(
    Gene       = genes_in_module,
    Module     = module,
    MM         = round(MM[genes_in_module], 4),
    GS         = round(GS[genes_in_module], 4),
    Best_Trait = best_trait_name,
    Trait_Cor  = best_trait_cor,
    stringsAsFactors = FALSE
  )

  df <- df[order(abs(df$MM), decreasing = TRUE), ]
  hub_genes_all[[module]] <- head(df, 20)
}

for (module in modules_significant) {
  cat("\n===", module, "— top 5 hub genes ===\n")
  print(hub_genes_all[[module]][1:5, c("Gene", "MM", "GS", "Best_Trait")])
}

all_hubs <- bind_rows(hub_genes_all)
write.csv(all_hubs, "results/WGCNA_hub_genes_all_modules.csv", row.names = FALSE)