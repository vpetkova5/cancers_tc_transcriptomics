# Analysis code: transcriptomic convergence and stromal remodelling in immune-enriched papillary thyroid carcinoma

R scripts accompanying the manuscript:

> Petkova V, Dodova R, Vidinov K, Kaneva R. Transcriptomic convergence reveals
> tumour-specific stromal remodelling in immune-enriched papillary thyroid
> carcinoma. Cancers 2026, 18. https://doi.org/10.3390/xxxxx

The study is a pilot, hypothesis-generating bulk RNA-seq analysis of 19 surgical
specimens (10 papillary thyroid carcinoma, 9 nodular goitre), followed by
evaluation of the derived 441-gene stromal programme in TCGA-THCA.

## Repository layout

| Folder | Content |
|---|---|
| `scripts/01_deconvolution` | Immune/stromal deconvolution (ESTIMATE, quanTIseq, EPIC, MCP-counter, TIMER, CIBERSORT) |
| `scripts/02_tme_classification` | Median-split TME subgroups (TC-IE, TC-ID, NG-IH, NG-IL) and quanTIseq statistics |
| `scripts/03_descriptive_stats` | Descriptive statistics, Kruskal-Wallis / Dunn tests, EPIC fractions |
| `scripts/04_caf_analysis` | CAF fractions, correlation with ESTIMATE scores, MCP-counter cross-check |
| `scripts/05_differential_expression` | DESeq2 (six contrasts), PCA, volcano plots, up/down split, intersection analysis |
| `scripts/06_ecm_stromal` | Collagen-containing ECM genes, GO/Reactome enrichment of the stromal and shared gene sets |
| `scripts/07_wgcna` | WGCNA network construction, module-trait correlations, generic module enrichment |
| `scripts/08_tcga_validation` | Stromal score vs. clinical risk categories in TCGA-THCA |
| `scripts/09_tcga_validation_figures` | ssGSEA stromal score, TDS, TCGA-THCA panels and main figures |
| `scripts/10_roc` | Reserved for ROC  |

## How to run

The scripts share R objects (e.g. `dds`, `tme_df`, `epic_long`, `tc_unique`), so
they are written to be run in folder order, within one R session, from the
repository root. Files inside a folder are run in the order of their numeric
prefix where one exists. 

Scripts read from `data/` and write to `results/` and `figures/`; these folders
must exist (they are created empty in the repository).

## Requirements

- R 4.2.2
- DESeq2 1.38.3, WGCNA 1.74
- Other packages used: omnideconv, immunedeconv, dplyr, tidyr, tidyverse,
  openxlsx, ggplot2, ggpubr, rstatix, patchwork, ggrepel, reshape2, Hmisc,
  pheatmap, EnhancedVolcano, enrichR, epitools

## Data

- Raw and processed RNA-seq data of the discovery cohort are not publicly
  available because they are part of an ongoing study. Requests should be
  directed to Radka Kaneva (see the Data Availability Statement of the paper).
- The TCGA-THCA data (456 papillary thyroid carcinomas) were retrieved with
  TCGAbiolinks (see Methods 2.12).
- Read processing and quantification (FastQC, Trimmomatic, STAR, featureCounts)
  are described in the manuscript, Methods 2.3-2.4.

## External dependency: CIBERSORT

CIBERSORT is available only under an academic licence from Stanford and is
not included in this repository. Place `CIBERSORT.R` and `LM22.txt` in
`external/` after obtaining a licence.

## Contact

Veronika Petkova, Molecular Medicine Center, Medical University of Sofia -
vpetkova@mmcbg.org

## Licence

To be added.
