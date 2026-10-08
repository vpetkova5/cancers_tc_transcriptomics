library(omnideconv)
library(immunedeconv)
library(dplyr)
library(tidyr)
library(openxlsx)

# CIBERSORT must be downlowaded form Stanford — CIBERSORT.R and LM22.txt matrix
source("external/CIBERSORT.R")
set_cibersort_mat("external/LM22.txt")

# read the raw count data - must be in TPM normalisation
expr <- read.delim("data/count_table_19_Samples_tpm_data.txt")
# rownames = gene symbols, colnames = samples

# deconvoltion — immunedeconv  
res_quanti <- deconvolute(expr, "quantiseq")
res_epic   <- deconvolute(expr, "epic")
res_mcp    <- deconvolute(expr, "mcp_counter")
res_est    <- deconvolute(expr, "estimate")
res_timer  <- deconvolute_timer(expr, indications = rep("thca", ncol(expr)))

# CIBERSORT absolute — immunedeconv 
ciber_raw <- CIBERSORT(
  sig_matrix   = "external/LM22_temp.txt",
  mixture_file = "data/expr_for_cibersort.txt",
  perm         = 0,
  QN           = FALSE,
  saveLocation = "results/cibersort_results.txt"
)

# write the count tables
write.csv(as.data.frame(res_quanti), "results/quantiseq_results.csv", row.names = FALSE)
write.csv(as.data.frame(res_epic),   "results/epic_results.csv",      row.names = FALSE)
write.csv(as.data.frame(res_mcp),    "results/mcp_counter_results.csv", row.names = FALSE)
write.csv(as.data.frame(res_est),    "results/estimate_results.csv",  row.names = FALSE)
write.csv(as.data.frame(res_timer),  "results/timer_results.csv",     row.names = FALSE)

wb <- createWorkbook()
addWorksheet(wb, "quanTIseq");   writeData(wb, "quanTIseq",   as.data.frame(res_quanti))
addWorksheet(wb, "EPIC");        writeData(wb, "EPIC",        as.data.frame(res_epic))
addWorksheet(wb, "MCP-counter"); writeData(wb, "MCP-counter", as.data.frame(res_mcp))
addWorksheet(wb, "ESTIMATE");    writeData(wb, "ESTIMATE",    as.data.frame(res_est))
addWorksheet(wb, "TIMER");       writeData(wb, "TIMER",       as.data.frame(res_timer))
addWorksheet(wb, "CIBERSORT");   writeData(wb, "CIBERSORT",   as.data.frame(ciber_raw))
saveWorkbook(wb, "results/deconvolution_results.xlsx", overwrite = TRUE)
