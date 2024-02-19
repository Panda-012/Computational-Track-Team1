install.packages(c("BiocManager", "dplyr", "glmnet", "data.table", "survival", "preprocessCore", "limma", "pracma", "umap"))
BiocManager::install(c("IlluminaHumanMethylation450kanno.ilmn12.hg19", "TCGAbiolinks"))

packageVersion("BiocManager")
packageVersion("dplyr")
packageVersion("glmnet")
packageVersion("data.table")
packageVersion("survival")
packageVersion("preprocessCore")
packageVersion("limma")
packageVersion("IlluminaHumanMethylation450kanno.ilmn12.hg19")
packageVersion("TCGAbiolinks")

subset_matrix <- Sisi[, 1:5]
write.table(subset_matrix, file = "GE.tsv", sep = "\t", quote = FALSE, row.names = TRUE)

# Assuming `run.intend.using.tcga.lvo.preiction.model.customized` is the function you plan to use
# Load your preprocessed data
exp_data <- read.delim("GE.tsv", row.names = 1)
meth_data <- read.delim("DM.tsv", row.names = 1)

# Combine into a list expected by the function
exp_met_data <- list(expression = exp_data, methylation = meth_data)

library(parallel)
setwd("C:/Users/Panda/Desktop/Graduation/Scripts/INTEND-main")
source('intend.R')

# Example execution - adjust parameters as needed
run.intend.using.tcga.lvo.preiction.model.customized(
  exp.met.data = exp_met_data, 
  lvo.subtype = "Breast",  # Specify your subtype of interest
  custom.subtype = "Luminal",
  verbose = TRUE
)

