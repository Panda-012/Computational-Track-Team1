# Load packages
library(affy)
library(limma)
library(hgu133plus2.db)
library(AnnotationDbi)
library(stringr)
library(affyPLM)
library(sva)

# Set working directory
setwd(normalizePath(dirname(rstudioapi::getSourceEditorContext()$path)))
getwd()

# Read CEL files
data <- ReadAffy(celfile.path = "./GSE45827_RAW")

# Background correction and normalization
bg.corrected <- affy::rma(data)
exprSet <- exprs(bg.corrected)

# Quality assessment
affy::plotDensity(exprSet)
plot(density(apply(exprSet, 2, mean, na.rm = TRUE)),main="omics",cex.axis=0.5) #1 rows, 2 columns

# Batch effect detection
boxplot(exprSet)
plotMDS(exprSet)

# Convert probe IDs to gene names
Gene_Names <- unlist(mget(rownames(exprSet), hgu133plus2SYMBOL))
rownames(exprSet) <- Gene_Names

# Remove NA rows
NewSet <- exprSet[!str_detect(rownames(exprSet), "^NA\\."), ]
NewSet <- na.omit(NewSet)

# Handle duplicate gene names
anyDuplicated(rownames(NewSet))
aggregatedNewSet <- aggregate(NewSet, by = list(rownames(NewSet)), FUN = mean)  # Use mean for aggregation
rownames(aggregatedNewSet) <- as.character(aggregatedNewSet[,1])
aggregatedNewSet = aggregatedNewSet[,-1]
anyDuplicated(rownames(aggregatedNewSet))

# Reassess quality after summarization
affy::plotDensity(aggregatedNewSet)
plot(density(apply(aggregatedNewSet, 2, mean, na.rm = TRUE)),main="omics",cex.axis=0.5) #1 rows, 2 columns
boxplot(aggregatedNewSet)
plotMDS(aggregatedNewSet)

# Create a model matrix for the different conditions
conditions <- factor(c(rep("Luminal_A", 29), rep("Luminal_B", 30), rep("Normal", 11)))
mod <- model.matrix(~ conditions)
mod0 <- model.matrix(~ 1, data = aggregatedNewSet) # Null model without the conditions
sva.obj <- sva(aggregatedNewSet, mod, mod0) # Run SVA to find surrogate variables (SVs)
svs <- sva.obj$sv # Extract the surrogate variables (SVs)
aggregatedNewSet.adj <- removeBatchEffect(aggregatedNewSet, covariates = mod, batch = svs) # Adjust the expression data for the surrogate variables

# Reassess quality after batch effect removal
affy::plotDensity(aggregatedNewSet.adj)
plot(density(apply(aggregatedNewSet.adj, 2, mean, na.rm = TRUE)),main="Post-SVA",cex.axis=0.5) #1 rows, 2 columns
boxplot(aggregatedNewSet.adj)
plotMDS(aggregatedNewSet.adj)

# Save the expression matrix
colnames(aggregatedNewSet.adj) <- sub("\\.CEL", "", colnames(aggregatedNewSet.adj))
write.table(aggregatedNewSet.adj, file = "./Microarray_Data.tsv", sep = "\t", row.names = TRUE)