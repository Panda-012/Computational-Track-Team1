# Load packages
library("affy")
library("limma")
library("hgu133plus2.db")
library("AnnotationDbi")
library("stringr")
library("affyPLM")

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

# Save the expression matrix
colnames(aggregatedNewSet) <- sub("\\.CEL", "", colnames(aggregatedNewSet))
write.table(aggregatedNewSet, file = "./Microarray_Data.tsv", sep = "\t", row.names = TRUE)