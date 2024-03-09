# Load necessary libraries
library(readr)
library(limma)
library(Rtsne)
library(edgeR)
library(glmmTMB)
library(data.table)
library(reshape2) # For melt function
library(dplyr) # For data manipulation
library(jsonlite)  # For loading metadata and clinical data

# Set working directory and load the data
setwd(normalizePath(dirname(rstudioapi::getSourceEditorContext()$path)))
getwd()

# Load RNA-Seq count matrix
count_matrix <- read.delim("RNA-Seq_Data.tsv")
colnames(count_matrix) <- gsub("X", "#", gsub("\\.", "-", gsub("^X", "", colnames(count_matrix))))
rownames(count_matrix) <- as.character(count_matrix[,1])
count_matrix <- count_matrix[,-1]
count_matrix <- as.matrix(count_matrix)  # Convert to numeric matrix

# Load metadata and clinical data
metadata <- fromJSON("RNA-Seq_Metadata.json")
clinical_data <- fromJSON("RNA-Seq_Clinical_Data.json")

# Create an empty column for case_id in the metadata
metadata <- cbind(metadata[, 1:3], case_id = NA, metadata[, 4:ncol(metadata)])

# Loop through each entry in the metadata to find its case id
for (i in 1:nrow(metadata)) {
  matching_record <- metadata[grep(metadata$file_name[[i]], metadata$file_name), ]
  metadata$case_id[i] <- matching_record$associated_entities[[1]]$case_id
}

# --- Quality control and filtering ---
rm (matching_record,i)

# Check for duplicate gene names
sum(duplicated(rownames(count_matrix)))

# Identify missing values
sum(is.na(count_matrix))
sum(count_matrix == 0)

# Create dge, estimate dispersions, and normalize the data using TMM
dge <- DGEList(counts = count_matrix)
dge <- estimateDisp(dge)
dge <- calcNormFactors(dge, method = "TMM")  # Perform TMM normalization within DGEList

# Subtyping of samples via external R script
Environment_Data <- ls()
Sisi <- cpm(dge, normalized.lib.size = TRUE)
source("Breast_Cancer_Subtyping.R")
table(PAM50Preds$subtype)
rm(list=setdiff(ls(), c("LumA", "LumB", "Basal", "Her2", "Normal", Environment_Data)))

# Create dge, estimate dispersions, and normalize the data using TMM
dge <- DGEList(counts = count_matrix[, c(LumA, LumB), drop = FALSE])
dge <- estimateDisp(dge)
dge <- calcNormFactors(dge, method = "TMM")  # Perform TMM normalization within DGEList

# Attach sample information dataframe as an additional component to the dge list
dge$samples <- cbind(dge$samples, data.frame(sample = colnames(dge$counts),
                          subtype = rep(NA, length(colnames(dge$counts)))))

dge$samples$subtype <- ifelse(colnames(dge$counts) %in% LumA, "Luminal_A",
                       ifelse(colnames(dge$counts) %in% LumB, "Luminal_B", NA))

##################################################################################################
# # Fit the data to the Zero-Inflated Negative Binomial (ZINB) Model requirements
# count_matrix <- as.data.table(round(dge$counts))
# count_matrix[, gene := rownames(dge$counts)]
# 
# total_counts_raw <- colSums(dge$counts)
# names(total_counts_raw) <- colnames(dge$counts)
# 
# count_matrix_long <- setDT(melt(count_matrix, id.vars = "gene", variable.name = "sample", value.name = "gene_expression"))
# count_matrix_long[, total_counts := total_counts_raw[sample], by = sample]
# 
# final_df <- merge(count_matrix_long, dge$samples [, c("sample" , "subtype")], by = "sample", all.x = TRUE)
# count_matrix <- as.data.frame (count_matrix)
# rownames (count_matrix) <- count_matrix$gene
# count_matrix <- as.data.frame (count_matrix [, -122])
# 
# 
# # Fit the Zero-Inflated Negative Binomial (ZINB) Model to the data
# zinb_model <- glmmTMB(gene_expression ~ subtype + offset(log(total_counts)),
#                       zi = ~ subtype,
#                       family = nbinom2,  # This specifies a negative binomial distribution
#                       data = final_df)
# 
# final_df$predicted_counts <- predict(zinb_model, type = "response")
# final_df$zero_inflation_prob <- predict(zinb_model, type = "zprob")
##################################################################################################

# Identify outliers (visually and statistically)
norm_counts <- cpm(dge, normalized.lib.size = TRUE)
length(boxplot.stats(norm_counts)$out)
boxplot(norm_counts)

# Quality check
affy::plotDensity(norm_counts)
plot(density(apply(norm_counts, 2, mean, na.rm = TRUE)),main="omics",cex.axis=0.5) #1 rows, 2 columns
plotMDS(norm_counts)


# Filter lowly expressed genes
keep <- rowSums(norm_counts >= 10) >= 5  # Adjust thresholds as needed
norm_counts <- norm_counts[keep, ]


# Reassess quality after normalization
boxplot(count_matrix_norm)
plotMDS(count_matrix_norm)

# --- Dimensionality reduction ---

# Perform PCA
pca_results <- prcomp(t(data_matrix_normalized), center = TRUE, scale. = TRUE)

# Visualize PCA results
plot(pca_results$x, col = sample_colors)
legend("topright", legend = c( "Luminal_B", "Luminal_A", "Normal"), col = unique(sample_colors), pch = 1, cex = 0.8)


# Create Elbow Graph
variance_explained <- pca_results$sdev^2 / sum(pca_results$sdev^2) * 100
plot(1:length(variance_explained), variance_explained, type = "b", pch = 19)
abline(h = 3, col = "red", lty = 2)  # Optional threshold for selecting PCs

# Perform t-SNE
set.seed(123)  # Set seed for reproducibility
tsne_results <- Rtsne(pca_results$x[, 1:5], perplexity = 23)  # Adjust perplexity if needed

# Visualize t-SNE results
plot(tsne_results$Y, col = sample_colors)
legend("bottomleft", legend = c( "Luminal_B", "Luminal_A", "Normal"), col = unique(sample_colors), pch = 1, cex = 0.8)

