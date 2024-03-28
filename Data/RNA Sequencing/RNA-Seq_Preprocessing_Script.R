# Load necessary libraries
library(readr)
library(limma)
library(Rtsne)
library(edgeR)
library(glmmTMB)
library(e1071)
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
rm (matching_record,i)

# Check for duplicate gene names
sum(duplicated(rownames(count_matrix)))

# Identify missing values
sum(is.na(count_matrix))
sum(count_matrix == 0)

# Create dge and normalize the data using TMM
dge <- DGEList(counts = count_matrix)
dge <- calcNormFactors(dge, method = "TMM")  # Perform TMM normalization within DGEList

# Subtyping of samples via external R script
Environment_Data <- ls()
Sisi <- cpm(dge, normalized.lib.size = TRUE)
source("Breast_Cancer_Subtyping.R")
table(PAM50Preds$subtype)
save(LumA, LumB, Her2, Basal, Normal, file = "../luminal_samples.RData")


# Filter the selected samples and start renaming them
data_matrix <- count_matrix[, c(LumA, LumB), drop = FALSE]
combined_samples <- c(LumA, LumB)
new_sample_names <- character(length(combined_samples))

# Loop through the combined sample names to generate new names
for (i in seq_along(combined_samples)) {
  prefix <- ifelse(combined_samples[i] %in% LumA, "LumA_", "LumB_")
  seq_num <- sum(startsWith(new_sample_names, prefix)) + 1
  new_sample_names[i] <- paste0(prefix, sprintf("%02d", seq_num))
}

colnames(data_matrix) <- new_sample_names
sample_colors <- ifelse(startsWith(colnames(data_matrix), "LumB"), "red", "orange")

# Assess quality before processing
affy::plotDensity(data_matrix)
plot(density(apply(data_matrix, 2, mean, na.rm = TRUE)),main="omics",cex.axis=0.5) #1 rows, 2 columns
plotMDS(data_matrix)

# Calculate median expression level across genes and median absolute deviation (MAD) for each gene
median_expression <- apply(data_matrix, 1, median)
mad_expression <- apply(data_matrix, 1, function(x) mad(x, constant = 1))

# Dynamically set mad_multiplier based on skewness
data_skewness <- skewness(log1p(rowMedians(data_matrix)))
mad_multiplier <- ifelse(data_skewness <= 1, 1.5, ifelse(data_skewness <= 2, 2, 2.5))

# Combine median and MAD for a robust measure of gene expression variability
threshold_expression <- median_expression + (mad_multiplier * mad_expression)

# Dynamically set samples_threshold based on overall expression density using a gradient approach
overall_non_zero_proportion <- mean(rowSums(data_matrix > 0)) / ncol(data_matrix)
samples_threshold <- 0.05 + (0.15 * (1 - overall_non_zero_proportion))

# Identify genes to keep based on the threshold and filter data based on the keep_genes logical vector
keep_genes <- rowSums(data_matrix >= threshold_expression) >= (samples_threshold * ncol(data_matrix))
filtered_data_matrix <- data_matrix[keep_genes, ]

# After filtering, clean up the environment from temporary variables
rm(list=setdiff(ls(), c("LumA", "LumB", "Basal", "Her2", "Normal", Environment_Data, "data_matrix", "filtered_data_matrix", "sample_colors")))

# Reassess quality after filtration
affy::plotDensity(filtered_data_matrix)
plot(density(apply(filtered_data_matrix, 2, mean, na.rm = TRUE)),main="omics",cex.axis=0.5) #1 rows, 2 columns
plotMDS(filtered_data_matrix)

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
# boxplot(filtered_data_matrix)
length(boxplot.stats(filtered_data_matrix)$out)

# Create dge and normalize the data using TMM
dge <- DGEList(counts = filtered_data_matrix)
dge <- calcNormFactors(dge, method = "TMM")  # Perform TMM normalization within DGEList

# Apply log2 transformation to stabilize variance
log_data_matrix <- cpm(dge, log = TRUE, prior.count = 1)
length(boxplot.stats(log_data_matrix)$out)

# Reassess quality after log transformation
affy::plotDensity(log_data_matrix)
plot(density(apply(log_data_matrix, 2, mean, na.rm = TRUE)),main="omics",cex.axis=0.5) #1 rows, 2 columns
plotMDS(log_data_matrix)

# Function to adjust outliers directly within a single vector of gene expression values
Environment_Data <- ls()
adjust_outliers_directly <- function(gene_expression) {
  # Constants for outlier detection
  k <- 3 # Adjust k as needed for stricter or looser criteria
  
  # Compute lower and upper bounds for this gene
  gene_mad <- mad(gene_expression, constant = 1)
  gene_median <- median(gene_expression)
  lb <- gene_median - k * gene_mad
  ub <- gene_median + k * gene_mad
  
  # Adjust outliers to the nearest value within the non-outlier range
  gene_expression[gene_expression < lb] <- lb
  gene_expression[gene_expression > ub] <- ub
  
  return(gene_expression)
}

# Apply the direct adjustment of outliers across all genes
adjusted_log_data_matrix <- t(apply(log_data_matrix, 1, adjust_outliers_directly))
length(boxplot.stats(adjusted_log_data_matrix)$out)
rm(list=setdiff(ls(), c(Environment_Data, "adjust_outliers_directly", "adjusted_log_data_matrix")))

# Reassess quality after handling outliers
affy::plotDensity(adjusted_log_data_matrix)
plot(density(apply(adjusted_log_data_matrix, 2, mean, na.rm = TRUE)),main="omics",cex.axis=0.5) #1 rows, 2 columns
plotMDS(adjusted_log_data_matrix)

# Perform PCA
pca_results <- prcomp(t(adjusted_log_data_matrix), center = TRUE, scale. = FALSE)

# Visualize PCA results
plot(pca_results$x, col = sample_colors)
legend("topright", legend = c("Luminal_A", "Luminal_B"), col = unique(sample_colors), pch = 1, cex = 0.8)

# Create Elbow Graph
variance_explained <- pca_results$sdev^2 / sum(pca_results$sdev^2) * 100
plot(1:length(variance_explained), variance_explained, type = "b", pch = 19)
abline(h = 1, col = "red", lty = 2)  # Optional threshold for selecting PCs

# Perform t-SNE
set.seed(123)  # Set seed for reproducibility
tsne_results <- Rtsne(pca_results$x[, 1:15], perplexity = 20)  # Adjust perplexity if needed

# Visualize t-SNE results
plot(tsne_results$Y, col = sample_colors)
legend("topleft", legend = c("Luminal_A", "Luminal_B"), col = unique(sample_colors), pch = 1, cex = 0.8)
