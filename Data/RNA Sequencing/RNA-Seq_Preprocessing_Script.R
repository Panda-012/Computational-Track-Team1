# Load necessary libraries
library(readr)
library(limma)
library(edgeR)
library(e1071)
library(data.table)
library(Rtsne)
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
# load("../luminal_samples.RData")

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

# Normalize the data using the TMM method from edgeR
dge <- DGEList(counts = data_matrix)
dge <- calcNormFactors(dge, method = "TMM")
log_data_matrix <- cpm(dge, log = TRUE, prior.count = 1)

# Calculate median expression level across genes and median absolute deviation (MAD) for each gene
median_expression <- apply(log_data_matrix, 1, median)
mad_expression <- apply(log_data_matrix, 1, mad)

# Dynamically set mad_multiplier based on skewness
data_skewness <- skewness(median_expression)
mad_multiplier <- ifelse(data_skewness <= 1, 1.5, ifelse(data_skewness <= 2, 2, 2.5))

# Combine median and MAD for a robust measure of gene expression variability
threshold_expression <- median_expression + (mad_multiplier * mad_expression)

# Filter genes with low expression variability
keep_genes <- mad_expression > quantile(mad_expression, 0.25)  # Keep genes with highest variability
filtered_data_matrix <- log_data_matrix[keep_genes, ]

# After filtering, clean up the environment from temporary variables
rm(list=setdiff(ls(), c("LumA", "LumB", "Basal", "Her2", "Normal", Environment_Data, "data_matrix", "filtered_data_matrix", "sample_colors")))
Environment_Data <- ls()

# Reassess quality after filtration
affy::plotDensity(filtered_data_matrix)
plot(density(apply(filtered_data_matrix, 2, mean, na.rm = TRUE)),main="omics",cex.axis=0.5) #1 rows, 2 columns
plotMDS(filtered_data_matrix)

# Function to adjust outliers directly within a single vector of gene expression values
adjust_outliers <- function(filtered_data_matrix) {
  calculate_outlier_bounds <- function(filtered_data_matrix, k) {
    Q1 <- quantile(filtered_data_matrix, 0.25, na.rm = TRUE)
    Q3 <- quantile(filtered_data_matrix, 0.75, na.rm = TRUE)
    IQR <- Q3 - Q1
    lower_bound <- Q1 - k * IQR
    upper_bound <- Q3 + k * IQR
    return(c(lower_bound, upper_bound))
  }
  
  # Initial outlier detection with standard k=1.5
  k <- 1.5
  bounds <- calculate_outlier_bounds(filtered_data_matrix, k)
  outlier_proportion <- mean(filtered_data_matrix < bounds[1] | filtered_data_matrix > bounds[2], na.rm = TRUE)
  
  # Adjust k based on the outlier proportion
  target_proportion <- 0.05 # Targeting 5% outliers
  while (outlier_proportion > target_proportion && k < 10) {
    k <- k + 0.5
    bounds <- calculate_outlier_bounds(filtered_data_matrix, k)
    outlier_proportion <- mean(filtered_data_matrix < bounds[1] | filtered_data_matrix > bounds[2], na.rm = TRUE)
  }
  
  # Adjust the data points considered as outliers
  adjusted_matrix <- ifelse(filtered_data_matrix < bounds[1], bounds[1], ifelse(filtered_data_matrix > bounds[2], bounds[2], filtered_data_matrix))
  
  return(adjusted_matrix)
}

# Apply the direct adjustment of outliers across all genes
adjusted_data_matrix <- t(apply(filtered_data_matrix, 1, adjust_outliers))
rm(list=setdiff(ls(), c(Environment_Data, "adjust_outliers_directly", "adjusted_data_matrix")))

# Reassess quality after handling outliers
affy::plotDensity(adjusted_data_matrix)
plot(density(apply(adjusted_data_matrix, 2, mean, na.rm = TRUE)),main="omics",cex.axis=0.5) #1 rows, 2 columns
plotMDS(adjusted_data_matrix)

# Perform PCA
pca_results <- prcomp(t(adjusted_data_matrix), center = TRUE, scale. = FALSE)

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
