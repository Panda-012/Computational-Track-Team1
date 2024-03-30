# Load necessary libraries
library(readr)
library(limma)
library(Rtsne)
library(edgeR)
library(glmmTMB)
library(e1071)
library(data.table)
library(impute)
library(ggplot2)
library(minfi) # For DNA methylation data handling
library(reshape2) # For melt function
library(dplyr) # For data manipulation
library(jsonlite)  # For loading metadata and clinical data

# Set working directory and load the data
setwd(normalizePath(dirname(rstudioapi::getSourceEditorContext()$path)))
getwd()

# Load DNA methylation count matrix
load("DNA_Methylation_Data.RData")
rownames(count_matrix) <- as.character(count_matrix[,1])
count_matrix <- count_matrix[,-1]
count_matrix <- as.matrix(count_matrix)  # Convert to numeric matrix
sample_colors <- ifelse(startsWith(colnames(count_matrix), "LumB"), "red", "orange")
Environment_Data <- ls()

# Check for duplicate CpG Sites
sum(duplicated(rownames(count_matrix)))

# Identify missing values
sum(is.na(count_matrix))
sum(count_matrix == 0, na.rm = TRUE)

# Assess the missing data
missing_per_site <- colSums(is.na(count_matrix)) / nrow(count_matrix)
missing_per_sample <- rowSums(is.na(count_matrix)) / ncol(count_matrix)

# Dynamic thresholding based on the distribution of missing values
site_missing_threshold <- mean(missing_per_site) + sd(missing_per_site)
sample_missing_threshold <- mean(missing_per_sample) + sd(missing_per_sample)

# Filter sites and samples based on the calculated thresholds
filtered_matrix <- count_matrix[rowSums(is.na(count_matrix)) / ncol(count_matrix) < sample_missing_threshold, ]
filtered_matrix <- filtered_matrix[, colSums(is.na(filtered_matrix)) / nrow(filtered_matrix) < site_missing_threshold]

# Calculate median methylation and MAD before imputation, for filtering CpG Sites
median_methylation <- apply(filtered_matrix, 1, median, na.rm = TRUE)
mad_methylation <- apply(filtered_matrix, 1, function(x) mad(x, constant = 1, na.rm = TRUE))

# Dynamically set mad_multiplier based on skewness
data_skewness <- skewness(log1p(rowMedians(filtered_matrix, na.rm = TRUE)))
mad_multiplier <- ifelse(data_skewness <= 1, 1.5, ifelse(data_skewness <= 2, 2, 2.5))

# Combine median and MAD for a robust measure of DNA Methylation variability
threshold_methylation <- median_methylation + (mad_multiplier * mad_methylation)

# Dynamically set samples_threshold based on overall methylation density
overall_non_zero_proportion <- mean(rowSums(filtered_matrix > 0, na.rm = TRUE)) / ncol(filtered_matrix)
samples_threshold <- 0.05 + (0.15 * (1 - overall_non_zero_proportion))

# Identify and keep CpG Sites based on the threshold
keep_sites <- rowSums(filtered_matrix >= threshold_methylation, na.rm = TRUE) >= (samples_threshold * ncol(filtered_matrix))
filtered_data_matrix <- filtered_matrix[keep_sites, ]

# Now, impute missing values in the filtered matrix
k <- min(10, ncol(filtered_data_matrix) - 1) # Set neighbors for KNN imputation
imputed_matrix <- impute.knn(as.matrix(filtered_data_matrix), k = k)$data

# Assess quality before processing
affy::plotDensity(imputed_matrix)
plot(density(apply(imputed_matrix, 2, mean, na.rm = TRUE)),main="omics",cex.axis=0.5) #1 rows, 2 columns
plotMDS(imputed_matrix)

# After filtering, clean up the environment from temporary variables
rm(list=setdiff(ls(), c(Environment_Data, "imputed_matrix")))

# Identify outliers (visually and statistically)
boxplot(imputed_matrix)
length(boxplot.stats(imputed_matrix)$out)

# Create dge and normalize the data using TMM
dge <- DGEList(counts = imputed_matrix)
dge <- calcNormFactors(dge, method = "TMM")  # Perform TMM normalization within DGEList

# Apply log2 transformation to stabilize variance
log_data_matrix <- cpm(dge, log = TRUE, prior.count = 1)
length(boxplot.stats(log_data_matrix)$out)

# Reassess quality after log transformation
affy::plotDensity(log_data_matrix)
plot(density(apply(log_data_matrix, 2, mean, na.rm = TRUE)),main="omics",cex.axis=0.5) #1 rows, 2 columns
plotMDS(log_data_matrix)

# Function to adjust outliers directly within a single vector of CpG Methylation values
Environment_Data <- ls()
adjust_outliers_directly <- function(CpG_Methylation) {
  # Constants for outlier detection
  k <- 3 # Adjust k as needed for stricter or looser criteria
  
  # Compute lower and upper bounds for this CpG site
  CpG_mad <- mad(CpG_Methylation, constant = 1)
  CpG_median <- median(CpG_Methylation)
  lb <- CpG_median - k * CpG_mad
  ub <- CpG_median + k * CpG_mad
  
  # Adjust outliers to the nearest value within the non-outlier range
  CpG_Methylation[CpG_Methylation < lb] <- lb
  CpG_Methylation[CpG_Methylation > ub] <- ub
  
  return(CpG_Methylation)
}

# Apply the direct adjustment of outliers across all CpG sites
adjusted_log_data_matrix <- t(apply(log_data_matrix, 1, adjust_outliers_directly))
length(boxplot.stats(adjusted_log_data_matrix)$out)
rm(list=setdiff(ls(), c(Environment_Data, "adjust_outliers_directly", "adjusted_log_data_matrix")))

# Reassess quality after handling outliers
affy::plotDensity(adjusted_log_data_matrix)
plot(density(apply(adjusted_log_data_matrix, 2, mean, na.rm = TRUE)),main="omics",cex.axis=0.5) #1 rows, 2 columns
plotMDS(adjusted_log_data_matrix)

# Perform PCA
pca_results <- prcomp(t(adjusted_log_data_matrix), center = TRUE, scale. = TRUE)

# Visualize PCA results
plot(pca_results$x, col = sample_colors)
legend("topright", legend = c("Luminal_A", "Luminal_B"), col = unique(sample_colors), pch = 1, cex = 0.8)

# Create Elbow Graph
variance_explained <- pca_results$sdev^2 / sum(pca_results$sdev^2) * 100
plot(1:length(variance_explained), variance_explained, type = "b", pch = 19)
abline(h = 2.5, col = "red", lty = 2)  # Optional threshold for selecting PCs

# Perform t-SNE
set.seed(123)  # Set seed for reproducibility
tsne_results <- Rtsne(pca_results$x[, 1:10], perplexity = 10)  # Adjust perplexity if needed

# Visualize t-SNE results
plot(tsne_results$Y, col = sample_colors)
legend("topleft", legend = c("Luminal_A", "Luminal_B"), col = unique(sample_colors), pch = 1, cex = 0.8)
