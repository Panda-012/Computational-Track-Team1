# Load necessary libraries
library(limma)
library(Rtsne)
library(edgeR)
library(e1071)
library(impute)
library(preprocessCore)
library(minfi) # For DNA methylation data handling
library(dplyr) # For data manipulation

# Set working directory and load the data
setwd(normalizePath(dirname(rstudioapi::getSourceEditorContext()$path)))
getwd()

# Load DNA methylation count matrix
count_matrix <- readRDS("DNA_Methylation_Data.rds")
rownames(count_matrix) <- as.character(count_matrix[,1])
count_matrix <- count_matrix[,-1]
count_matrix <- as.matrix(count_matrix)  # Convert to numeric matrix
sample_colors <- ifelse(startsWith(colnames(count_matrix), "LumB"), "red", "orange")
Environment_Data <- ls()

# Check for duplicate CpG Sites
sum(duplicated(rownames(count_matrix)))

# Identify missing values
sum(count_matrix == 0, na.rm = TRUE)
sum(is.na(count_matrix))

# Assess quality before processing
affy::plotDensity(count_matrix)
plot(density(apply(count_matrix, 2, mean, na.rm = TRUE)),main="omics",cex.axis=0.5) #1 rows, 2 columns
plotMDS(count_matrix)

# Assess the missing data
missing_per_site <- colSums(is.na(count_matrix)) / nrow(count_matrix)
missing_per_sample <- rowSums(is.na(count_matrix)) / ncol(count_matrix)

# Dynamic thresholding based on the distribution of missing values
site_missing_threshold <- mean(missing_per_site) + sd(missing_per_site)
sample_missing_threshold <- mean(missing_per_sample) + sd(missing_per_sample)

# Filter sites and samples based on the calculated thresholds
filtered_matrix <- count_matrix[rowSums(is.na(count_matrix)) / ncol(count_matrix) < sample_missing_threshold, ]
filtered_matrix <- filtered_matrix[, colSums(is.na(filtered_matrix)) / nrow(filtered_matrix) < site_missing_threshold]

# Assess quality after initial filtration
affy::plotDensity(filtered_matrix)
plot(density(apply(filtered_matrix, 2, mean, na.rm = TRUE)),main="omics",cex.axis=0.5) #1 rows, 2 columns
plotMDS(filtered_matrix)

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

# Assess quality after filtration
affy::plotDensity(filtered_data_matrix)
plot(density(apply(filtered_data_matrix, 2, mean, na.rm = TRUE)),main="omics",cex.axis=0.5) #1 rows, 2 columns
plotMDS(filtered_data_matrix)

# Now, impute missing values in the filtered matrix
k <- min(10, ncol(filtered_data_matrix) - 1) # Set neighbors for KNN imputation
imputed_matrix <- impute.knn(as.matrix(filtered_data_matrix), k = k)$data

# Assess quality after imputation
affy::plotDensity(imputed_matrix)
plot(density(apply(imputed_matrix, 2, mean, na.rm = TRUE)),main="omics",cex.axis=0.5) #1 rows, 2 columns
plotMDS(imputed_matrix)

# After filtering, clean up the environment from temporary variables
rm(list=setdiff(ls(), c(Environment_Data, "imputed_matrix", "filtered_matrix","filtered_data_matrix")))

# Apply Quantile Normalization
normalized_data_matrix <- as.data.frame(normalize.quantiles(as.matrix(imputed_matrix)))

# Restore rownames and colnames
rownames(normalized_data_matrix) <- rownames(imputed_matrix)
colnames(normalized_data_matrix) <- colnames(imputed_matrix)

# Convert back to matrix
normalized_data_matrix <- as.matrix(normalized_data_matrix)

# Reassess quality after Quantile Normalization
affy::plotDensity(normalized_data_matrix)
plot(density(apply(normalized_data_matrix, 2, mean, na.rm = TRUE)),main="omics",cex.axis=0.5) #1 rows, 2 columns
plotMDS(normalized_data_matrix)

# Function to adjust outliers directly within a single vector of CpG Methylation values
Environment_Data <- ls()
adjust_outliers_dynamically <- function(data) {
  calculate_outlier_bounds <- function(data, k) {
    Q1 <- quantile(data, 0.25, na.rm = TRUE)
    Q3 <- quantile(data, 0.75, na.rm = TRUE)
    IQR <- Q3 - Q1
    lower_bound <- Q1 - k * IQR
    upper_bound <- Q3 + k * IQR
    return(c(lower_bound, upper_bound))
  }
  
  # Initial outlier detection with standard k=1.5
  k <- 1.5
  bounds <- calculate_outlier_bounds(data, k)
  outlier_proportion <- mean(data < bounds[1] | data > bounds[2], na.rm = TRUE)
  
  # Adjust k based on the outlier proportion
  target_proportion <- 0.05 # Targeting 5% outliers
  while (outlier_proportion > target_proportion && k < 10) {
    k <- k + 0.5
    bounds <- calculate_outlier_bounds(data, k)
    outlier_proportion <- mean(data < bounds[1] | data > bounds[2], na.rm = TRUE)
  }
  
  # Adjust the data points considered as outliers
  adjusted_data <- ifelse(data < bounds[1], bounds[1], ifelse(data > bounds[2], bounds[2], data))
  return(adjusted_data)
}

# Apply the dynamic adjustment of outliers across all CpG sites
adjusted_data_matrix <- t(apply(normalized_data_matrix, 1, adjust_outliers_dynamically))
rm(list=setdiff(ls(), c(Environment_Data, "adjust_outliers_dynamically", "adjusted_data_matrix")))

# Reassess quality after handling outliers
affy::plotDensity(adjusted_data_matrix)
plot(density(apply(adjusted_data_matrix, 2, mean, na.rm = TRUE)),main="omics",cex.axis=0.5) #1 rows, 2 columns
plotMDS(adjusted_data_matrix)

# Perform PCA
pca_results <- prcomp(t(adjusted_data_matrix), center = TRUE, scale. = TRUE)

# Visualize PCA results
plot(pca_results$x, col = sample_colors)
legend("topright", legend = c("Luminal_A", "Luminal_B"), col = unique(sample_colors), pch = 1, cex = 0.8)

# Create Elbow Graph
variance_explained <- pca_results$sdev^2 / sum(pca_results$sdev^2) * 100
plot(1:length(variance_explained), variance_explained, type = "b", pch = 19)
abline(h = 2.3, col = "red", lty = 2)  # Optional threshold for selecting PCs

# Perform t-SNE
set.seed(123)  # Set seed for reproducibility
tsne_results <- Rtsne(pca_results$x[, 1:15], perplexity = 10)  # Adjust perplexity if needed

# Visualize t-SNE results
plot(tsne_results$Y, col = sample_colors)
legend("topleft", legend = c("Luminal_A", "Luminal_B"), col = unique(sample_colors), pch = 1, cex = 0.8)
