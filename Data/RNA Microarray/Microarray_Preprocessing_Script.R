# Load necessary libraries
library(limma)
library(DescTools)
library(preprocessCore)
library(outliers)
library(genefilter)
library(Rtsne)

# Set working directory
setwd(normalizePath(dirname(rstudioapi::getSourceEditorContext()$path)))
getwd()

# Read TSV file
data <- read.delim("Microarray_Data.tsv")
sample_colors <- ifelse(startsWith(colnames(data), "Normal"), "blue", 
                        ifelse(startsWith(colnames(data), "Luminal_B"), "red", "orange"))

# check for duplicate gene names
any(duplicated(rownames(data)))

# Identify missing values
any(is.na(data))
any(data == 0)

# Identify outliers
data_matrix <- as.matrix(data)  # Convert to numeric matrix
length(boxplot.stats(data_matrix)$out)

# Visualize the data distribution
boxplot(data_matrix)
affy::plotDensity(data_matrix)
plot(density(apply(data_matrix, 2, mean, na.rm = TRUE)),main="omics",cex.axis=0.5) #1 rows, 2 columns
plotMDS(data_matrix)

# Check for zeros and add a small constant if necessary
data_matrix[data_matrix == 0] <- min(data_matrix[data_matrix > 0]) / 100

# Apply log2 transformation to stabilize variance
log_data_matrix <- log2(data_matrix)
length(boxplot.stats(log_data_matrix)$out)

# Outlier detection using multiple methods
consensus_outliers <- apply(log_data_matrix, 1, function(x) {
  grubbs_out <- tryCatch(outliers::grubbs.test(x, opposite = TRUE)$p.value, error = function(e) NA)
  gesd_out <- tryCatch(genefilter::rowQ(x, probs = c(0.25, 0.75), na.rm = TRUE), error = function(e) rep(NA, 2))
  iqr_out <- IQR(x)
  
  # Ensure 'gesd_out' has non-NA values before proceeding
  if (!is.na(grubbs_out) && !any(is.na(gesd_out)) && grubbs_out < 0.05) {
    # Get the 75th percentile value from 'gesd_out'
    upper_gesd <- gesd_out[2]
    # Compare the maximum value to the upper threshold (75th percentile + 1.5 * IQR)
    if(max(x) > (upper_gesd + 1.5 * iqr_out)) {
      return(max(x))  # Return the outlier value
    }
  }
  return(NA)  # If there's no consensus, or an NA was encountered, return NA
})

# Apply Winsorization to the log-transformed data using consensus outliers
data_matrix_winsorized <- t(apply(log_data_matrix, 1, function(x) {
  # If there is a consensus outlier for the row, perform Winsorization
  if (any(!is.na(consensus_outliers))) {
    lower_bound <- quantile(x, probs = 0.05, na.rm = TRUE)
    upper_bound <- quantile(x, probs = 0.95, na.rm = TRUE)
    return(Winsorize(x, minval = lower_bound, maxval = upper_bound))
  } else {
    return(x)  # If no outliers, return the row as-is
  }
}))
length(boxplot.stats(data_matrix_winsorized)$out)

# Reassess quality after Winsorizing
boxplot(data_matrix_winsorized)
affy::plotDensity(data_matrix_winsorized)
plot(density(apply(data_matrix_winsorized, 2, mean, na.rm = TRUE)),main="omics",cex.axis=0.5) #1 rows, 2 columns
plotMDS(data_matrix_winsorized)

# Apply quantile normalization to your expression matrix
data_matrix_normalized <- normalize.quantiles(data_matrix_winsorized)
rownames(data_matrix_normalized) <- rownames(data_matrix_winsorized)
colnames(data_matrix_normalized) <- colnames(data_matrix_winsorized)

# Reassess quality after quantile normalization
boxplot(data_matrix_winsorized)
affy::plotDensity(data_matrix_normalized)
plot(density(apply(data_matrix_normalized, 2, mean, na.rm = TRUE)),main="omics",cex.axis=0.5) #1 rows, 2 columns
plotMDS(data_matrix_winsorized)

# Perform PCA
pca_results <- prcomp(t(data_matrix_normalized), center = TRUE, scale. = TRUE)

# Visualize PCA results
plot(pca_results$x, col = sample_colors)
legend("bottomright", legend = c( "Luminal_B", "Luminal_A", "Normal"), col = unique(sample_colors), pch = 1, cex = 0.8)

# Create Elbow Graph
variance_explained <- pca_results$sdev^2 / sum(pca_results$sdev^2) * 100
plot(1:length(variance_explained), variance_explained, type = "b", pch = 19)
abline(h = 2, col = "red", lty = 2)  # Optional threshold for selecting PCs

# Perform t-SNE
set.seed(123)  # Set seed for reproducibility
tsne_results <- Rtsne(pca_results$x[, 1:10], perplexity = 5)  # Adjust perplexity if needed

# Visualize t-SNE results
plot(tsne_results$Y, col = sample_colors)
legend("bottomleft", legend = c( "Luminal_B", "Luminal_A", "Normal"), col = unique(sample_colors), pch = 1, cex = 0.8)
