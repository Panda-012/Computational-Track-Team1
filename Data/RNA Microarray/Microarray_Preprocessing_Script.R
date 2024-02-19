# Load necessary libraries
library(limma)
library(DescTools)
library(preprocessCore)
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

# Handle outliers by applying Winsorizing to the expression matrix
data_matrix_winsorized <- Winsorize(data_matrix, probs = c(0.01, 0.99))
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
