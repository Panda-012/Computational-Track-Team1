# Load necessary libraries
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(minfi)
library(impute)
library(dplyr)
library(e1071)
library(preprocessCore)

# Main function to preprocess DNA methylation data
preprocess_methylation_matrix <- function(file_path) {
  # Load DNA methylation beta value matrix
  beta_matrix <- load_beta_value_matrix(file_path)
  
  # Check for and handle duplicates and missing values
  imputed_filtered_beta_matrix <- filter_and_impute(beta_matrix)
  
  # Apply quantile normalization directly to the beta value matrix
  normalized_beta_matrix <- normalize(imputed_filtered_beta_matrix)
  
  # Adjust for outliers
  final_beta_matrix <- t(apply(normalized_beta_matrix, 1, adjust_outliers))
  
  return(final_beta_matrix)
}

# Function to load DNA methylation data
load_beta_value_matrix <- function(file_path) {
  file_extension <- tools::file_ext(file_path)
  
  if (file_extension == "csv") {
    beta_matrix <- read.csv(file_path, stringsAsFactors = FALSE, check.names = FALSE)
  } else if (file_extension == "tsv" || file_extension == "txt") {
    beta_matrix <- read.delim(file_path, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    stop("Unsupported file type. Please provide a CSV or TSV file.")
  }
  
  # Handling duplicates by averaging beta values
  colnames(beta_matrix)[1] <- "CpG_Site"
  if (anyDuplicated(beta_matrix[, "CpG_Site"])) {
    beta_matrix <- beta_matrix %>%
      group_by(CpG_Site) %>%
      summarise(across(where(is.numeric), mean, na.rm = TRUE), .groups = 'drop')
  }
  
  # Convert the first column into row names
  beta_matrix <- as.data.frame(beta_matrix)
  row.names(beta_matrix) <- beta_matrix$CpG_Site
  beta_matrix <- as.matrix(beta_matrix[, -1])
  
  return(beta_matrix)
}

# Function to filter CpG Sites with low variability and handle missing values
filter_and_impute <- function(beta_matrix) {
  # Assess the missing data
  missing_per_site <- colSums(is.na(beta_matrix)) / nrow(beta_matrix)
  missing_per_sample <- rowSums(is.na(beta_matrix)) / ncol(beta_matrix)
  
  # Dynamic thresholding based on the distribution of missing values
  site_missing_threshold <- mean(missing_per_site) + sd(missing_per_site)
  sample_missing_threshold <- mean(missing_per_sample) + sd(missing_per_sample)
  
  # Filter sites and samples based on the calculated thresholds
  filtered_beta_matrix <- beta_matrix[rowSums(is.na(beta_matrix)) / ncol(beta_matrix) < sample_missing_threshold, ]
  filtered_beta_matrix <- filtered_beta_matrix[, colSums(is.na(filtered_beta_matrix)) / nrow(filtered_beta_matrix) < site_missing_threshold]
  
  # Calculate median methylation and MAD before imputation, for filtering CpG Sites
  median_methylation <- apply(filtered_beta_matrix, 1, median, na.rm = TRUE)
  mad_methylation <- apply(filtered_beta_matrix, 1, function(x) mad(x, constant = 1, na.rm = TRUE))
  
  # Dynamically set mad_multiplier based on skewness
  data_skewness <- skewness(log1p(rowMedians(filtered_beta_matrix, na.rm = TRUE)))
  mad_multiplier <- ifelse(data_skewness <= 1, 1.5, ifelse(data_skewness <= 2, 2, 2.5))
  
  # Combine median and MAD for a robust measure of DNA Methylation variability
  threshold_methylation <- median_methylation + (mad_multiplier * mad_methylation)
  
  # Dynamically set samples_threshold based on overall methylation density
  overall_non_zero_proportion <- mean(rowSums(filtered_beta_matrix > 0, na.rm = TRUE)) / ncol(filtered_beta_matrix)
  samples_threshold <- 0.05 + (0.15 * (1 - overall_non_zero_proportion))
  
  # Identify and keep CpG Sites based on the threshold
  keep_sites <- rowSums(filtered_beta_matrix >= threshold_methylation, na.rm = TRUE) >= (samples_threshold * ncol(filtered_beta_matrix))
  filtered_beta_matrix <- filtered_beta_matrix[keep_sites, ]
  
  # Now, impute missing values in the filtered matrix
  k <- min(10, ncol(filtered_beta_matrix) - 1) # Set neighbors for KNN imputation
  imputed_filtered_beta_matrix <- impute.knn(as.matrix(filtered_beta_matrix), k = k)$data
  
  return(imputed_filtered_beta_matrix)
}

# Function to perform Quantile Normalization
normalize <- function(imputed_filtered_beta_matrix){
  # Apply Quantile Normalization
  normalized_beta_matrix <- as.data.frame(normalize.quantiles(as.matrix(imputed_filtered_beta_matrix)))
  
  # Restore rownames and colnames
  rownames(normalized_beta_matrix) <- rownames(imputed_filtered_beta_matrix)
  colnames(normalized_beta_matrix) <- colnames(imputed_filtered_beta_matrix)
  
  return(as.matrix(normalized_beta_matrix))
}

# Function to adjust for outliers
adjust_outliers <- function(normalized_beta_matrix) {
  calculate_outlier_bounds <- function(normalized_beta_matrix, k) {
    Q1 <- quantile(normalized_beta_matrix, 0.25, na.rm = TRUE)
    Q3 <- quantile(normalized_beta_matrix, 0.75, na.rm = TRUE)
    IQR <- Q3 - Q1
    lower_bound <- Q1 - k * IQR
    upper_bound <- Q3 + k * IQR
    return(c(lower_bound, upper_bound))
  }
  
  # Initial outlier detection with standard k=1.5
  k <- 1.5
  bounds <- calculate_outlier_bounds(normalized_beta_matrix, k)
  outlier_proportion <- mean(normalized_beta_matrix < bounds[1] | normalized_beta_matrix > bounds[2], na.rm = TRUE)
  
  # Adjust k based on the outlier proportion
  target_proportion <- 0.05 # Targeting 5% outliers
  while (outlier_proportion > target_proportion && k < 10) {
    k <- k + 0.5
    bounds <- calculate_outlier_bounds(normalized_beta_matrix, k)
    outlier_proportion <- mean(normalized_beta_matrix < bounds[1] | normalized_beta_matrix > bounds[2], na.rm = TRUE)
  }
  
  # Adjust the data points considered as outliers
  adjusted_matrix <- ifelse(normalized_beta_matrix < bounds[1], bounds[1], ifelse(normalized_beta_matrix > bounds[2], bounds[2], normalized_beta_matrix))
  
  return(adjusted_matrix)
}

# Apply the functions on your datasets
setwd(normalizePath(dirname(rstudioapi::getSourceEditorContext()$path)))
sample_DNA_Methylation <- preprocess_methylation_matrix("MethylationSample.csv")
