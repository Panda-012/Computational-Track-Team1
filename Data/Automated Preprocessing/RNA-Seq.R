# Load necessary libraries
library(readr)
library(limma)
library(edgeR)
library(biomaRt)
library(dplyr)

# Function to preprocess count matrix: loads data, maps Ensembl IDs, aggregates duplicates
preprocess_count_matrix <- function(file_path, dataset = "hsapiens_gene_ensembl") {
  # Determine file extension to decide on read function
  file_extension <- tools::file_ext(file_path)
  
  # Load the RNA-Seq count matrix based on the file type
  if (file_extension == "csv") {
    count_matrix <- read.csv(file_path, stringsAsFactors = FALSE)
  } else if (file_extension == "tsv") {
    count_matrix <- read.delim(file_path, stringsAsFactors = FALSE)
  } else {
    stop("Unsupported file type. Please provide a CSV or TSV file.")
  }
  
  # Clean the Ensembl IDs by removing version numbers
  ensembl_ids <- sub("\\..*$", "", count_matrix[, 1]) 
  
  # Connect to the Ensembl database
  ensembl <- biomaRt::useMart("ensembl", dataset = dataset)
  
  # Fetch gene symbols for the cleaned Ensembl IDs
  gene_symbols <- biomaRt::getBM(attributes = c('ensembl_gene_id', 'external_gene_name'), 
                                 filters = 'ensembl_gene_id', 
                                 values = ensembl_ids, 
                                 mart = ensembl)
  
  # Replace missing gene symbols with NA for easier filtering
  gene_symbols$external_gene_name[gene_symbols$external_gene_name == ""] <- NA
  
  # Map Ensembl IDs to gene symbols
  id_to_symbol_map <- setNames(gene_symbols$external_gene_name, gene_symbols$ensembl_gene_id)
  count_matrix[, 1] <- id_to_symbol_map[ensembl_ids]
  
  # Remove rows with missing gene symbols
  count_matrix <- count_matrix[!(is.na(count_matrix[, 1])), ]
  
  # Aggregate duplicate gene symbols by summing their counts
  if (anyDuplicated(count_matrix[, 1])) {
    count_matrix <- count_matrix %>%
      group_by(Gene = count_matrix[, 1]) %>%
      summarise(across(where(is.numeric), sum, na.rm = TRUE), .groups = 'drop') %>%
      as.data.frame()
  }
  
  # Convert to a matrix for further analysis
  rownames(count_matrix) <- as.character(count_matrix[, 1])
  count_matrix <- count_matrix[, -1]
  count_matrix <- as.matrix(count_matrix)
  
  # apply normalization and outliers handling function
  adjusted_log_count_matrix <- normalize_and_adjust_outliers(count_matrix)
  return (adjusted_log_count_matrix)
}

# Function to normalize count matrix and adjust for outliers
normalize_and_adjust_outliers <- function(count_matrix) {
  # Create a DGEList object for normalization
  dge <- DGEList(counts = count_matrix)
  
  # Normalize the data using the TMM method
  dge <- calcNormFactors(dge, method = "TMM")
  
  # Log-transform the counts for variance stabilization
  log_count_matrix <- cpm(dge, log = TRUE, prior.count = 1)
  
  # Function to adjust outliers within a single vector of gene expression values
  adjust_outliers_directly <- function(gene_expression) {
    # Constants for outlier detection
    k <- 3
    gene_mad <- mad(gene_expression, constant = 1)
    gene_median <- median(gene_expression)
    lb <- gene_median - k * gene_mad
    ub <- gene_median + k * gene_mad
    
    # Adjust outliers to the nearest value within the non-outlier range
    gene_expression[gene_expression < lb] <- lb
    gene_expression[gene_expression > ub] <- ub
    return(gene_expression)
  }
  
  # Apply the outlier adjustment across all genes
  adjusted_log_count_matrix <- t(apply(log_count_matrix, 1, adjust_outliers_directly))
  
  return(adjusted_log_count_matrix)
}

# Apply the functions on your datasets
setwd(normalizePath(dirname(rstudioapi::getSourceEditorContext()$path)))
SKCM_RNA_count_matrix <- preprocess_count_matrix("SKCM_RNA.csv")
