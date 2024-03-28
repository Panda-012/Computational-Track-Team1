# Load necessary libraries
library(data.table)
library(jsonlite) 

# Set the current working directory to the script directory
setwd(normalizePath(dirname(rstudioapi::getSourceEditorContext()$path)))
getwd()

# Read the metadata file content and create an empty column for case_id
metadata <- fromJSON(readLines("../RNA-Seq_Metadata.json"))
metadata <- cbind(metadata[, 1:3], case_id = NA, metadata[, 4:ncol(metadata)])

# Loop through each entry in the metadata to find its case id
for (i in 1:nrow(metadata)) {
  matching_record <- metadata[grep(metadata$file_name[[i]], metadata$file_name), ]
  metadata$case_id[i] <- matching_record$associated_entities[[1]]$case_id
}

# Prepare the replicates samples for later name handling
replicates <- data.frame (
  id = unique(metadata$case_id[which(duplicated(metadata$case_id))]),
  count = rep(0, length(unique(metadata$case_id[which(duplicated(metadata$case_id))])))
  )

# Check if the result file exists and delete it if it does
result_file = "RNA-Seq_Data.tsv"
if (file.exists(result_file)) {
  file.remove(result_file)
}

# List all TSV files in the current directory
file_list <- list.files(pattern = "\\.tsv", recursive = TRUE)

# Initialize an empty character vector to store all identifiers
max_rownames <- character()

# Loop through each file to collect all identifiers
for (file in file_list) {
  # Read the current file into a data.table and skip unnecessary rows
  current_df <- fread(file, skip = 6, header = FALSE, fill = TRUE, select = c(2,4), col.names = c("gene_name", "expression"))
  
  # Add unique identifiers to the set
  max_rownames <- union(max_rownames, as.character(current_df$gene_name))
}

# Create a dataframe from the stored data
result_df <- data.frame("gene_name" = max_rownames)
counter <- 1

# Loop through each file
for (file in file_list) {
  # Read the current file into a data.table and aggregate duplicate gene names
  current_df <- fread(file, skip = 6, header = FALSE, fill = TRUE, select = c(2,4), col.names = c("gene_name", "expression"))
  current_df <- current_df[, .(expression = mean(expression)), by = gene_name]
  extracted_filename <- gsub("^.*/([^/]+)$", "\\1", file)
  
  # Extract the colname and the values column
  current_colname <- as.character(metadata$case_id[metadata$file_name == extracted_filename])
  current_rownames <- current_df$gene_name
  current_col <- current_df$expression
  
  # Initialize mapped_values with NA for all rows
  mapped_values <- rep(NA, length(max_rownames))
  
  # Find matching indices between current_rownames and max_rownames
  match_indices <- match(current_rownames, max_rownames)
  
  # Update mapped_values with matching values
  mapped_values[match_indices[!is.na(match_indices)]] <- current_col[!is.na(match_indices)]

  # Store mapped values for this file
  if (current_colname %in% replicates$id){
    replicates$count [replicates$id == current_colname] <- replicates$count [replicates$id == current_colname] + 1
    current_colname <- paste0 (current_colname, "X", replicates$count [replicates$id == current_colname])
  }
  result_df[[current_colname]] <- mapped_values
  
  # Print a progress message
  cat(paste0("[", counter, "]"), file, "\n")
  counter <- counter + 1
}

# Write the result dataframe to a TSV file without quoting string values
write.table(result_df, file = result_file, sep = "\t", row.names = FALSE, quote = FALSE)
