# This script should be adjusted easily to work with any file types
# by simply altering the file extension. Further modifications may 
# be required regarding the tabular data structure and desired outcome,
# including the output rows and columns.

# Set the current working directory to the script directory
setwd(normalizePath(dirname(rstudioapi::getSourceEditorContext()$path)))
getwd()

# Check if the result file exists and delete it if it does
result_file = "Sample_DNA_Methylation_Data.txt"
if (file.exists(result_file)) {
  file.remove(result_file)
}

# List all text files in the current directory
file_list <- list.files(pattern = "\\.txt", recursive = TRUE)

# Initialize an empty set to store all identifiers
max_rownames <- c()

# Loop through each file to collect all identifiers
for (file in file_list) {
  current_df <- read.table(file, header = FALSE, sep = "\t", fill = TRUE)
  
  # Collect identifiers from the current file
  current_rownames <- current_df[, 1]
  
  # Add unique identifiers to the set
  max_rownames <- union(max_rownames, current_rownames)
}

# Create a dataframe from the stored data
result_df <- data.frame("Methylation_ID" = max_rownames)

# Loop through each file
for (file in file_list) {
  # Read the current file into a dataframe without considering the first row as header
  current_df <- read.table(file, header = FALSE, sep = "\t", fill = TRUE)
  
  # Extract the colname and the values column
  current_rownames <- current_df[, 1]
  current_colname <- gsub("/.*", "", file)
  current_col <- current_df[, 2]
  
  # Initialize mapped_values with NA for all rows
  mapped_values <- rep(NA, length(max_rownames))
  
  # Find matching indices between current_rownames and max_rownames
  match_indices <- match(current_rownames, max_rownames)
  
  # Update mapped_values with matching values
  mapped_values[match_indices[!is.na(match_indices)]] <- current_col[!is.na(match_indices)]

  
  # Store mapped values for this file
  result_df[[current_colname]] <- mapped_values
}

# Write the result dataframe to a text file without quoting string values
write.table(result_df, file = result_file, sep = "\t", row.names = FALSE, quote = FALSE)