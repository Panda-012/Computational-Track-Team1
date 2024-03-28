# Load necessary libraries
library(data.table)
library(jsonlite) 

# Set the current working directory to the script directory
setwd(normalizePath(dirname(rstudioapi::getSourceEditorContext()$path)))
getwd()

# Read the metadata file content and create an empty column for case_id
metadata <- fromJSON(readLines("../DNA_Methylation_Metadata.json"))
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
result_file = "DNA_Methylation_Data.txt"
if (file.exists(result_file)) {
  file.remove(result_file)
}

# List all TXT files in the current directory
file_list <- list.files(pattern = "\\.txt", recursive = TRUE)

# Initialize an empty character vector to store all identifiers
max_rownames <- character()

# Loop through each file to collect all identifiers
for (file in file_list) {
  # Read the current file into a data.table
  current_df <- fread(file, header = FALSE, fill = TRUE, select = c(1,2), col.names = c("CpG_site", "methylation"))
  
  # Add unique identifiers to the set
  max_rownames <- union(max_rownames, as.character(current_df$CpG_site))
}

# Create a dataframe from the stored data
result_df <- data.frame("CpG_site" = max_rownames)
counter <- 1

# Loop through each file
for (file in file_list) {
  # Read the current file into a dataframe without considering the first row as header
  current_df <- fread(file, header = FALSE, fill = TRUE, select = c(1,2), col.names = c("CpG_site", "methylation"))
  current_df <- current_df[, .(methylation = mean(methylation)), by = CpG_site]
  extracted_filename <- gsub("^.*/([^/]+)$", "\\1", file)
  
  # Extract the colname and the values column
  current_colname <- as.character(metadata$case_id[metadata$file_name == extracted_filename])
  current_rownames <- current_df$CpG_site
  current_col <- current_df$methylation
  
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

# Truncate sample names in LumA and LumB to first 36 characters and remove duplicates
load("../../luminal_samples.RData")
LumA_truncated <- unique(substr(LumA, 1, 36))
LumB_truncated <- unique(substr(LumB, 1, 36))

# Function to find and rename samples based on truncated IDs
find_and_rename_samples <- function(df, truncated_ids, prefix) {
  # Initialize a list to store matching column names for each ID
  matched_columns_list <- list()
  
  # Loop through each truncated ID to find and store matching columns
  for (id in truncated_ids) {
    matches <- grep(paste0("^", id), names(df), value = TRUE)
    if (length(matches) > 0) {
      matched_columns_list[[id]] <- matches
    }
  }
  
  # Flatten the list to a vector of all matched column names, maintaining original order
  matched_columns <- unlist(matched_columns_list)
  
  if (length(matched_columns) == 0) {
    message("No matching columns found for ", prefix)
    return(data.frame())
  }
  
  # Subset the original dataframe to include only matched columns
  df_subset <- df[, matched_columns, drop = FALSE]
  
  # Rename columns in the subset dataframe based on their sequence and group
  new_names <- sprintf("%s_%02d", prefix, seq_along(matched_columns))
  colnames(df_subset) <- new_names
  
  return(df_subset)
}

# Apply function to LumA and LumB, and combine the results
LumA_df <- find_and_rename_samples(result_df, LumA_truncated, "LumA")
LumB_df <- find_and_rename_samples(result_df, LumB_truncated, "LumB")

# Combine LumA and LumB data frames while preserving row names from the original data frame
count_matrix <- cbind(result_df$CpG_site, LumA_df, LumB_df)
colnames(count_matrix)[1] <- "CpG_site"
save(count_matrix, file = "../DNA_Methylation_Data.RData")

# Write the result dataframe to a text file without quoting string values
write.table(final_df, file = result_file, sep = "\t", row.names = FALSE, quote = FALSE)
