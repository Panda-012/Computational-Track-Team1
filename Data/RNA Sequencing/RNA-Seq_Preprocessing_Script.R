# Load necessary libraries
library(readr)
library(limma)
library(Rtsne)
library(DESeq2)
library(edgeR)
library(pscl)
library(glmmTMB)
library(data.table)
library(reshape2) # For melt function
library(dplyr) # For data manipulation
library(MCMCglmm)  # For fitting hurdle models
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

# --- Quality control and filtering ---
rm (matching_record,i)

# Check for duplicate gene names
sum(duplicated(rownames(count_matrix)))

# Identify missing values
sum(is.na(count_matrix))
sum(count_matrix == 0)

# Create dge, estimate dispersions, and normalize the data using TMM
dge <- DGEList(counts = count_matrix)
dge <- estimateDisp(dge)
dge <- calcNormFactors(dge, method = "TMM")  # Perform TMM normalization within DGEList

# Subtyping of samples via external R script
Environment_Data <- ls()
Sisi <- cpm(dge, normalized.lib.size = TRUE)
source("Breast_Cancer_Subtyping.R")
table(PAM50Preds$subtype)
rm(list=setdiff(ls(), c("LumA", "LumB", "Basal", "Her2", "Normal", Environment_Data)))

# Create dge, estimate dispersions, and normalize the data using TMM
dge <- DGEList(counts = count_matrix[, c(LumA, LumB, Normal), drop = FALSE])
dge <- estimateDisp(dge)
dge <- calcNormFactors(dge, method = "TMM")  # Perform TMM normalization within DGEList

# Attach sample information dataframe as an additional component to the dge list
dge$samples <- cbind(dge$samples, data.frame(sample_name = colnames(dge$counts),
                          subtype = rep(NA, length(colnames(dge$counts)))))

dge$samples$subtype <- ifelse(colnames(dge$counts) %in% LumA, "Luminal_A",
                       ifelse(colnames(dge$counts) %in% LumB, "Luminal_B",
                       ifelse(colnames(dge$counts) %in% Normal, "Normal", NA)))

###############################################################

# Fit the data to the Zero-Inflated Negative Binomial (ZINB) Model requirements
count_matrix_dt <- as.data.table(dge$counts)
count_matrix_dt[, gene := rownames(dge$counts)]

set.seed(123) # Setting seed for reproducibility
count_matrix_dt <- count_matrix_dt[sample(nrow(count_matrix_dt), 100), ]

total_counts_raw <- colSums(dge$counts)
names(total_counts_raw) <- colnames(dge$counts)

count_matrix_long <- setDT(melt(count_matrix_dt, id.vars = "gene", variable.name = "sample", value.name = "gene_expression"))
count_matrix_long[, total_counts := total_counts_raw[sample], by = sample]

subtype_info_dt <- setDT(as.data.frame(cbind(sample = dge$samples$sample_name, subtype = dge$samples$subtype)))
final_df <- merge(count_matrix_long, subtype_info_dt, by = "sample", all.x = TRUE)



# Fit the Zero-Inflated Negative Binomial (ZINB) Model to the data
zinb_model <- glmmTMB(gene_expression ~ subtype + offset(log(total_counts)),
                      zi = ~ subtype,
                      family = nbinom2,  # This specifies a negative binomial distribution
                      data = final_df)

final_df$predicted_counts <- predict(zinb_model, type = "response")
rownames(count_matrix_dt) <- count_matrix_dt$gene
count_matrix_dt <- as.matrix(count_matrix_dt[ , -122])


# For all samples in the count matrix, replace the zeros with the predicted counts
for (sample_id in colnames(count_matrix_dt)) {
  predicted_count <- unique(final_df$predicted_counts[final_df$sample == sample_id])
  zero_rows <- which(count_matrix_dt[[sample_id]] == 0)
  count_matrix_dt[zero_rows, (sample_id) := predicted_count]  
}

# Check the missing values
sum(is.na(count_matrix_dt))
sum(count_matrix_dt == 0)

# Create dge, estimate dispersions, and normalize the data using TMM
dge <- DGEList(counts = count_matrix_dt)
dge <- estimateDisp(dge)
dge <- calcNormFactors(dge, method = "TMM")  # Perform TMM normalization within DGEList



# Extract normalized counts from the model object
norm_counts <- hurdle_fit$Sol  
sum(norm_counts == 0)
length(boxplot.stats(norm_counts)$out)



################################################################



# Identify outliers (visually and statistically)
boxplot(norm_counts)
affy::plotDensity(norm_counts)
plot(density(apply(norm_counts, 2, mean, na.rm = TRUE)),main="omics",cex.axis=0.5) #1 rows, 2 columns
plotMDS(norm_counts)



# Filter lowly expressed genes
keep <- rowSums(norm_counts >= 10) >= 5  # Adjust thresholds as needed
norm_counts <- norm_counts[keep, ]

# --- Normalization ---

# Calculate normalization factors (e.g., using edgeR)
norm_factors <- calcNormFactors(norm_counts)
count_matrix_norm <- t(t(norm_counts) / norm_factors)

# Reassess quality after normalization
boxplot(count_matrix_norm)
plotMDS(count_matrix_norm)

# --- Dimensionality reduction ---

# Perform PCA
pca_results <- prcomp(t(count_matrix_norm), center = TRUE, scale. = TRUE)

# Visualize PCA results
plot(pca_results$x)  # Color by metadata or clinical variables later

# Create Elbow Graph
variance_explained <- pca_results$sdev^2 / sum(pca_results$sdev^2) * 100
plot(1:length(variance_explained), variance_explained, type = "b", pch = 20)

# Perform t-SNE
set.seed(123)
tsne_results <- Rtsne(pca_results$x[, 1:5], perplexity = 23)  # Adjust perplexity

# Visualize t-SNE results
plot(tsne_results$Y)  # Color by metadata or clinical variables later
