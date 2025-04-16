#!/usr/bin/env Rscript

# Script to create a sample sheet for bulk RNA-seq data and process .quant files
# to match the structure of single-cell pseudobulk matrices

library(tidyverse)
library(fs)

setwd("/data/gpfs/projects/punim2251/Aim1_LongBench/longbench-analysis-pipeline")

# Define paths to directories containing bulk data from oarfish
ont_bulk_dir <- "/data/gpfs/projects/punim2251/rnaseq-minimap-oarfish/results_ontBulk/oarfish"
pb_bulk_dir <- "/data/gpfs/projects/punim2251/rnaseq-minimap-oarfish/results_pbBulk/oarfish"

# Create output directory
output_dir <- "bulk_processed_matrices"
if (!dir_exists(output_dir)) {
  dir_create(output_dir)
}

# Function to parse file name and extract metadata
parse_quant_file_name <- function(file_path) {
  file_name <- basename(file_path)
  
  # Extract technology from directory path
  technology <- if(grepl("ontBulk", file_path)) "ONT" else "PacBio"
  
  # Extract percentage from file name (e.g., 10per, 25per, etc.)
  percentage <- str_extract(file_name, "\\d+per") %>%
    str_replace("per", "")
  
  # Convert percentage to numeric
  percentage_num <- as.numeric(percentage)
  
  # Calculate target median based on percentage
  # This assumes the same mapping as in single-cell data
  target_median <- case_when(
    percentage_num == 10 ~ "1569",
    percentage_num == 25 ~ "3923",
    percentage_num == 50 ~ "7846",
    percentage_num == 75 ~ "11769",
    percentage_num == 100 ~ "15692",
    TRUE ~ NA_character_
  )
  
  # Determine data type - could be gene or transcript level
  # Assuming all are transcript level for now, as they're .quant files
  data_type <- "transcript"
  
  return(tibble(
    file_path = file_path,
    technology = technology,
    percentage = percentage_num,
    target_median = target_median,
    data_type = data_type
  ))
}

# Function to read .quant file and convert to desired format
process_quant_file <- function(file_path, sample_info) {
  message(sprintf("Processing %s", basename(file_path)))
  
  # Read the .quant file
  quant_data <- read_tsv(file_path, col_types = cols(
    tname = col_character(),
    len = col_integer(),
    num_reads = col_integer()
  ))
  
  # Extract gene/transcript information from tname
  transcript_info <- quant_data %>%
    separate(tname, 
             into = c("transcript_id", "gene_id", "other1", "other2", "transcript_name", "gene_name", "length", "biotype"),
             sep = "\\|", 
             fill = "right")
  
  # Create sample ID based on technology
  tech_prefix <- if(sample_info$technology == "ONT") "ont" else "pb"
  sampleID <- paste0(tech_prefix, "_Bulk")
  
  # Create a dataframe with the same structure as the pseudobulk matrices
  bulk_matrix <- tibble(
    # Use transcript_id as row names later
    transcript_id = transcript_info$transcript_id,
    # Use num_reads as counts
    counts = quant_data$num_reads,
    # Add metadata columns
    sampleID = sampleID,
    targetMedian = sample_info$target_median,
    rarefactionDesc = paste0(sampleID, "_", sample_info$percentage, "per_bulk"),
    # Add additional information that might be useful
    gene_id = transcript_info$gene_id,
    gene_name = transcript_info$gene_name,
    transcript_name = transcript_info$transcript_name,
    biotype = transcript_info$biotype
  )
  
  # Sort by counts descending, as done in pseudobulk matrices
  bulk_matrix <- bulk_matrix %>%
    arrange(desc(counts))
  
  return(bulk_matrix)
}

# Function to collapse transcript-level counts to gene-level counts
collapse_to_gene_level <- function(transcript_matrix) {
  message("Collapsing transcript-level counts to gene-level")
  
  # Check if gene_id is available
  if (!("gene_id" %in% colnames(transcript_matrix))) {
    stop("gene_id column not found in transcript matrix")
  }
  
  # Aggregate counts by gene_id
  gene_matrix <- transcript_matrix %>%
    group_by(gene_id) %>%
    summarize(
      counts = sum(counts),
      # Keep metadata columns (take first occurrence)
      sampleID = first(sampleID),
      targetMedian = first(targetMedian),
      rarefactionDesc = first(rarefactionDesc),
      gene_name = first(gene_name),
      .groups = "drop"
    ) %>%
    # Sort by counts descending
    arrange(desc(counts))
  
  # Count number of genes
  message(sprintf("Created gene-level matrix with %d genes", nrow(gene_matrix)))
  
  return(gene_matrix)
}

# Function to filter matrices based on count thresholds
filter_count_matrices <- function(matrix_list, count_threshold) {
  filtered_list <- list()
  
  for (name in names(matrix_list)) {
    original_df <- matrix_list[[name]]
    original_count <- nrow(original_df)
    
    # Filter by count threshold
    filtered_df <- original_df %>%
      filter(counts >= count_threshold)
    
    filtered_count <- nrow(filtered_df)
    
    message(sprintf("Filtered %s: %d features kept out of %d (%.1f%%)", 
                    name, filtered_count, original_count, 
                    100 * filtered_count / original_count))
    
    filtered_list[[name]] <- filtered_df
  }
  
  return(filtered_list)
}

# Find all .quant files in both directories
quant_files <- c(
  dir_ls(ont_bulk_dir, glob = "*.quant"),
  dir_ls(pb_bulk_dir, glob = "*.quant")
)

message(sprintf("Found %d .quant files", length(quant_files)))

# Create sample sheet by parsing file names
sample_sheet <- map_dfr(quant_files, parse_quant_file_name)

# Write sample sheet to file
write_csv(sample_sheet, file.path(output_dir, "bulk_sample_sheet.csv"))
message("Created bulk sample sheet")

# Process each quant file and store results
ont_bulk_tx_matrices <- list()
pb_bulk_tx_matrices <- list()
ont_bulk_gene_matrices <- list()
pb_bulk_gene_matrices <- list()

# Process ONT files
ont_files <- sample_sheet %>% filter(technology == "ONT")
for (i in 1:nrow(ont_files)) {
  file_info <- ont_files[i, ]
  
  # Process transcript-level data
  tx_result <- process_quant_file(file_info$file_path, file_info)
  tx_list_name <- paste0("All_Samples_ONT_tx_", file_info$target_median, "M")
  ont_bulk_tx_matrices[[tx_list_name]] <- tx_result
  
  # Collapse to gene-level
  gene_result <- collapse_to_gene_level(tx_result)
  gene_list_name <- paste0("All_Samples_ONT_gene_", file_info$target_median, "M")
  ont_bulk_gene_matrices[[gene_list_name]] <- gene_result
}

# Process PacBio files
pb_files <- sample_sheet %>% filter(technology == "PacBio")
for (i in 1:nrow(pb_files)) {
  file_info <- pb_files[i, ]
  
  # Process transcript-level data
  tx_result <- process_quant_file(file_info$file_path, file_info)
  tx_list_name <- paste0("All_Samples_PacBio_tx_", file_info$target_median, "M")
  pb_bulk_tx_matrices[[tx_list_name]] <- tx_result
  
  # Collapse to gene-level
  gene_result <- collapse_to_gene_level(tx_result)
  gene_list_name <- paste0("All_Samples_PacBio_gene_", file_info$target_median, "M")
  pb_bulk_gene_matrices[[gene_list_name]] <- gene_result
}

# Apply filtering to transcript-level and gene-level matrices
message("Filtering transcript matrices to keep transcripts with ≥10 counts")
ont_bulk_tx_matrices_filtered <- filter_count_matrices(ont_bulk_tx_matrices, 10)
pb_bulk_tx_matrices_filtered <- filter_count_matrices(pb_bulk_tx_matrices, 10)

message("Filtering gene matrices to keep genes with ≥5 counts")
ont_bulk_gene_matrices_filtered <- filter_count_matrices(ont_bulk_gene_matrices, 5)
pb_bulk_gene_matrices_filtered <- filter_count_matrices(pb_bulk_gene_matrices, 5)

# Save processed and filtered matrices
saveRDS(ont_bulk_tx_matrices_filtered, file.path(output_dir, "bulk_isoDfs_ont_filtered.rds"))
saveRDS(pb_bulk_tx_matrices_filtered, file.path(output_dir, "bulk_isoDfs_pb_filtered.rds"))
saveRDS(ont_bulk_gene_matrices_filtered, file.path(output_dir, "bulk_geneDfs_ont_filtered.rds"))
saveRDS(pb_bulk_gene_matrices_filtered, file.path(output_dir, "bulk_geneDfs_pb_filtered.rds"))

# Also save unfiltered matrices for reference
saveRDS(ont_bulk_tx_matrices, file.path(output_dir, "bulk_isoDfs_ont_unfiltered.rds"))
saveRDS(pb_bulk_tx_matrices, file.path(output_dir, "bulk_isoDfs_pb_unfiltered.rds"))
saveRDS(ont_bulk_gene_matrices, file.path(output_dir, "bulk_geneDfs_ont_unfiltered.rds"))
saveRDS(pb_bulk_gene_matrices, file.path(output_dir, "bulk_geneDfs_pb_unfiltered.rds"))

# Generate summary statistics for filtered transcript data
tx_summary_stats <- tibble(
  technology = character(),
  data_type = character(),
  target_median = character(),
  file_name = character(),
  total_features = integer(),
  features_with_counts = integer(),
  median_counts = numeric(),
  max_counts = integer(),
  filtered_features = integer(),
  percent_retained = numeric()
)

# Add ONT transcript stats
for (name in names(ont_bulk_tx_matrices_filtered)) {
  df_filtered <- ont_bulk_tx_matrices_filtered[[name]]
  df_original <- ont_bulk_tx_matrices[[name]]
  target_med <- unique(df_filtered$targetMedian)
  new_row <- tibble(
    technology = "ONT",
    data_type = "transcript",
    target_median = target_med,
    file_name = basename(ont_files$file_path[ont_files$target_median == target_med]),
    total_features = nrow(df_original),
    features_with_counts = sum(df_original$counts > 0),
    median_counts = median(df_filtered$counts),
    max_counts = max(df_filtered$counts),
    filtered_features = nrow(df_filtered),
    percent_retained = nrow(df_filtered) / nrow(df_original) * 100
  )
  tx_summary_stats <- bind_rows(tx_summary_stats, new_row)
}

# Add PacBio transcript stats
for (name in names(pb_bulk_tx_matrices_filtered)) {
  df_filtered <- pb_bulk_tx_matrices_filtered[[name]]
  df_original <- pb_bulk_tx_matrices[[name]]
  target_med <- unique(df_filtered$targetMedian)
  new_row <- tibble(
    technology = "PacBio",
    data_type = "transcript",
    target_median = target_med,
    file_name = basename(pb_files$file_path[pb_files$target_median == target_med]),
    total_features = nrow(df_original),
    features_with_counts = sum(df_original$counts > 0),
    median_counts = median(df_filtered$counts),
    max_counts = max(df_filtered$counts),
    filtered_features = nrow(df_filtered),
    percent_retained = nrow(df_filtered) / nrow(df_original) * 100
  )
  tx_summary_stats <- bind_rows(tx_summary_stats, new_row)
}

# Generate summary statistics for filtered gene data
gene_summary_stats <- tibble(
  technology = character(),
  data_type = character(),
  target_median = character(),
  file_name = character(),
  total_features = integer(),
  features_with_counts = integer(),
  median_counts = numeric(),
  max_counts = integer(),
  filtered_features = integer(),
  percent_retained = numeric()
)

# Add ONT gene stats
for (name in names(ont_bulk_gene_matrices_filtered)) {
  df_filtered <- ont_bulk_gene_matrices_filtered[[name]]
  df_original <- ont_bulk_gene_matrices[[name]]
  target_med <- unique(df_filtered$targetMedian)
  new_row <- tibble(
    technology = "ONT",
    data_type = "gene",
    target_median = target_med,
    file_name = basename(ont_files$file_path[ont_files$target_median == target_med]),
    total_features = nrow(df_original),
    features_with_counts = sum(df_original$counts > 0),
    median_counts = median(df_filtered$counts),
    max_counts = max(df_filtered$counts),
    filtered_features = nrow(df_filtered),
    percent_retained = nrow(df_filtered) / nrow(df_original) * 100
  )
  gene_summary_stats <- bind_rows(gene_summary_stats, new_row)
}

# Add PacBio gene stats
for (name in names(pb_bulk_gene_matrices_filtered)) {
  df_filtered <- pb_bulk_gene_matrices_filtered[[name]]
  df_original <- pb_bulk_gene_matrices[[name]]
  target_med <- unique(df_filtered$targetMedian)
  new_row <- tibble(
    technology = "PacBio",
    data_type = "gene",
    target_median = target_med,
    file_name = basename(pb_files$file_path[pb_files$target_median == target_med]),
    total_features = nrow(df_original),
    features_with_counts = sum(df_original$counts > 0),
    median_counts = median(df_filtered$counts),
    max_counts = max(df_filtered$counts),
    filtered_features = nrow(df_filtered),
    percent_retained = nrow(df_filtered) / nrow(df_original) * 100
  )
  gene_summary_stats <- bind_rows(gene_summary_stats, new_row)
}

# Combine and save summary stats
summary_stats <- bind_rows(tx_summary_stats, gene_summary_stats)
write_csv(summary_stats, file.path(output_dir, "bulk_processing_summary_filtered.csv"))

message("Processing complete! Filtered matrices and summary saved to: ", output_dir) 