# Script to convert single-cell RNA-seq matrices to pseudobulk matrices with filtering

library(Seurat)
library(tidyverse)
library(purrr)
library(readr)
library(ggplot2)
library(fs)

setwd("/data/gpfs/projects/punim2251/Aim1_LongBench/longbench-analysis-pipeline")

# Create directory for filtered outputs
output_dir <- "sc_processed_matrices"
if (!dir_exists(output_dir)) {
  dir_create(output_dir)
}

# Create dataframe to store filtering statistics
filtering_stats <- tibble(
  technology = character(),
  data_type = character(),
  target_median = character(),
  percentage = numeric(),
  cells_before = integer(),
  cells_after = integer(), 
  features_before = integer(),
  features_after = integer(),
  novel_transcripts_removed = integer()
)

## Define functions for cell barcode filtering and pseudobulk conversion ---------

# Function to read cell barcodes and create cell line annotation dataframe
read_cell_barcode_files <- function(barcode_dir = "/data/gpfs/projects/punim2251/LongBench_data/cell_line_bc_list/sc") {
  # Get all barcode files
  barcode_files <- list.files(barcode_dir, pattern = "BC_list.txt", full.names = TRUE)
  
  # Initialize empty dataframe to store all barcodes and cell line annotations
  all_barcodes <- data.frame(barcode = character(), cell_line = character(), stringsAsFactors = FALSE)
  
  # Process each barcode file
  for (file in barcode_files) {
    # Extract cell line from filename
    cell_line <- sub("_BC_list.txt", "", basename(file))
    
    # Read barcodes
    barcodes <- readLines(file)
    
    # Add to dataframe
    file_barcodes <- data.frame(barcode = barcodes, cell_line = cell_line, stringsAsFactors = FALSE)
    all_barcodes <- rbind(all_barcodes, file_barcodes)
  }
  
  message(sprintf("Loaded %d cell barcodes across %d cell lines", 
                 nrow(all_barcodes), 
                 length(unique(all_barcodes$cell_line))))
  
  return(all_barcodes)
}

# Function to safely read count matrices with proper error handling
read_count_matrix <- function(file_path, txMatrix = FALSE) {
  tryCatch({
    # Read the CSV file
    count_data <- read.csv(file_path, check.names = FALSE)
    
    # Check if first column can be used as row names
    if (ncol(count_data) > 1 && (is.character(count_data[, 1]) || is.factor(count_data[, 1]))) {
      # Set the first column as row names
      rownames(count_data) <- count_data[, 1]
      count_data <- count_data[, -1, drop = FALSE]
    }
    
    # Special handling for transcript matrices if needed
    if (txMatrix && ncol(count_data) > 0) {
      # Check for additional metadata columns in transcript matrices
      possible_meta_cols <- c("transcript_name", "gene_name", "transcript_id", "gene_id")
      meta_cols <- which(colnames(count_data) %in% possible_meta_cols)
      
      if (length(meta_cols) > 0) {
        # Keep only count columns (not metadata)
        count_cols <- setdiff(1:ncol(count_data), meta_cols)
        if (length(count_cols) > 0) {
          count_data <- count_data[, count_cols, drop = FALSE]
        }
      }
    }
    
    # Ensure rownames exist and are unique
    if (is.null(rownames(count_data)) || any(duplicated(rownames(count_data)))) {
      message("Warning: Matrix has missing or duplicated rownames, generating unique IDs")
      rownames(count_data) <- paste0("feature_", seq_len(nrow(count_data)))
    }
    
    message(sprintf("Successfully read count matrix with %d features and %d cells", 
                   nrow(count_data), 
                   ncol(count_data)))
    
    return(count_data)
  }, error = function(e) {
    message(sprintf("Error reading count matrix from %s: %s", file_path, e$message))
    stop(e)
  })
}

# Function to filter cells and identify high confidence genes/isoforms
filter_high_confidence_features <- function(countMatrix_Path, 
                                           cell_annotations,
                                           txMatrix = FALSE,
                                           gene_min_count = 5,
                                           isoform_min_count = 10,
                                           min_cell_lines = 2,
                                           sample_info = NULL) {
  # Read count matrix with safer function
  counts <- read_count_matrix(countMatrix_Path, txMatrix = txMatrix)
  
  message(sprintf("Count matrix dimensions before processing: %d features x %d cells", 
                 nrow(counts), 
                 ncol(counts)))
  
  # Create initial Seurat object
  seurat_obj <- CreateSeuratObject(counts = counts)
  
  # Get total cell count before filtering
  total_cells_before <- ncol(seurat_obj)
  total_features_before <- nrow(seurat_obj)
  
  # Filter cells to only include those in the barcode list
  valid_cells <- colnames(seurat_obj)[colnames(seurat_obj) %in% cell_annotations$barcode]
  
  if (length(valid_cells) == 0) {
    message("WARNING: No cells matched the whitelist. Check barcode formats.")
    message("Sample cell barcodes from the matrix:")
    print(head(colnames(seurat_obj), 5))
    message("Sample cell barcodes from the whitelist:")
    print(head(cell_annotations$barcode, 5))
    stop("No valid cells found after filtering with whitelist")
  }
  
  message(sprintf("Found %d/%d (%.1f%%) cells in the whitelist", 
                 length(valid_cells), 
                 total_cells_before,
                 length(valid_cells)/total_cells_before * 100))
  
  # Subset to valid cells
  seurat_obj <- subset(seurat_obj, cells = valid_cells)
  
  # Add cell line annotations to the cells
  cell_line_annotations <- cell_annotations %>%
    filter(barcode %in% valid_cells)
  
  # Create a lookup table for faster matching
  cell_line_lookup <- setNames(cell_line_annotations$cell_line, cell_line_annotations$barcode)
  seurat_obj$cell_line <- cell_line_lookup[colnames(seurat_obj)]
  
  # Count number of cells per cell line
  cells_per_line <- table(seurat_obj$cell_line)
  message("Cells per cell line after filtering:")
  print(cells_per_line)
  
  # Get expression matrix
  expr_matrix <- GetAssayData(seurat_obj, slot = "counts")
  
  # Calculate counts per feature per cell line
  cell_lines <- unique(seurat_obj$cell_line)
  feature_counts_by_line <- list()
  
  for(cl in cell_lines) {
    if (is.na(cl)) next
    cl_cells <- WhichCells(seurat_obj, expression = cell_line == cl)
    if (length(cl_cells) > 0) {
      feature_counts_by_line[[cl]] <- rowSums(expr_matrix[, cl_cells, drop = FALSE])
    }
  }
  
  if (length(feature_counts_by_line) == 0) {
    stop("No valid cell lines found after filtering")
  }
  
  # Convert to dataframe for easier filtering
  feature_counts_df <- bind_cols(feature_counts_by_line) %>%
    as.data.frame()
  
  # Ensure rownames are preserved
  if (nrow(feature_counts_df) > 0) {
    rownames(feature_counts_df) <- rownames(expr_matrix)
    feature_counts_df <- rownames_to_column(feature_counts_df, "feature")
  } else {
    stop("No features found after aggregating by cell line")
  }
  
  # Determine minimum count threshold based on data type
  min_count <- if(txMatrix) isoform_min_count else gene_min_count
  
  # Filter for high confidence features
  # For each feature, count how many cell lines have >= min_count counts
  high_conf_features <- feature_counts_df %>%
    pivot_longer(cols = -feature, names_to = "cell_line", values_to = "count") %>%
    filter(count >= min_count) %>%
    group_by(feature) %>%
    summarize(cell_lines_passing = n()) %>%
    filter(cell_lines_passing >= min_cell_lines) %>%
    pull(feature)
  
  # Statistics
  total_features <- nrow(expr_matrix)
  high_conf_count <- length(high_conf_features)
  
  feature_type <- if(txMatrix) "transcripts" else "genes"
  message(sprintf("High confidence %s: %d/%d (%.1f%%)", 
                 feature_type,
                 high_conf_count, 
                 total_features,
                 high_conf_count/total_features * 100))
  
  if (high_conf_count == 0) {
    stop("No high confidence features found with the current thresholds")
  }
  
  # Filter Seurat object to only include high confidence features
  seurat_filtered <- subset(seurat_obj, features = high_conf_features)
  
  # Count novel transcripts before removing them (if transcript data)
  novel_transcripts <- 0
  
  # Remove novel transcripts if dealing with transcripts
  if(txMatrix) {
    features_before <- nrow(seurat_filtered)
    non_novel_features <- rownames(seurat_filtered)[!grepl("Bambu", rownames(seurat_filtered))]
    
    if (length(non_novel_features) > 0) {
      novel_transcripts <- features_before - length(non_novel_features)
      seurat_filtered <- subset(seurat_filtered, features = non_novel_features)
      features_after <- nrow(seurat_filtered)
      message(sprintf("Removed %d novel transcripts", novel_transcripts))
    } else {
      message("No novel transcripts found to remove")
    }
  }
  
  # Add filtering statistics to global dataframe if sample_info provided
  if (!is.null(sample_info)) {
    # Extract values from sample info
    tech <- sample_info$technology
    data_type <- sample_info$data_type
    target_med <- sample_info$target_median
    pct <- sample_info$percentage
    
    # Create a new row for the filtering stats dataframe
    new_stats <- tibble(
      technology = tech,
      data_type = data_type,
      target_median = as.character(target_med),
      percentage = pct,
      cells_before = total_cells_before,
      cells_after = length(valid_cells),
      features_before = total_features_before,
      features_after = nrow(seurat_filtered),
      novel_transcripts_removed = novel_transcripts
    )
    
    # Add to global dataframe using <<- to modify the parent environment
    filtering_stats <<- bind_rows(filtering_stats, new_stats)
  }
  
  return(seurat_filtered)
}

# Function to create pseudobulk matrix from filtered Seurat object
create_pseudobulk <- function(seurat_obj, sample_name, target_median) {
  # Create pseudobulk
  tryCatch({
    # For debugging
    message(sprintf("Creating pseudobulk from Seurat object with %d features and %d cells", 
                   nrow(seurat_obj), ncol(seurat_obj)))
    
    # Check that the RNA assay exists
    if (!"RNA" %in% Assays(seurat_obj)) {
      stop("RNA assay not found in Seurat object")
    }
    
    # Get raw counts directly (safer approach)
    counts_matrix <- GetAssayData(seurat_obj, slot = "counts", assay = "RNA")
    
    # Aggregate manually
    pseudobulk_counts <- rowSums(counts_matrix)
    
    # Convert to dataframe
    pseudobulk <- data.frame(
      counts = pseudobulk_counts,
      sampleID = sample_name,
      targetMedian = target_median,
      rarefactionDesc = paste0(sample_name, "_", target_median, "M_sc")
    )
    
    # Make sure feature names are preserved
    rownames(pseudobulk) <- rownames(counts_matrix)
    
    # Sort by counts descending
    pseudobulk <- pseudobulk[order(pseudobulk$counts, decreasing = TRUE), ]
    
    message("Successfully created pseudobulk matrix")
    return(pseudobulk)
  }, error = function(e) {
    message(sprintf("Error in create_pseudobulk: %s", e$message))
    # Print detailed debugging info
    message("Seurat object structure:")
    print(str(seurat_obj))
    message("Assays in Seurat object:")
    print(Assays(seurat_obj))
    stop(e)
  })
}

# Function to process multiple count matrices using the sample sheet
process_count_matrices <- function(sample_sheet_filtered, sampleID, cell_annotations,
                                  gene_min_count = 5, isoform_min_count = 10, min_cell_lines = 2) {
  
  sampleID <- paste0(sampleID, "_sc")
  
  # Get file paths and target medians from the filtered sample sheet
  matrixFilePaths <- sample_sheet_filtered$file_path
  targetMedians <- as.character(sample_sheet_filtered$target_median)
  
  # Determine if we're handling gene or transcript data
  isGene <- unique(sample_sheet_filtered$data_type) == "gene"
  txMatrix <- !isGene
  
  # Process each matrix
  pseudobulk_list <- list()
  
  for (i in seq_along(matrixFilePaths)) {
    file_path <- matrixFilePaths[i]
    target_median <- targetMedians[i]
    
    # Get the complete row for this sample
    sample_info <- sample_sheet_filtered[i, ]
    
    message(sprintf("\nProcessing %s matrix with target median %s", 
                   if(isGene) "gene" else "transcript", 
                   target_median))
    
    # Try to process this file, continue to next file on error
    tryCatch({
      # Filter for high confidence features
      filtered_seurat <- filter_high_confidence_features(
        countMatrix_Path = file_path,
        cell_annotations = cell_annotations,
        txMatrix = txMatrix,
        gene_min_count = gene_min_count,
        isoform_min_count = isoform_min_count,
        min_cell_lines = min_cell_lines,
        sample_info = sample_info  # Pass sample info for tracking
      )
      
      # Create pseudobulk matrix
      pseudobulk <- create_pseudobulk(filtered_seurat, sampleID, target_median)
      
      # Store in list
      feature_type <- if(isGene) "gene" else "tx"
      list_name <- paste0(sampleID, "_", feature_type, "_", target_median, "M")
      pseudobulk_list[[list_name]] <- pseudobulk
      
      message(sprintf("Successfully processed %s", file_path))
    }, error = function(e) {
      message(sprintf("Error processing %s: %s", file_path, e$message))
      # Print stack trace for debugging
      message("Stack trace:")
      print(rlang::last_trace())
      message("Skipping to next file...")
    })
  }
  
  if (length(pseudobulk_list) == 0) {
    warning("No matrices were successfully processed")
  } else {
    message(sprintf("Successfully processed %d/%d matrices", 
                   length(pseudobulk_list), 
                   length(matrixFilePaths)))
  }
  
  return(pseudobulk_list)
}

# Function to generate filtering report
generate_filtering_report <- function(filtering_stats, output_dir) {
  # Make sure filtering_stats is not empty
  if (nrow(filtering_stats) == 0) {
    message("No filtering stats available to generate report")
    return(NULL)
  }
  
  # Calculate percentage of cells and features lost
  filtering_stats <- filtering_stats %>%
    mutate(
      cells_lost_pct = 100 * (cells_before - cells_after) / cells_before,
      features_lost_pct = 100 * (features_before - features_after) / features_before,
      novel_txs_pct = ifelse(data_type == "transcript", 
                            100 * novel_transcripts_removed / features_before, 
                            0)
    )
  
  # Create summary tables
  tech_summary <- filtering_stats %>%
    group_by(technology, data_type) %>%
    summarize(
      avg_cells_lost_pct = mean(cells_lost_pct),
      avg_features_lost_pct = mean(features_lost_pct),
      avg_novel_txs_pct = mean(novel_txs_pct),
      .groups = "drop"
    )
  
  # Save full stats table as CSV
  filtering_stats %>%
    write_csv(file.path(output_dir, "filtering_statistics.csv"))
  
  # Generate summary table as CSV
  tech_summary %>%
    write_csv(file.path(output_dir, "filtering_summary.csv"))
  
  # Create visualization: Cells lost by technology and rarefaction
  p1 <- ggplot(filtering_stats, aes(x = percentage, y = cells_lost_pct, color = technology)) +
    geom_point() +
    geom_line(aes(group = paste(technology, data_type))) +
    facet_wrap(~data_type) +
    labs(
      title = "Percentage of Cells Filtered Out",
      x = "Rarefaction Percentage",
      y = "% Cells Filtered"
    ) +
    theme_minimal()
  
  # Create visualization: Features lost by technology and rarefaction
  p2 <- ggplot(filtering_stats, aes(x = percentage, y = features_lost_pct, color = technology)) +
    geom_point() +
    geom_line(aes(group = paste(technology, data_type))) +
    facet_wrap(~data_type) +
    labs(
      title = "Percentage of Features Filtered Out",
      x = "Rarefaction Percentage",
      y = "% Features Filtered"
    ) +
    theme_minimal()
  
  # Save plots
  ggsave(file.path(output_dir, "cells_filtered_plot.pdf"), p1, width = 8, height = 6)
  ggsave(file.path(output_dir, "features_filtered_plot.pdf"), p2, width = 8, height = 6)
  
  # If transcript data exists, create transcript-specific plot
  if (any(filtering_stats$data_type == "transcript")) {
    tx_data <- filtering_stats %>% 
      filter(data_type == "transcript")
    
    p3 <- ggplot(tx_data, aes(x = percentage, y = novel_txs_pct, color = technology)) +
      geom_point() +
      geom_line(aes(group = technology)) +
      labs(
        title = "Percentage of Novel Transcripts Removed",
        x = "Rarefaction Percentage",
        y = "% Novel Transcripts"
      ) +
      theme_minimal()
    
    ggsave(file.path(output_dir, "novel_transcripts_plot.pdf"), p3, width = 8, height = 6)
  }
  
  message("Filtering report generated successfully")
  return(TRUE)
}

# Step 1: Read in cell barcode files to create whitelist
cell_annotations <- read_cell_barcode_files()

# Step 2: Read in sample sheet
sample_sheet <- read_csv("sample_sheet.csv")

# Step 3: Process ONT data
# Filter sample sheet for ONT gene data
ont_gene_samples <- sample_sheet %>% 
  filter(technology == "ONT", data_type == "gene")

# Filter sample sheet for ONT transcript data
ont_tx_samples <- sample_sheet %>% 
  filter(technology == "ONT", data_type == "transcript")

# Create filtered pseudobulk matrices for ONT data
message("\n--- Processing ONT gene data ---")
gene.pseudobulk.dfs <- process_count_matrices(
  sample_sheet_filtered = ont_gene_samples,
  sampleID = "All_Samples_ONT",
  cell_annotations = cell_annotations,
  gene_min_count = 5,
  min_cell_lines = 2
)

message("\n--- Processing ONT transcript data ---")
tx.pseudobulk.dfs <- process_count_matrices(
  sample_sheet_filtered = ont_tx_samples,
  sampleID = "All_Samples_ONT",
  cell_annotations = cell_annotations,
  isoform_min_count = 10,
  min_cell_lines = 2
)

# Step 4: Process PacBio data
# Filter sample sheet for PacBio gene data
pb_gene_samples <- sample_sheet %>% 
  filter(technology == "PacBio", data_type == "gene")

# Filter sample sheet for PacBio transcript data
pb_tx_samples <- sample_sheet %>% 
  filter(technology == "PacBio", data_type == "transcript")

# Create filtered pseudobulk matrices for PacBio data
message("\n--- Processing PacBio gene data ---")
gene.pseudobulk.dfs.pb <- process_count_matrices(
  sample_sheet_filtered = pb_gene_samples,
  sampleID = "All_Samples_PacBio",
  cell_annotations = cell_annotations,
  gene_min_count = 5,
  min_cell_lines = 2
)

message("\n--- Processing PacBio transcript data ---")
tx.pseudobulk.dfs.pb <- process_count_matrices(
  sample_sheet_filtered = pb_tx_samples,
  sampleID = "All_Samples_PacBio",
  cell_annotations = cell_annotations,
  isoform_min_count = 10,
  min_cell_lines = 2
)

# Step 5: Save filtered pseudobulk matrices to the output directory
# Save ONT results if available
if (length(gene.pseudobulk.dfs) > 0) {
  saveRDS(gene.pseudobulk.dfs, file.path(output_dir, "filtered_genePseudobulkDfs_ont.rds"))
  message("Saved ONT gene pseudobulk matrices")
}
if (length(tx.pseudobulk.dfs) > 0) {
  saveRDS(tx.pseudobulk.dfs, file.path(output_dir, "filtered_isoPseudobulkDfs_ont.rds"))
  message("Saved ONT isoform pseudobulk matrices")
}

# Save PacBio results if available
if (length(gene.pseudobulk.dfs.pb) > 0) {
  saveRDS(gene.pseudobulk.dfs.pb, file.path(output_dir, "filtered_genePseudobulkDfs_pb.rds"))
  message("Saved PacBio gene pseudobulk matrices")
}
if (length(tx.pseudobulk.dfs.pb) > 0) {
  saveRDS(tx.pseudobulk.dfs.pb, file.path(output_dir, "filtered_isoPseudobulkDfs_pb.rds"))
  message("Saved PacBio isoform pseudobulk matrices")
}

# Step 6: Generate filtering report
generate_filtering_report(filtering_stats, output_dir)
message("Analysis complete! Filtered matrices and reports saved to: ", output_dir)

