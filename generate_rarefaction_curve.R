# Script to make rarefaction curve for matched bulk and sc PB and ONT data using filtered matrices

library(Seurat)
library(tidyverse)
library(purrr)
library(readr)

# Set working directory to the pipeline location
setwd("/data/gpfs/projects/punim2251/Aim1_LongBench/longbench-analysis-pipeline")

theme_set(theme_minimal())

# Define paths to filtered matrices
BULK_MATRICES_DIR <- "./bulk_processed_matrices"
SC_MATRICES_DIR <- "./sc_processed_matrices"

## Step 1: Load the filtered matrices -------
# Load filtered single-cell matrices
gene_pseudobulk_ont <- readRDS(file.path(SC_MATRICES_DIR, "filtered_genePseudobulkDfs_ont.rds"))
iso_pseudobulk_ont <- readRDS(file.path(SC_MATRICES_DIR, "filtered_isoPseudobulkDfs_ont.rds"))
gene_pseudobulk_pb <- readRDS(file.path(SC_MATRICES_DIR, "filtered_genePseudobulkDfs_pb.rds"))
iso_pseudobulk_pb <- readRDS(file.path(SC_MATRICES_DIR, "filtered_isoPseudobulkDfs_pb.rds"))

# Load filtered bulk matrices
bulk_gene_ont <- readRDS(file.path(BULK_MATRICES_DIR, "bulk_geneDfs_ont_filtered.rds"))
bulk_iso_ont <- readRDS(file.path(BULK_MATRICES_DIR, "bulk_isoDfs_ont_filtered.rds"))
bulk_gene_pb <- readRDS(file.path(BULK_MATRICES_DIR, "bulk_geneDfs_pb_filtered.rds"))
bulk_iso_pb <- readRDS(file.path(BULK_MATRICES_DIR, "bulk_isoDfs_pb_filtered.rds"))

# Load sample sheets for metadata
bulk_sample_sheet <- read.csv(file.path(BULK_MATRICES_DIR, "bulk_sample_sheet.csv"))
target_medians <- unique(bulk_sample_sheet$target_median)

# Function to generate feature number summary from RDS data frame list
generateFeatureNumSummary <- function(df_list, tech_name, data_type, source_type) {
  feature_nums <- data.frame()
  
  # Extract median targets from names if they exist in the format
  for (df_name in names(df_list)) {
    current_df <- df_list[[df_name]]
    
    # Try to extract target median from the df_name if available
    target_median <- NA
    if (grepl("_\\d+M$", df_name)) {
      target_median <- as.numeric(gsub(".*_(\\d+)M$", "\\1", df_name))
    }
    
    # Count features
    if (is.data.frame(current_df)) {
      feature_count <- nrow(current_df)
      
      # Create a standardized sample ID format
      sample_id <- if(source_type == "SC") {
        paste0(tech_name, " Single-Cell ", data_type)
      } else {
        paste0(tech_name, " ", source_type, " ", data_type)
      }
      
      new_row <- data.frame(
        targetMedian = target_median,
        featureNum = feature_count,
        sampleID = sample_id,
        type = paste0(data_type, " Discovery")
      )
      
      feature_nums <- rbind(feature_nums, new_row)
    }
  }
  
  return(feature_nums)
}

# Function to process bulk data frames which might have a different structure
processBulkDfs <- function(bulk_dfs, target_medians, tech_name, data_type) {
  feature_nums <- data.frame()
  
  if (length(bulk_dfs) == length(target_medians)) {
    for (i in seq_along(bulk_dfs)) {
      current_df <- bulk_dfs[[i]]
      target_median <- target_medians[i]
      
      # Count features
      if (is.data.frame(current_df)) {
        feature_count <- nrow(current_df)
        
        # Use a standardized sample ID format
        sample_id <- paste0(tech_name, " Bulk ", data_type)
        
        new_row <- data.frame(
          targetMedian = target_median,
          featureNum = feature_count,
          sampleID = sample_id,
          type = paste0(data_type, " Discovery")
        )
        
        feature_nums <- rbind(feature_nums, new_row)
      }
    }
  }
  
  return(feature_nums)
}

## Step 2: Prepare Rarefaction Curve Dataframe Inputs -------
# Process single-cell data
gene_sc_featureNums_ont <- generateFeatureNumSummary(gene_pseudobulk_ont, "ONT", "Genes", "SC")
tx_sc_featureNums_ont <- generateFeatureNumSummary(iso_pseudobulk_ont, "ONT", "Isoforms", "SC")
gene_sc_featureNums_pb <- generateFeatureNumSummary(gene_pseudobulk_pb, "PacBio", "Genes", "SC")
tx_sc_featureNums_pb <- generateFeatureNumSummary(iso_pseudobulk_pb, "PacBio", "Isoforms", "SC")

# Process bulk data
gene_bulk_featureNums_ont <- processBulkDfs(bulk_gene_ont, target_medians, "ONT", "Genes")
tx_bulk_featureNums_ont <- processBulkDfs(bulk_iso_ont, target_medians, "ONT", "Isoforms")
gene_bulk_featureNums_pb <- processBulkDfs(bulk_gene_pb, target_medians, "PacBio", "Genes")
tx_bulk_featureNums_pb <- processBulkDfs(bulk_iso_pb, target_medians, "PacBio", "Isoforms")

## Step 3: Plotting function ------
plotRarefactionCurve <- function(df_list_input,
                                 plt_title = "Rarefaction Curve",
                                 export_plt = FALSE,
                                 file_prefix = "unnamed_plt",
                                 txt_size = 14,
                                 file_width = 10,
                                 file_height = 6,
                                 removeLegend = FALSE,
                                 output_dir = "rarefaction_plots") {
  # Combine all data frames into one
  rarefactionCurveData <- do.call(rbind, df_list_input)
  
  # Print unique sample IDs for debugging
  print("Unique sample IDs in plot data:")
  print(unique(rarefactionCurveData$sampleID))
  
  # Define a color palette for different sample types - ensure names match exactly with the data
  color_palette <- c(
    "ONT Single-Cell Genes" = "#619CFF",
    "ONT Single-Cell Isoforms" = "#00BA38",
    "PacBio Single-Cell Genes" = "#B79F00",
    "PacBio Single-Cell Isoforms" = "#F8766D",
    "ONT Bulk Genes" = "#00BFC4",
    "ONT Bulk Isoforms" = "#FF61CC",
    "PacBio Bulk Genes" = "#E76BF3",
    "PacBio Bulk Isoforms" = "maroon"
  )
  
  # Create the rarefaction plot
  rarefaction_plt <- ggplot(rarefactionCurveData, aes(x = targetMedian, y = featureNum, color = sampleID)) +
    geom_line() +
    geom_point() +
    labs(x = "Target Median Reads",
         y = "# Unique Features Detected") +
    ggtitle(plt_title) +
    theme(text = element_text(size = txt_size),
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text(angle = 45, hjust = 1),
          legend.title = element_blank()) +
    scale_y_continuous(n.breaks = 10, labels = scales::comma)
  
  # Make sure all unique sample IDs have a color
  unique_samples <- unique(rarefactionCurveData$sampleID)
  if (any(!unique_samples %in% names(color_palette))) {
    missing_samples <- unique_samples[!unique_samples %in% names(color_palette)]
    warning("Missing color definitions for: ", paste(missing_samples, collapse=", "))
    
    # Generate additional colors for missing samples
    additional_colors <- scales::hue_pal()(length(missing_samples))
    names(additional_colors) <- missing_samples
    color_palette <- c(color_palette, additional_colors)
  }
  
  # Apply color scale with the palette - using the subset of colors actually in the data
  needed_colors <- color_palette[unique_samples]
  rarefaction_plt <- rarefaction_plt + 
    scale_color_manual(values = needed_colors)
  
  # Remove legend if specified
  if (removeLegend) {
    rarefaction_plt <- rarefaction_plt + theme(legend.position = "none")
  }
  
  # Save plot if requested
  if (export_plt) {
    # Create directory if it doesn't exist
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    filename <- file.path(output_dir, paste0(file_prefix, ".png"))
    ggsave(filename, plot = rarefaction_plt, 
           width = file_width, height = file_height, 
           dpi = 300)
  }
  
  return(rarefaction_plt)
}

# Create directories for plots
if (!dir.exists("rarefaction_plots")) {
  dir.create("rarefaction_plots")
}
if (!dir.exists("threshold_analysis_plots")) {
  dir.create("threshold_analysis_plots")
}

## Step 4: Create and display rarefaction plots ------
# Combine all feature number data frames
all_feature_nums <- list(
  gene_sc_featureNums_ont,
  tx_sc_featureNums_ont,
  gene_sc_featureNums_pb,
  tx_sc_featureNums_pb,
  gene_bulk_featureNums_ont,
  tx_bulk_featureNums_ont,
  gene_bulk_featureNums_pb,
  tx_bulk_featureNums_pb
)

# Create combined rarefaction plot
rarefaction_plot <- plotRarefactionCurve(
  df_list_input = all_feature_nums,
  plt_title = "Rarefaction Curve: Features Detected in Filtered Data",
  export_plt = TRUE,
  file_prefix = "filtered_rarefaction_curve",
  txt_size = 16,
  file_height = 8,
  file_width = 12,
  output_dir = "rarefaction_plots"
)

# Create separate plots for genes and isoforms
gene_feature_nums <- list(
  gene_sc_featureNums_ont,
  gene_sc_featureNums_pb,
  gene_bulk_featureNums_ont,
  gene_bulk_featureNums_pb
)

gene_rarefaction_plot <- plotRarefactionCurve(
  df_list_input = gene_feature_nums,
  plt_title = "Rarefaction Curve: Genes Detected in Filtered Data",
  export_plt = TRUE,
  file_prefix = "filtered_gene_rarefaction_curve",
  txt_size = 16,
  output_dir = "rarefaction_plots"
)

isoform_feature_nums <- list(
  tx_sc_featureNums_ont,
  tx_sc_featureNums_pb,
  tx_bulk_featureNums_ont,
  tx_bulk_featureNums_pb
)

isoform_rarefaction_plot <- plotRarefactionCurve(
  df_list_input = isoform_feature_nums,
  plt_title = "Rarefaction Curve: Isoforms Detected in Filtered Data",
  export_plt = TRUE,
  file_prefix = "filtered_isoform_rarefaction_curve",
  txt_size = 16,
  output_dir = "rarefaction_plots"
)

# Create technology-specific plots
# ONT only (single-cell vs bulk)
ont_feature_nums <- list(
  gene_sc_featureNums_ont,
  tx_sc_featureNums_ont,
  gene_bulk_featureNums_ont,
  tx_bulk_featureNums_ont
)

ont_rarefaction_plot <- plotRarefactionCurve(
  df_list_input = ont_feature_nums,
  plt_title = "ONT Rarefaction Curve: Bulk vs Single-Cell",
  export_plt = TRUE,
  file_prefix = "ont_rarefaction_curve",
  txt_size = 16,
  file_height = 8,
  file_width = 10,
  output_dir = "rarefaction_plots"
)

# PacBio only (single-cell vs bulk)
pb_feature_nums <- list(
  gene_sc_featureNums_pb,
  tx_sc_featureNums_pb,
  gene_bulk_featureNums_pb,
  tx_bulk_featureNums_pb
)

pb_rarefaction_plot <- plotRarefactionCurve(
  df_list_input = pb_feature_nums,
  plt_title = "PacBio Rarefaction Curve: Bulk vs Single-Cell",
  export_plt = TRUE,
  file_prefix = "pb_rarefaction_curve",
  txt_size = 16,
  file_height = 8,
  file_width = 10,
  output_dir = "rarefaction_plots"
)

# Create data source-specific plots
# Only single-cell data (both technologies)
sc_feature_nums <- list(
  gene_sc_featureNums_ont,
  tx_sc_featureNums_ont,
  gene_sc_featureNums_pb,
  tx_sc_featureNums_pb
)

sc_rarefaction_plot <- plotRarefactionCurve(
  df_list_input = sc_feature_nums,
  plt_title = "Single-Cell Rarefaction Curve: ONT vs PacBio",
  export_plt = TRUE,
  file_prefix = "sc_rarefaction_curve",
  txt_size = 16,
  file_height = 8,
  file_width = 10,
  output_dir = "rarefaction_plots"
)

# Only bulk data (both technologies)
bulk_feature_nums <- list(
  gene_bulk_featureNums_ont,
  tx_bulk_featureNums_ont,
  gene_bulk_featureNums_pb,
  tx_bulk_featureNums_pb
)

bulk_rarefaction_plot <- plotRarefactionCurve(
  df_list_input = bulk_feature_nums,
  plt_title = "Bulk Rarefaction Curve: ONT vs PacBio",
  export_plt = TRUE,
  file_prefix = "bulk_rarefaction_curve",
  txt_size = 16,
  file_height = 8,
  file_width = 10,
  output_dir = "rarefaction_plots"
)

# Display the plots
print(rarefaction_plot)
print(gene_rarefaction_plot)
print(isoform_rarefaction_plot)
print(ont_rarefaction_plot)
print(pb_rarefaction_plot)
print(sc_rarefaction_plot)
print(bulk_rarefaction_plot)

# Save the feature number data for future use
saveRDS(all_feature_nums, "rarefaction_plots/rarefaction_feature_numbers.rds")

### Step 5: Investigating why there are so many more isoforms in pb versus ont:
# Load unfiltered data if not already loaded
iso_ont_unfiltered <- readRDS(file.path(BULK_MATRICES_DIR, "bulk_isoDfs_ont_unfiltered.rds"))
iso_pb_unfiltered <- readRDS(file.path(BULK_MATRICES_DIR, "bulk_isoDfs_pb_unfiltered.rds"))

# Create a function to count rows with counts greater than or equal to threshold
count_rows_above_threshold <- function(df, threshold) {
  df %>% 
    filter(counts >= threshold) %>% 
    nrow()
}

# Define thresholds to test
thresholds <- 1:10

# Initialize data frames to store results
plot_data_all <- tibble()
ratio_data_all <- tibble()

# Process all ONT and PacBio dataframes
for (i in seq_along(iso_ont_unfiltered)) {
  ont_df_name <- names(iso_ont_unfiltered)[i]
  pb_df_name <- names(iso_pb_unfiltered)[i]
  
  # Extract target median from the name if available
  target_median <- NA
  if (grepl("_\\d+M$", ont_df_name)) {
    target_median <- as.numeric(gsub(".*_(\\d+)M$", "\\1", ont_df_name))
  }
  
  # Count isoforms above each threshold
  ont_counts <- map_dbl(thresholds, ~count_rows_above_threshold(
    iso_ont_unfiltered[[i]], .x))
  
  pb_counts <- map_dbl(thresholds, ~count_rows_above_threshold(
    iso_pb_unfiltered[[i]], .x))
  
  # Create dataframe for plotting counts
  temp_plot_data <- tibble(
    threshold = rep(thresholds, 2),
    isoforms = c(ont_counts, pb_counts),
    dataset = rep(c("ONT", "PacBio"), each = length(thresholds)),
    target_median = rep(target_median, 2 * length(thresholds)),
    sample = rep(c(ont_df_name, pb_df_name), each = length(thresholds))
  )
  
  # Create dataframe for plotting ratios
  temp_ratio_data <- tibble(
    threshold = thresholds,
    ont = ont_counts,
    pb = pb_counts,
    ratio = pb_counts / ont_counts,
    target_median = target_median,
    sample_pair = paste0(ont_df_name, " vs ", pb_df_name)
  )
  
  # Append to combined dataframes
  plot_data_all <- bind_rows(plot_data_all, temp_plot_data)
  ratio_data_all <- bind_rows(ratio_data_all, temp_ratio_data)
}

# Calculate average ratio
avg_ratio <- mean(ratio_data_all$ratio, na.rm = TRUE)

# Create the combined plot with all dataframes
count_plot <- ggplot(plot_data_all, aes(x = threshold, y = isoforms, color = dataset, 
                               group = interaction(dataset, target_median),
                               shape = factor(target_median))) +
  geom_line(size = 1) +
  geom_point(size = 2, color = "black") +
  scale_x_continuous(breaks = thresholds) +
  scale_y_continuous(labels = scales::comma, n.breaks = 10) +
  labs(
    title = "Number of Isoforms vs Count Threshold Across All Samples",
    subtitle = paste0("Comparing ONT and PacBio datasets (Avg PacBio:ONT ratio: ", 
                     round(avg_ratio, 2), "x)"),
    x = "Minimum Count Threshold (≥)",
    y = "Number of Isoforms",
    color = "Dataset",
    shape = "Matched Target Median Reads/Cell"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

# Create a ratio plot with all dataframe pairs
ratio_plot <- ggplot(ratio_data_all, aes(x = threshold, y = ratio, 
                                group = target_median, 
                                color = factor(target_median))) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  geom_hline(yintercept = avg_ratio, linetype = "dashed", color = "darkgrey") +
  annotate("text", x = max(thresholds)*0.8, y = avg_ratio*1.05, 
           label = paste0("Avg ratio: ", round(avg_ratio, 2)), color = "darkgrey") +
  scale_x_continuous(breaks = thresholds) +
  scale_y_continuous(n.breaks = 10) +
  labs(
    title = "Magnitude Difference (PacBio/ONT) vs Count Threshold Across All Samples",
    x = "Minimum Count Threshold (≥)",
    y = "Ratio (PacBio/ONT)",
    color = "Matched Target Median Reads/Cell"
  ) +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank()
  )

# Display and save the plots
print(count_plot)
print(ratio_plot)

# Save the threshold analysis plots in their own directory
ggsave("threshold_analysis_plots/isoform_counts_vs_threshold_bulk_data.png", count_plot, width = 12, height = 8)
ggsave("threshold_analysis_plots/isoform_ratio_vs_threshold_bulk_data.png", ratio_plot, width = 12, height = 8)