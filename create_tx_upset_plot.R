library(ggplot2)
library(ComplexUpset)

#######################################################
# Define directories and setup
#######################################################

setwd("/data/gpfs/projects/punim2251/Aim1_LongBench/longbench-analysis-pipeline")

# Input directories with processed matrices
SC_INPUT_DIR <- "sc_processed_matrices_with_metadata"
BULK_INPUT_DIR <- "bulk_processed_matrices_with_metadata"

# Output directory for upset plots
OUTPUT_DIR <- "upset_plots"

# Create output directory if it doesn't exist
if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR)
}

#######################################################
# Load transcript data from different modalities
#######################################################

message("Loading transcript data from all modalities...")

# Load single-cell transcript data
sc_tx_ont <- readRDS(file.path(SC_INPUT_DIR, "filtered_isoPseudobulkDfs_ont_with_metadata.rds"))
sc_tx_pb <- readRDS(file.path(SC_INPUT_DIR, "filtered_isoPseudobulkDfs_pb_with_metadata.rds"))

# Load bulk transcript data
bulk_tx_ont <- readRDS(file.path(BULK_INPUT_DIR, "bulk_isoDfs_ont_filtered_with_metadata.rds"))
bulk_tx_pb <- readRDS(file.path(BULK_INPUT_DIR, "bulk_isoDfs_pb_filtered_with_metadata.rds"))
######################################################
######################################################

sequencing_depths = c("1569", 
                      "3923",
                      "7846",
                      "11769",
                      "15692")

for(read_depth in sequencing_depths){
  message(paste0("Processing read depth: ", read_depth, "M..."))
  
  # Find the corresponding matrix from each dataset
  ont_sc_matrix <- sc_tx_ont[[grep(paste0("_", read_depth, "M"), names(sc_tx_ont))[1]]]
  pb_sc_matrix <- sc_tx_pb[[grep(paste0("_", read_depth, "M"), names(sc_tx_pb))[1]]]
  ont_bulk_matrix <- bulk_tx_ont[[grep(paste0("_", read_depth, "M"), names(bulk_tx_ont))[1]]]
  pb_bulk_matrix <- bulk_tx_pb[[grep(paste0("_", read_depth, "M"), names(bulk_tx_pb))[1]]]
  
  # Extract transcript IDs from each matrix
  ont_sc_tx_ids <- ont_sc_matrix$tx_id
  pb_sc_tx_ids <- pb_sc_matrix$tx_id
  ont_bulk_tx_ids <- ont_bulk_matrix$tx_id
  pb_bulk_tx_ids <- pb_bulk_matrix$tx_id
  
  # Create a data frame with all unique transcript IDs
  all_tx_ids <- unique(c(ont_sc_tx_ids, pb_sc_tx_ids, ont_bulk_tx_ids, pb_bulk_tx_ids))
  
  # Create a more comprehensive presence/absence matrix with additional metadata
  # First, collect all necessary data
  all_data <- data.frame(
    tx_id = all_tx_ids,
    ONT_SC = all_tx_ids %in% ont_sc_tx_ids,
    PB_SC = all_tx_ids %in% pb_sc_tx_ids,
    ONT_Bulk = all_tx_ids %in% ont_bulk_tx_ids,
    PB_Bulk = all_tx_ids %in% pb_bulk_tx_ids,
    tx_len = NA,
    counts = NA,
    stringsAsFactors = FALSE
  )
  
  # Modified approach to calculate average counts across modalities
  for(i in 1:nrow(all_data)) {
    tx <- all_data$tx_id[i]
    counts_values <- c()
    
    # Collect counts from all available modalities
    if(all_data$ONT_SC[i]) {
      idx <- which(ont_sc_matrix$tx_id == tx)
      if(length(idx) > 0) {
        all_data$tx_len[i] <- ont_sc_matrix$tx_len[idx[1]]  # Still use tx_len from first source
        counts_values <- c(counts_values, ont_sc_matrix$counts[idx[1]])
      }
    }
    
    if(all_data$PB_SC[i]) {
      idx <- which(pb_sc_matrix$tx_id == tx)
      if(length(idx) > 0) {
        if(is.na(all_data$tx_len[i])) all_data$tx_len[i] <- pb_sc_matrix$tx_len[idx[1]]
        counts_values <- c(counts_values, pb_sc_matrix$counts[idx[1]])
      }
    }
    
    if(all_data$ONT_Bulk[i]) {
      idx <- which(ont_bulk_matrix$tx_id == tx)
      if(length(idx) > 0) {
        if(is.na(all_data$tx_len[i])) all_data$tx_len[i] <- ont_bulk_matrix$tx_len[idx[1]]
        counts_values <- c(counts_values, ont_bulk_matrix$counts[idx[1]])
      }
    }
    
    if(all_data$PB_Bulk[i]) {
      idx <- which(pb_bulk_matrix$tx_id == tx)
      if(length(idx) > 0) {
        if(is.na(all_data$tx_len[i])) all_data$tx_len[i] <- pb_bulk_matrix$tx_len[idx[1]]
        counts_values <- c(counts_values, pb_bulk_matrix$counts[idx[1]])
      }
    }
    
    # Calculate average counts if values were found
    if(length(counts_values) > 0) {
      all_data$counts[i] <- mean(counts_values)
    }
  }
  
  # Set row names for upset function
  rownames(all_data) <- all_data$tx_id
  
  # Create upset plot with annotations for tx_len and counts
  upset_plot <- upset(
    all_data,
    c("ONT_SC", "PB_SC", "ONT_Bulk", "PB_Bulk"),
    name = 'sequencing_modality',
    annotations = list(
      # Transcript length distribution (log10 scale in bp)
      'Transcript Length (log10 bp)' = list(
        aes = aes(x = intersection, y = log10(tx_len + 1)),  # Log10 transformation of bp values
        geom = list(
          geom_boxplot(na.rm = TRUE, fill = "lightblue", alpha = 0.7)  # Show outliers
        )
      ),
      # Read counts distribution (log10 scale)
      'Read Counts (log10)' = list(
        aes = aes(x = intersection, y = log10(counts + 1)),  # log10 scale
        geom = list(
          geom_violin(na.rm = TRUE, fill = "lightgreen", alpha = 0.7),
          geom_jitter(alpha = 0.2, width = 0.2, na.rm = TRUE, size = 0.5, alpha = 0.4)
        )
      )
    ),
    min_size = 5,  # Only show intersections with at least 5 transcripts
    width_ratio = 0.1,
    sort_sets = FALSE,  # Keep the original order of sets
    stripes = "white"
  )
  
  # Save the plot with higher resolution and dimensions for the annotations
  output_file <- file.path(OUTPUT_DIR, paste0("transcript_overlap_", read_depth, "M.png"))
  png(output_file, width = 1800, height = 1200, res = 150)
  print(upset_plot)
  dev.off()
  
  # Save the data for further analysis
  saveRDS(all_data, file.path(OUTPUT_DIR, paste0("transcript_intersection_data_", read_depth, "M.rds")))
  
  message(paste0("Completed upset plot for ", read_depth, "M sequencing depth"))
}

# Create a summary table of transcript counts per modality and shared across modalities
summary_data <- data.frame(
  Depth = character(),
  ONT_SC_Total = integer(),
  PB_SC_Total = integer(),
  ONT_Bulk_Total = integer(),
  PB_Bulk_Total = integer(),
  All_Modalities_Shared = integer(),
  stringsAsFactors = FALSE
)

# Add additional summary data about transcript properties
transcript_properties <- data.frame(
  Depth = character(),
  Intersection = character(),
  Count = integer(),
  Median_Length = numeric(),
  Median_Counts = numeric(),
  stringsAsFactors = FALSE
)

# Create a nested list to store transcript IDs for each intersection set
tx_intersection_ids <- list()

for(read_depth in sequencing_depths) {
  # Read the saved intersection data
  intersection_data <- readRDS(file.path(OUTPUT_DIR, paste0("transcript_intersection_data_", read_depth, "M.rds")))
  
  # Create a list for this depth
  tx_intersection_ids[[paste0(read_depth, "M")]] <- list()
  
  # Print summary statistics for transcript lengths
  message(paste0("\nSummary statistics for transcript lengths at ", read_depth, "M sequencing depth:"))
  print(summary(intersection_data$tx_len))
  message(paste0("99th percentile: ", round(quantile(intersection_data$tx_len, 0.99, na.rm = TRUE))))
  message(paste0("Number of transcripts: ", nrow(intersection_data)))
  
  # Count transcripts in each modality
  ont_sc_count <- sum(intersection_data$ONT_SC)
  pb_sc_count <- sum(intersection_data$PB_SC)
  ont_bulk_count <- sum(intersection_data$ONT_Bulk)
  pb_bulk_count <- sum(intersection_data$PB_Bulk)
  
  # Count transcripts shared across all modalities
  shared_count <- sum(intersection_data$ONT_SC & intersection_data$PB_SC & 
                        intersection_data$ONT_Bulk & intersection_data$PB_Bulk)
  
  # Add to summary data
  summary_data <- rbind(summary_data, data.frame(
    Depth = paste0(read_depth, "M"),
    ONT_SC_Total = ont_sc_count,
    PB_SC_Total = pb_sc_count,
    ONT_Bulk_Total = ont_bulk_count,
    PB_Bulk_Total = pb_bulk_count, 
    All_Modalities_Shared = shared_count,
    stringsAsFactors = FALSE
  ))
  
  # Store transcript IDs for each modality
  tx_intersection_ids[[paste0(read_depth, "M")]][["ONT_SC_Total"]] <- intersection_data$tx_id[intersection_data$ONT_SC]
  tx_intersection_ids[[paste0(read_depth, "M")]][["PB_SC_Total"]] <- intersection_data$tx_id[intersection_data$PB_SC]
  tx_intersection_ids[[paste0(read_depth, "M")]][["ONT_Bulk_Total"]] <- intersection_data$tx_id[intersection_data$ONT_Bulk]
  tx_intersection_ids[[paste0(read_depth, "M")]][["PB_Bulk_Total"]] <- intersection_data$tx_id[intersection_data$PB_Bulk]
  tx_intersection_ids[[paste0(read_depth, "M")]][["All_Modalities_Shared"]] <- 
    intersection_data$tx_id[intersection_data$ONT_SC & intersection_data$PB_SC & 
                            intersection_data$ONT_Bulk & intersection_data$PB_Bulk]
  
  # Calculate properties for each intersection
  # Define all possible intersections
  intersections <- list(
    "ONT_SC_only" = intersection_data$ONT_SC & !intersection_data$PB_SC & 
      !intersection_data$ONT_Bulk & !intersection_data$PB_Bulk,
    "PB_SC_only" = !intersection_data$ONT_SC & intersection_data$PB_SC & 
      !intersection_data$ONT_Bulk & !intersection_data$PB_Bulk,
    "ONT_Bulk_only" = !intersection_data$ONT_SC & !intersection_data$PB_SC & 
      intersection_data$ONT_Bulk & !intersection_data$PB_Bulk,
    "PB_Bulk_only" = !intersection_data$ONT_SC & !intersection_data$PB_SC & 
      !intersection_data$ONT_Bulk & intersection_data$PB_Bulk,
    "All_modalities" = intersection_data$ONT_SC & intersection_data$PB_SC & 
      intersection_data$ONT_Bulk & intersection_data$PB_Bulk,
    "SC_only" = (intersection_data$ONT_SC | intersection_data$PB_SC) & 
      !intersection_data$ONT_Bulk & !intersection_data$PB_Bulk,
    "Bulk_only" = !intersection_data$ONT_SC & !intersection_data$PB_SC & 
      (intersection_data$ONT_Bulk | intersection_data$PB_Bulk),
    "ONT_only" = (intersection_data$ONT_SC | intersection_data$ONT_Bulk) & 
      !intersection_data$PB_SC & !intersection_data$PB_Bulk,
    "PB_only" = !intersection_data$ONT_SC & !intersection_data$ONT_Bulk & 
      (intersection_data$PB_SC | intersection_data$PB_Bulk),
    "ONT_SC_and_PB_SC" = intersection_data$ONT_SC & intersection_data$PB_SC & 
      !intersection_data$ONT_Bulk & !intersection_data$PB_Bulk,
    "ONT_Bulk_and_PB_Bulk" = !intersection_data$ONT_SC & !intersection_data$PB_SC & 
      intersection_data$ONT_Bulk & intersection_data$PB_Bulk,
    "ONT_SC_and_ONT_Bulk" = intersection_data$ONT_SC & !intersection_data$PB_SC & 
      intersection_data$ONT_Bulk & !intersection_data$PB_Bulk,
    "PB_SC_and_PB_Bulk" = !intersection_data$ONT_SC & intersection_data$PB_SC & 
      !intersection_data$ONT_Bulk & intersection_data$PB_Bulk
  )
  
  for(int_name in names(intersections)) {
    subset_data <- intersection_data[intersections[[int_name]], ]
    if(nrow(subset_data) > 0) {
      transcript_properties <- rbind(transcript_properties, data.frame(
        Depth = paste0(read_depth, "M"),
        Intersection = int_name,
        Count = nrow(subset_data),
        Median_Length = median(subset_data$tx_len, na.rm = TRUE),
        Median_Counts = median(subset_data$counts, na.rm = TRUE),
        stringsAsFactors = FALSE
      ))
      
      # Store transcript IDs for this intersection
      tx_intersection_ids[[paste0(read_depth, "M")]][[int_name]] <- subset_data$tx_id
    }
  }
}

# Save summary tables
write.csv(summary_data, file.path(OUTPUT_DIR, "transcript_counts_summary.csv"), row.names = FALSE)
write.csv(transcript_properties, file.path(OUTPUT_DIR, "transcript_properties_summary.csv"), row.names = FALSE)

# Save the transcript IDs for each intersection set
saveRDS(tx_intersection_ids, file.path(OUTPUT_DIR, "transcript_intersection_ids.rds"))
message("Saved transcript IDs for all intersection sets to transcript_intersection_ids.rds")

# Create summary plots
# 1. Bar plot of transcript counts by modality and depth
summary_long <- tidyr::pivot_longer(
  summary_data,
  cols = c("ONT_SC_Total", "PB_SC_Total", "ONT_Bulk_Total", "PB_Bulk_Total", "All_Modalities_Shared"),
  names_to = "Modality",
  values_to = "Count"
)

# Order the Depth factor from lowest to highest sequencing depth
summary_long$Depth <- factor(summary_long$Depth, 
                           levels = paste0(c("1569", "3923", "7846", "11769", "15692"), "M"),
                           ordered = TRUE)

counts_plot <- ggplot(summary_long, aes(x = Depth, y = Count, fill = Modality)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(
    title = "Transcript Counts by Sequencing Modality and Depth",
    x = "Sequencing Depth",
    y = "Number of Transcripts"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save the summary plot
ggsave(file.path(OUTPUT_DIR, "transcript_counts_summary.png"), 
       counts_plot, width = 10, height = 7, dpi = 150)

message("All upset plots, summary data, and summary plots have been created successfully!")