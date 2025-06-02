# Script to Investigate Transcript BioType Distributions
library(gridExtra)
library(grid)
library(tidyverse)
library(patchwork)
library(purrr)

#######################################################
# Define directories and setup
#######################################################

setwd("/data/gpfs/projects/punim2251/Aim1_LongBench/longbench-analysis-pipeline")

# Input directories with processed matrices
SC_INPUT_DIR <- "sc_processed_matrices_with_metadata"
BULK_INPUT_DIR <- "bulk_processed_matrices_with_metadata"

# Output directory for upset plots
OUTPUT_DIR <- "biotype_distribution_plots"

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

# Load the intersection IDs
tx_intersection_ids <- readRDS("upset_plots/transcript_intersection_ids.rds")

# For example, to filter a count matrix for transcripts in ONT_SC_only at 7846M depth:
#ont_sc_only_ids <- tx_ids[["7846M"]][["ONT_SC_only"]]
#filtered_matrix <- your_matrix[your_matrix$tx_id %in% ont_sc_only_ids, ]

######################################################
######################################################

sequencing_depths = c("1569", 
                      "3923",
                      "7846",
                      "11769",
                      "15692")

######################################################
# Step 1: create cleaner tx biotype col
######################################################

# Function to create a cleaner tx_biotpe col in tx-level dataframe
clean_tx_biotype <- function(df) {
  df$tx_biotype_clean <- df$tx_biotype
  
  # Define replacement rules
  df$tx_biotype_clean[grepl("protein_coding|TR|IG", df$tx_biotype_clean)] <- "Protein Coding"
  df$tx_biotype_clean[df$tx_biotype_clean == "TEC"] <- "TEC Putative Protein Coding"
  df$tx_biotype_clean[grepl("pseudogene", df$tx_biotype_clean)] <- "Pseudogene"
  df$tx_biotype_clean[df$tx_biotype_clean == "nonsense_mediated_decay"] <- "Nonsense Mediated Decay"
  df$tx_biotype_clean[df$tx_biotype_clean == "non_stop_decay"] <- "Non-STOP Decay"
  df$tx_biotype_clean[df$tx_biotype_clean == "retained_intron"] <- "Retained Intron"
  df$tx_biotype_clean[df$tx_biotype_clean %in% c("lncRNA", "processed_transcript")] <- "Long Non-Coding RNA"
  df$tx_biotype_clean[grepl("rRNA|tRNA|misc_|vault_|snoRNA|scRNA|snRNA|scaRNA", df$tx_biotype_clean)] <- "Short Non-Coding RNA"
  df$tx_biotype_clean[df$tx_biotype_clean == "artifact"] <- "Artifact"
  df$tx_biotype_clean[df$tx_biotype_clean == "ribozyme"] <- "Ribozyme"
  
  return(df)
}
# Function to clean biotypes across a list of dataframes
clean_biotypes_in_list <- function(df_list) {
  lapply(df_list, function(df) {
    if (is.data.frame(df) && "tx_biotype" %in% colnames(df)) {
      return(clean_tx_biotype(df))
    }
    return(df)
  })
}

# Create cleaner tx_biotype col to all dfs in lists
sc_tx_ont <- clean_biotypes_in_list(sc_tx_ont)
sc_tx_pb <- clean_biotypes_in_list(sc_tx_pb)
bulk_tx_ont <- clean_biotypes_in_list(bulk_tx_ont)
bulk_tx_pb <- clean_biotypes_in_list(bulk_tx_pb)

######################################################
# Step 2: visualize biotype distribution of any intersection set
######################################################

# Function that takes transcript IDs from an intersection set and retrieves their biotype distribution across all depths
analyze_intersection_biotypes_across_depths <- function(intersection_name,
                                                        ont_sc_matrices = sc_tx_ont,
                                                        pb_sc_matrices = sc_tx_pb,
                                                        ont_bulk_matrices = bulk_tx_ont,
                                                        pb_bulk_matrices = bulk_tx_pb,
                                                        intersection_ids = tx_intersection_ids,
                                                        depths = sequencing_depths) {
  
  # Initialize combined data frame
  combined_data <- data.frame()
  
  # Process each depth
  for(depth_str in paste0(depths, "M")) {
    # Extract transcript IDs from the specified intersection at this depth
    tx_ids <- intersection_ids[[depth_str]][[intersection_name]]
    
    if(length(tx_ids) == 0) {
      message(paste("No transcripts found in the", intersection_name, "intersection at", depth_str))
      next
    }
    
    # Find corresponding matrix for the specified depth
    depth_num <- gsub("M", "", depth_str)
    ont_sc_matrix <- ont_sc_matrices[[grep(paste0("_", depth_num, "M"), names(ont_sc_matrices))[1]]]
    pb_sc_matrix <- pb_sc_matrices[[grep(paste0("_", depth_num, "M"), names(pb_sc_matrices))[1]]]
    ont_bulk_matrix <- ont_bulk_matrices[[grep(paste0("_", depth_num, "M"), names(ont_bulk_matrices))[1]]]
    pb_bulk_matrix <- pb_bulk_matrices[[grep(paste0("_", depth_num, "M"), names(pb_bulk_matrices))[1]]]
    
    # Function to extract biotype data from a matrix if it contains the transcript
    extract_biotype_data <- function(matrix_df, source_name) {
      if(is.null(matrix_df)) return(data.frame())
      
      filtered_df <- matrix_df[matrix_df$tx_id %in% tx_ids, ]
      if(nrow(filtered_df) == 0) return(data.frame())
      
      return(data.frame(
        tx_id = filtered_df$tx_id,
        tx_biotype_clean = filtered_df$tx_biotype_clean,
        counts = filtered_df$counts,
        source = source_name,
        stringsAsFactors = FALSE
      ))
    }
    
    # Gather biotype data from all available sources
    biotype_data <- rbind(
      extract_biotype_data(ont_sc_matrix, "ONT SC"),
      extract_biotype_data(pb_sc_matrix, "PacBio SC"),
      extract_biotype_data(ont_bulk_matrix, "ONT Bulk"),
      extract_biotype_data(pb_bulk_matrix, "PacBio Bulk")
    )
    
    # If a transcript appears in multiple sources, use the one with the highest counts
    if(nrow(biotype_data) > 0) {
      biotype_data <- biotype_data %>%
        group_by(tx_id) %>%
        slice_max(order_by = counts, n = 1) %>%
        ungroup() %>%
        # Add depth information
        mutate(depth = depth_str)
      
      combined_data <- rbind(combined_data, biotype_data)
    }
  }
  
  return(combined_data)
}

# Function to plot biotype distribution across all depths for a single intersection
plot_biotypes_across_depths <- function(intersection_name,
                                        plot_title = NULL,
                                        as_percentage = TRUE,
                                        ont_sc_matrices = sc_tx_ont,
                                        pb_sc_matrices = sc_tx_pb,
                                        ont_bulk_matrices = bulk_tx_ont,
                                        pb_bulk_matrices = bulk_tx_pb,
                                        intersection_ids = tx_intersection_ids,
                                        depths = sequencing_depths) {
  
  # Get biotype data across all depths
  all_depth_data <- analyze_intersection_biotypes_across_depths(
    intersection_name,
    ont_sc_matrices,
    pb_sc_matrices,
    ont_bulk_matrices,
    pb_bulk_matrices,
    intersection_ids,
    depths
  )
  
  if(nrow(all_depth_data) == 0) {
    message("No data available to plot")
    return(NULL)
  }
  
  # Create summary data for plotting
  plot_data <- all_depth_data %>%
    group_by(depth, tx_biotype_clean) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(depth) %>%
    mutate(percentage = (n / sum(n)) * 100,
           total_transcripts = sum(n)) %>%
    ungroup()
  
  # Order depths correctly
  plot_data$depth <- factor(plot_data$depth, 
                            levels = paste0(sort(as.numeric(gsub("M", "", unique(plot_data$depth)))), "M"))
  
  # Define a standardized color palette for transcript biotypes
  biotype_colors <- c(
    "Protein Coding" = "#4DAF4A",         # Green
    "Long Non-Coding RNA" = "#377EB8",    # Blue
    "Pseudogene" = "#E41A1C",             # Red
    "Retained Intron" = "#FF7F00",        # Orange
    "Nonsense Mediated Decay" = "#984EA3", # Purple
    "Non-STOP Decay" = "#A65628",         # Brown
    "Short Non-Coding RNA" = "#F781BF",   # Pink
    "TEC Putative Protein Coding" = "#999999", # Gray
    "Artifact" = "#FFFF33",               # Yellow
    "Ribozyme" = "#A6CEE3"                # Light blue
  )
  
  # Calculate the total number of transcripts across all depths
  total_tx_count <- sum(plot_data$n)
  
  # Default title if not provided
  if(is.null(plot_title)) {
    plot_title <- paste0("Transcript Biotype Distribution Across Sequencing Depths\nIntersection: ", 
                         intersection_name)
  }
  
  # Select y-axis based on preference
  y_var <- if(as_percentage) "percentage" else "n"
  y_lab <- if(as_percentage) "Percentage (%)" else "Number of Transcripts"
  
  # Create the plot with standardized colors and subtitle
  p <- ggplot(plot_data, aes(x = depth, y = get(y_var), fill = tx_biotype_clean)) +
    geom_bar(stat = "identity", position = "stack") +
    labs(
      title = plot_title,
      x = "Sequencing Depth",
      y = y_lab,
      fill = "Transcript Biotype"
    ) +
    theme_minimal(base_size = 20) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    ) +
    # Use the standardized color palette
    scale_fill_manual(values = biotype_colors)
  
  return(p)
}

# Function to create biotype distribution summary across all depths and key intersections
create_biotype_summary_plots <- function(output_dir = OUTPUT_DIR,
                                         as_percentage = TRUE) {
  
  # Set of important intersection groups
  key_intersections <- list(
    tech = c("ONT_only", "PB_only"),
    format = c("SC_only", "Bulk_only"),
    specific = c("All_modalities", "ONT_SC_only", "PB_SC_only", "ONT_Bulk_only", "PB_Bulk_only")
  )
  
  # Important biotypes to focus on
  key_biotypes <- c(
    "Protein Coding", "Long Non-Coding RNA", "Pseudogene", 
    "Retained Intron", "Nonsense Mediated Decay", "Short Non-Coding RNA"
  )
  
  # Create plots for each group
  for(group_name in names(key_intersections)) {
    # Create stacked bar plot faceted by intersection
    stacked_plot <- compare_intersections_across_depths(
      key_intersections[[group_name]],
      plot_title = paste("Biotype Distribution by Depth:", group_name),
      as_percentage = as_percentage,
      plot_type = "stacked"
    )
    
    if(!is.null(stacked_plot)) {
      filename <- file.path(output_dir, 
                            paste0("biotype_all_depths_", group_name, "_stacked.png"))
      ggsave(filename, stacked_plot, width = 14, height = 10, dpi = 300)
    }
    
    # Create grouped bar plot faceted by biotype (limited to key biotypes)
    grouped_biotype_plot <- compare_intersections_across_depths(
      key_intersections[[group_name]],
      plot_title = paste("Comparison by Biotype:", group_name),
      as_percentage = as_percentage,
      plot_type = "grouped_biotype",
      biotype_filter = key_biotypes
    )
    
    if(!is.null(grouped_biotype_plot)) {
      filename <- file.path(output_dir, 
                            paste0("biotype_all_depths_", group_name, "_by_biotype.png"))
      ggsave(filename, grouped_biotype_plot, width = 16, height = 12, dpi = 300)
    }
  }
  
  # Create individual plots for important intersections
  important_single_intersections <- c("All_modalities", "ONT_only", "PB_only", "SC_only", "Bulk_only")
  
  for(intsct in important_single_intersections) {
    single_plot <- plot_biotypes_across_depths(
      intsct,
      plot_title = paste("Biotype Distribution Across Depths:", intsct),
      as_percentage = as_percentage
    )
    
    if(!is.null(single_plot)) {
      filename <- file.path(output_dir, 
                            paste0("biotype_all_depths_", intsct, ".png"))
      ggsave(filename, single_plot, width = 12, height = 8, dpi = 300)
    }
  }
  
  message("All summary plots created and saved")
}

# Function to compare multiple intersection sets across all depths
compare_intersections_across_depths <- function(intersection_names,
                                                plot_title = NULL,
                                                as_percentage = TRUE,
                                                ont_sc_matrices = sc_tx_ont,
                                                pb_sc_matrices = sc_tx_pb,
                                                ont_bulk_matrices = bulk_tx_ont,
                                                pb_bulk_matrices = bulk_tx_pb,
                                                intersection_ids = tx_intersection_ids,
                                                depths = sequencing_depths,
                                                plot_type = "stacked",
                                                biotype_filter = NULL) {
  
  # Initialize combined data frame
  combined_data <- data.frame()
  
  # Collect data for each intersection
  for(intsct in intersection_names) {
    # Get biotype data across all depths
    intsct_data <- analyze_intersection_biotypes_across_depths(
      intsct,
      ont_sc_matrices,
      pb_sc_matrices,
      ont_bulk_matrices,
      pb_bulk_matrices,
      intersection_ids,
      depths
    )
    
    if(nrow(intsct_data) > 0) {
      # Add intersection name and combine
      intsct_data$intersection <- intsct
      combined_data <- rbind(combined_data, intsct_data)
    }
  }
  
  if(nrow(combined_data) == 0) {
    message("No data available to plot")
    return(NULL)
  }
  
  # Create summary data for plotting
  plot_data <- combined_data %>%
    group_by(depth, intersection, tx_biotype_clean) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(depth, intersection) %>%
    mutate(percentage = (n / sum(n)) * 100,
           total_transcripts = sum(n)) %>%
    ungroup()
  
  # Filter biotypes if requested
  if(!is.null(biotype_filter)) {
    plot_data <- plot_data %>%
      filter(tx_biotype_clean %in% biotype_filter)
  }
  
  # Order depths correctly
  plot_data$depth <- factor(plot_data$depth, 
                            levels = paste0(sort(as.numeric(gsub("M", "", unique(plot_data$depth)))), "M"))
  
  # Default title if not provided
  if(is.null(plot_title)) {
    plot_title <- "Comparison of Transcript Biotype Distributions Across Sequencing Depths"
  }
  
  # Select y-axis based on preference
  y_var <- if(as_percentage) "percentage" else "n"
  y_lab <- if(as_percentage) "Percentage (%)" else "Number of Transcripts"
  
  # Create summary data for total counts per intersection and depth
  total_counts_data <- plot_data %>%
    group_by(depth, intersection) %>%
    summarise(
      total = sum(get(y_var)),
      total_tx = max(total_transcripts),
      .groups = "drop"
    )
  
  # Choose the plot type
  if(plot_type == "stacked") {
    # Create faceted stacked bar plot by intersection
    p <- ggplot(plot_data, aes(x = depth, y = get(y_var), fill = tx_biotype_clean)) +
      geom_bar(stat = "identity", position = "stack") +
      facet_wrap(~ intersection, scales = "free_y") +
      labs(title = plot_title,
           x = "Sequencing Depth",
           y = y_lab,
           fill = "Transcript Biotype") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(hjust = 0.5)) +
      scale_fill_brewer(palette = "Set3")
    
    # Add total count annotation correctly
    p <- p + geom_text(
      data = total_counts_data,
      aes(x = depth, y = total, label = paste0("n=", format(total_tx, big.mark = ","))),
      vjust = -0.5,
      size = 3
    )
  } else if(plot_type == "grouped_biotype") {
    # Create grouped bar plot by biotype
    p <- ggplot(plot_data, aes(x = depth, y = get(y_var), fill = intersection)) +
      geom_bar(stat = "identity", position = "dodge") +
      facet_wrap(~ tx_biotype_clean, scales = "free_y") +
      labs(title = plot_title,
           x = "Sequencing Depth",
           y = y_lab,
           fill = "Intersection") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(hjust = 0.5)) +
      scale_fill_brewer(palette = "Set1")
  } else {
    # Create grouped bar plot by depth
    p <- ggplot(plot_data, aes(x = tx_biotype_clean, y = get(y_var), fill = intersection)) +
      geom_bar(stat = "identity", position = "dodge") +
      facet_wrap(~ depth, scales = "free_y") +
      labs(title = plot_title,
           x = "Transcript Biotype",
           y = y_lab,
           fill = "Intersection") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(hjust = 0.5)) +
      scale_fill_brewer(palette = "Set1")
  }
  
  return(p)
}

######################################################
# Step 3: Use functions to analyse biotype distributions for intersection sets of interest
######################################################
names(tx_intersection_ids[[1]])
unique(sc_tx_pb$All_Samples_PacBio_sc_tx_15692M$tx_biotype_clean)

# Plot biotype distribution for a single intersection across all depths
plot_biotypes_across_depths("All_Modalities_Shared") %>% 
  ggsave(
    filename = "tx_biotypes_all_modalities_shared.png",
    plot = .,
    path = OUTPUT_DIR,
    width = 14, 
    height = 10, 
    dpi = 600
  )
plot_biotypes_across_depths("PB_SC_only") %>% 
  ggsave(
    filename = "tx_biotypes_pb_sc.png",
    plot = .,
    path = OUTPUT_DIR,
    width = 14, 
    height = 10, 
    dpi = 600
  )
plot_biotypes_across_depths("ONT_SC_only") %>% 
  ggsave(
    filename = "tx_biotypes_ont_sc.png",
    plot = .,
    path = OUTPUT_DIR,
    width = 14, 
    height = 10, 
    dpi = 600
  )
plot_biotypes_across_depths("PB_Bulk_only") %>% 
  ggsave(
    filename = "tx_biotypes_pb_bulk.png",
    plot = .,
    path = OUTPUT_DIR,
    width = 14, 
    height = 10, 
    dpi = 600
  )
plot_biotypes_across_depths("ONT_Bulk_only") %>% 
  ggsave(
    filename = "tx_biotypes_ont_bulk.png",
    plot = .,
    path = OUTPUT_DIR,
    width = 14, 
    height = 10, 
    dpi = 600
  )

# Function to export multiple biotype distribution plots to a PDF
export_biotype_plots_to_pdf <- function(intersection_names,
                                        output_dir = OUTPUT_DIR,
                                        file_name = "biotype_distributions.pdf",
                                        as_percentage = TRUE,
                                        nrow = 1,
                                        ncol = 1) {
  
  # Make sure the output directory exists
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Generate all the plots
  plot_list <- list()
  for (i in seq_along(intersection_names)) {
    plot_list[[i]] <- plot_biotypes_across_depths(
      intersection_names[i],
      as_percentage = as_percentage
    )
  }
  
  # Full file path for the PDF
  pdf_path <- file.path(output_dir, file_name)
  
  # Calculate plots per page and page dimensions
  plots_per_page <- nrow * ncol
  page_width <- 10 * ncol  # Adjust width based on number of columns
  page_height <- 7 * nrow  # Adjust height based on number of rows
  
  # Open PDF device
  pdf(pdf_path, width = page_width, height = page_height)
  
  # If only one plot per page
  if (plots_per_page == 1) {
    # Print each plot on its own page
    for (plot in plot_list) {
      print(plot)  # Use print() instead of grid.draw()
      if (which(plot_list == plot) < length(plot_list)) {
        grid::grid.newpage()
      }
    }
  } else {
    # Use gridExtra to arrange multiple plots per page
    total_plots <- length(plot_list)
    total_pages <- ceiling(total_plots / plots_per_page)
    
    for (page in 1:total_pages) {
      # Calculate which plots go on this page
      start_idx <- (page - 1) * plots_per_page + 1
      end_idx <- min(page * plots_per_page, total_plots)
      
      # Handle case where we don't have enough plots to fill the page
      if (start_idx > total_plots) break
      plots_for_page <- plot_list[start_idx:end_idx]
      
      # For a partial last page
      if (length(plots_for_page) < plots_per_page) {
        # Just use these plots with auto layout
        do.call(gridExtra::grid.arrange, 
                c(plots_for_page, list(ncol = min(ncol, length(plots_for_page)))))
      } else {
        # For full pages use the specified layout
        do.call(gridExtra::grid.arrange, 
                c(plots_for_page, list(ncol = ncol, nrow = nrow)))
      }
      
      # If we're not on the last page, create a new page
      if (page < total_pages) {
        grid::grid.newpage()
      }
    }
  }
  
  # Close the PDF device
  dev.off()
  
  message(paste("PDF saved to:", pdf_path))
}

# Export with 3 rows and 2 columns per page (5 plots per page)
export_biotype_plots_to_pdf(
  intersection_names = c(
    "All_Modalities_Shared",
    "PB_SC_only",
    "ONT_SC_only",
    "PB_Bulk_only",
    "ONT_Bulk_only"
  ),
  file_name = "key_intersection_biotype_distributions_2x2.pdf",
  nrow = 3,
  ncol = 2
)

# Compare specific intersections with a faceted view by biotype
compare_intersections_across_depths(
  c("All_Modalities_Shared", "ONT_SC_only", "PB_SC_only"), 
  plot_type = "grouped_biotype", 
  biotype_filter = c("Protein Coding", "Long Non-Coding RNA", "Retained Intron")
) %>% 
  ggsave(
    filename = "key_biotypes_comparison_ONT_SC_only_and_PB_SC_only.png",
    plot = .,
    path = OUTPUT_DIR,
    width = 14, 
    height = 10, 
    dpi = 300
  )

compare_intersections_across_depths(
  c("All_Modalities_Shared", "ONT_Bulk_only", "PB_Bulk_only"),
  plot_type = "grouped_biotype",
  biotype_filter = c("Protein Coding", "Long Non-Coding RNA", "Retained Intron")
) %>% 
  ggsave(
    filename = "key_biotypes_comparison_ONT_Bulk_only_and_PB_Bulk_only.png",
    plot = .,
    path = OUTPUT_DIR,
    width = 14, 
    height = 10, 
    dpi = 300
  )

compare_intersections_across_depths(
  c("All_Modalities_Shared", "ONT_SC_and_PB_SC", "ONT_Bulk_and_PB_Bulk"),
  plot_type = "grouped_biotype",
  biotype_filter = c("Protein Coding", "Long Non-Coding RNA", "Retained Intron")
) %>% 
  ggsave(
    filename = "key_biotypes_comparison_SC_only_and_Bulk_only.png",
    plot = .,
    path = OUTPUT_DIR,
    width = 14, 
    height = 10, 
    dpi = 300
  )

# Generate all summary plots
create_biotype_summary_plots()
