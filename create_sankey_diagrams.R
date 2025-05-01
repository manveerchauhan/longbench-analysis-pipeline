# Script to create Sankey diagrams showing gene/transcript flow across sequencing depths
# Based on the filtered matrices used in generate_rarefaction_curve.R

library(tidyverse)
library(purrr)
library(readr)
library(networkD3)  # For creating cleaner Sankey diagrams
library(htmlwidgets)  # For saving interactive HTML widgets
library(htmltools)    # For HTML tags
library(viridis)      # For better color palettes
library(scales)       # For nice axis formatting
library(webshot)      # For capturing screenshots of HTML widgets

# Set working directory to the pipeline location
setwd("/data/gpfs/projects/punim2251/Aim1_LongBench/longbench-analysis-pipeline")

# Define custom colors for Sankey diagram
# These will be used consistently throughout all diagrams
COLOR_PRESERVED <- "#4585b4"  # Blue 
COLOR_GAINED <- "#10fd74"     # Teal/Green
COLOR_LOST <- "#fd7410"       # Red
COLOR_DEPTH_NODE <- "#888888" # Grey for depth nodes

# Define paths to filtered matrices
BULK_MATRICES_DIR <- "./bulk_processed_matrices"
SC_MATRICES_DIR <- "./sc_processed_matrices"

# Create output directory if it doesn't exist
if (!dir.exists("sankey_plots")) {
  dir.create("sankey_plots")
}

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

## Step 2: Function to create Sankey diagrams with networkD3 -------

# Function to get features from a dataframe based on source type
get_features <- function(df, feature_source = "column", feature_column = NULL) {
  if (feature_source == "rownames") {
    return(rownames(df))
  } else if (!is.null(feature_column)) {
    return(pull(df, feature_column))
  } else {
    stop("feature_column must be specified when feature_source is 'column'")
  }
}

# Function to create data for networkD3 Sankey diagram
create_sankey_data <- function(df_list, feature_source = "column", feature_column = NULL) {
  # Extract and sort median targets from names
  median_targets <- names(df_list) %>%
    str_extract("_\\d+M$") %>%
    str_replace("_", "") %>%
    str_replace("M$", "") %>%
    as.numeric() %>%
    sort()
  
  # If median targets couldn't be extracted, use numeric indices
  if (all(is.na(median_targets))) {
    median_targets <- seq_along(df_list)
    names(df_list) <- paste0("Depth_", median_targets)
  } else {
    # Ensure df_list is ordered by median target
    ordered_names <- paste0(str_replace(names(df_list)[1], "_\\d+M$", ""), "_", median_targets, "M")
    df_list <- df_list[ordered_names]
  }
  
  # Create depth labels with the new format
  depth_labels <- paste0(comma(median_targets), " Reads/Cell")
  # Keep a simple version for matching later
  simple_depth_labels <- paste0(median_targets, "M")
  
  # Extract features for each depth
  features_by_depth <- list()
  for (i in seq_along(df_list)) {
    features_by_depth[[i]] <- get_features(df_list[[i]], feature_source, feature_column)
  }
  
  # Create nodes data frame with unique IDs for each node
  nodes <- data.frame(
    name = character(),
    depth = integer(),
    type = character(),
    stringsAsFactors = FALSE
  )
  
  # Add depth nodes
  for (i in seq_along(depth_labels)) {
    nodes <- rbind(nodes, data.frame(
      name = depth_labels[i],
      depth = i,
      type = "depth",
      stringsAsFactors = FALSE
    ))
  }
  
  # Create links data frame
  links <- data.frame(
    source = integer(),
    target = integer(),
    value = integer(),
    status = character(),
    stringsAsFactors = FALSE
  )
  
  # Create nodes and links
  for (i in 1:(length(df_list) - 1)) {
    # Get features in current and next depth
    current_features <- features_by_depth[[i]]
    next_features <- features_by_depth[[i+1]]
    
    # Calculate preserved, lost, and gained features
    preserved <- intersect(current_features, next_features)
    lost <- setdiff(current_features, next_features)
    gained <- setdiff(next_features, current_features)
    
    # Get indices for the depth nodes
    source_idx <- which(nodes$name == depth_labels[i]) - 1  # 0-indexed for networkD3
    target_idx <- which(nodes$name == depth_labels[i+1]) - 1  # 0-indexed for networkD3
    
    # Add preserved features flow
    if (length(preserved) > 0) {
      links <- rbind(links, data.frame(
        source = source_idx,
        target = target_idx,
        value = length(preserved),
        status = "Preserved",
        stringsAsFactors = FALSE
      ))
    }
    
    # Add lost features - create a "Lost" node if it doesn't exist
    if (length(lost) > 0) {
      lost_node_name <- paste0("Lost after ", comma(median_targets[i]), " Reads/Cell\n(", comma(length(lost)), ")")
      if (!lost_node_name %in% nodes$name) {
        nodes <- rbind(nodes, data.frame(
          name = lost_node_name,
          depth = i + 0.5,  # Position between depths
          type = "lost",
          stringsAsFactors = FALSE
        ))
      }
      lost_idx <- which(nodes$name == lost_node_name) - 1  # 0-indexed
      
      links <- rbind(links, data.frame(
        source = source_idx,
        target = lost_idx,
        value = length(lost),
        status = "Lost",
        stringsAsFactors = FALSE
      ))
    }
    
    # Add gained features - create a "Gained" node if it doesn't exist
    if (length(gained) > 0) {
      gained_node_name <- paste0("Gained at ", comma(median_targets[i+1]), " Reads/Cell\n(", comma(length(gained)), ")")
      if (!gained_node_name %in% nodes$name) {
        nodes <- rbind(nodes, data.frame(
          name = gained_node_name,
          depth = i + 0.5,  # Position between depths
          type = "gained",
          stringsAsFactors = FALSE
        ))
      }
      gained_idx <- which(nodes$name == gained_node_name) - 1  # 0-indexed
      
      links <- rbind(links, data.frame(
        source = gained_idx,
        target = target_idx,
        value = length(gained),
        status = "Gained",
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # Reindex nodes and links for networkD3 (0-indexed)
  nodes$id <- 0:(nrow(nodes) - 1)
  
  # Update link source and target to match node ids
  links$source <- as.integer(links$source)
  links$target <- as.integer(links$target)
  
  # Add node type to links for coloring
  links$source_type <- nodes$type[links$source + 1]
  links$target_type <- nodes$type[links$target + 1]
  
  # Add label with count for links
  links$label <- paste0(links$status, ": ", comma(links$value))
  
  # Improve node names with counts for depth nodes only (gained/lost already have counts)
  for (i in 1:nrow(nodes)) {
    if (nodes$type[i] == "depth") {
      # For depth nodes, add the count of features
      depth_idx <- which(depth_labels == nodes$name[i])
      feature_count <- length(features_by_depth[[depth_idx]])
      nodes$name[i] <- paste0(nodes$name[i], "\n(", comma(feature_count), ")")
    }
  }
  
  return(list(nodes = nodes, links = links))
}

# Function to create and save a Sankey diagram using networkD3
create_sankey_diagram <- function(df_list, feature_source = "column", feature_column = NULL, 
                                  title = "Gene Flow Across Sequencing Depths",
                                  file_prefix = "unnamed_sankey", 
                                  width = 1000, height = 600) {
  
  # Create Sankey data
  sankey_data <- create_sankey_data(df_list, feature_source, feature_column)
  
  # Use the custom colors defined at the top of the script
  status_colors <- c(
    "Preserved" = COLOR_PRESERVED,
    "Gained" = COLOR_GAINED,
    "Lost" = COLOR_LOST
  )
  
  # Prepare node colors based on type using the same colors
  node_colors <- c(
    "depth" = COLOR_DEPTH_NODE,
    "lost" = COLOR_LOST,
    "gained" = COLOR_GAINED
  )
  
  # Assign colors to nodes
  nodes_color <- node_colors[sankey_data$nodes$type]
  
  # Assign colors to links based on status
  links_color <- status_colors[sankey_data$links$status]
  
  # Create the Sankey diagram
  sankey <- sankeyNetwork(
    Links = sankey_data$links,
    Nodes = sankey_data$nodes,
    Source = "source",
    Target = "target",
    Value = "value",
    NodeID = "name",
    NodeGroup = "type",
    LinkGroup = "status",
    colourScale = JS(paste0('d3.scaleOrdinal().domain(["Preserved", "Gained", "Lost"]).range(["', 
                            paste(status_colors, collapse = '", "'), '"])')),
    nodeWidth = 30,
    nodePadding = 10,
    fontSize = 12,
    fontFamily = "Arial",
    iterations = 64,
    sinksRight = FALSE,  # Allow manual node placement
    width = width,
    height = height
  )
  
  # Add labels to the flows
  sankey <- htmlwidgets::onRender(
    sankey,
    '
    function(el, x) {
      // Add labels to links
      d3.select(el).select("svg")
        .selectAll(".link")
        .append("title")
        .text(function(d) { return d.source.name + " → " + d.target.name + "\\n" + d.value; });
      
      // Add value labels on links
      var links = d3.select(el).select("svg").selectAll(".link");
      
      links.each(function(d) {
        var link = d3.select(this);
        var linkData = link.datum();
        
        // Calculate position for the text
        var x1 = (linkData.source.x + linkData.target.x) / 2;
        var y1 = linkData.y0 + (linkData.y1 - linkData.y0) / 2;
        
        // Append text element
        d3.select(el).select("svg")
          .append("text")
          .attr("x", x1)
          .attr("y", y1)
          .attr("dy", ".35em")
          .attr("text-anchor", "middle")
          .attr("font-size", "10px")
          .attr("fill", "black") // Black for better visibility on all backgrounds
          .attr("pointer-events", "none")
          .text(linkData.value);
      });
    }
    '
  )
  
  # Add title to the plot using htmltools::tags
  sankey <- htmlwidgets::prependContent(
    sankey,
    htmltools::tags$h3(
      title,
      style = "text-align: center; font-family: Arial; margin-top: 0;"
    )
  )
  
  # Save as HTML (interactive)
  html_path <- file.path("sankey_plots", paste0(file_prefix, ".html"))
  saveWidget(sankey, html_path, selfcontained = TRUE)
  
  # Check if webshot package is available, if not, try to install it
  if (!requireNamespace("webshot", quietly = TRUE)) {
    try(install.packages("webshot", repos = "https://cran.rstudio.com/"))
    try(webshot::install_phantomjs())
  }
  
  # Try to save static image for reports, with error handling
  webshot_path <- file.path("sankey_plots", paste0(file_prefix, ".png"))
  tryCatch({
    webshot::webshot(html_path, webshot_path, delay = 1, vwidth = width, vheight = height + 50)
    cat(paste0("Saved ", file_prefix, " as HTML and PNG\n"))
  }, error = function(e) {
    cat(paste0("Saved ", file_prefix, " as HTML. Error generating PNG: ", e$message, "\n"))
  })
  
  return(sankey)
}

## Step 3: Create Sankey diagrams for each dataset -------

# ONT Single-Cell Genes
ont_sc_gene_sankey <- create_sankey_diagram(
  gene_pseudobulk_ont,
  feature_source = "rownames",
  title = "ONT Single-Cell Gene Flow Across Sequencing Depths",
  file_prefix = "ont_sc_gene_sankey"
)

# ONT Single-Cell Isoforms
ont_sc_iso_sankey <- create_sankey_diagram(
  iso_pseudobulk_ont,
  feature_source = "rownames",
  title = "ONT Single-Cell Isoform Flow Across Sequencing Depths",
  file_prefix = "ont_sc_iso_sankey"
)

# PacBio Single-Cell Genes
pb_sc_gene_sankey <- create_sankey_diagram(
  gene_pseudobulk_pb,
  feature_source = "rownames",
  title = "PacBio Single-Cell Gene Flow Across Sequencing Depths",
  file_prefix = "pb_sc_gene_sankey"
)

# PacBio Single-Cell Isoforms
pb_sc_iso_sankey <- create_sankey_diagram(
  iso_pseudobulk_pb,
  feature_source = "rownames",
  title = "PacBio Single-Cell Isoform Flow Across Sequencing Depths",
  file_prefix = "pb_sc_iso_sankey"
)

# ONT Bulk Genes
ont_bulk_gene_sankey <- create_sankey_diagram(
  bulk_gene_ont,
  feature_source = "column",
  feature_column = "gene_id",
  title = "ONT Bulk Gene Flow Across Sequencing Depths",
  file_prefix = "ont_bulk_gene_sankey"
)

# ONT Bulk Isoforms
ont_bulk_iso_sankey <- create_sankey_diagram(
  bulk_iso_ont,
  feature_source = "column",
  feature_column = "transcript_id",
  title = "ONT Bulk Isoform Flow Across Sequencing Depths",
  file_prefix = "ont_bulk_iso_sankey"
)

# PacBio Bulk Genes
pb_bulk_gene_sankey <- create_sankey_diagram(
  bulk_gene_pb,
  feature_source = "column",
  feature_column = "gene_id",
  title = "PacBio Bulk Gene Flow Across Sequencing Depths",
  file_prefix = "pb_bulk_gene_sankey"
)

# PacBio Bulk Isoforms
pb_bulk_iso_sankey <- create_sankey_diagram(
  bulk_iso_pb,
  feature_source = "column",
  feature_column = "transcript_id",
  title = "PacBio Bulk Isoform Flow Across Sequencing Depths",
  file_prefix = "pb_bulk_iso_sankey"
)

## Step 4: Create summary statistics -------
# Create a function to summarize feature retention and discovery across depths
summarize_feature_flows <- function() {
  # Initialize empty dataframe for summary stats
  summary_stats <- tibble(
    dataset = character(),
    feature_type = character(),
    total_unique_features = numeric(),
    features_preserved_all_depths = numeric(),
    percent_preserved = numeric(),
    new_features_highest_depth = numeric()
  )
  
  # Function to get features based on dataset type
  get_features_for_dataset <- function(df, dataset_type) {
    if (grepl("SC", dataset_type)) {
      return(rownames(df))
    } else if (grepl("Genes", dataset_type)) {
      return(pull(df, "gene_id"))
    } else {
      return(pull(df, "transcript_id"))
    }
  }
  
  # Define datasets to process
  datasets <- list(
    list(name = "ONT_SC_Genes", data = gene_pseudobulk_ont, col = NULL),
    list(name = "ONT_SC_Isoforms", data = iso_pseudobulk_ont, col = NULL),
    list(name = "PacBio_SC_Genes", data = gene_pseudobulk_pb, col = NULL),
    list(name = "PacBio_SC_Isoforms", data = iso_pseudobulk_pb, col = NULL),
    list(name = "ONT_Bulk_Genes", data = bulk_gene_ont, col = "gene_id"),
    list(name = "ONT_Bulk_Isoforms", data = bulk_iso_ont, col = "transcript_id"),
    list(name = "PacBio_Bulk_Genes", data = bulk_gene_pb, col = "gene_id"),
    list(name = "PacBio_Bulk_Isoforms", data = bulk_iso_pb, col = "transcript_id")
  )
  
  # Process each dataset
  for (ds in datasets) {
    if (length(ds$data) < 2) next
    
    # Get all features across all depths
    all_features <- map(ds$data, function(df) {
      if (grepl("SC", ds$name)) {
        return(rownames(df))
      } else if (grepl("Genes", ds$name)) {
        return(pull(df, "gene_id"))
      } else {
        return(pull(df, "transcript_id"))
      }
    }) %>%
      unlist() %>%
      unique()
    
    # Get features in the first depth
    first_features <- if (grepl("SC", ds$name)) {
      rownames(ds$data[[1]])
    } else if (grepl("Genes", ds$name)) {
      pull(ds$data[[1]], "gene_id")
    } else {
      pull(ds$data[[1]], "transcript_id")
    }
    
    # Get features in the last depth
    last_features <- if (grepl("SC", ds$name)) {
      rownames(ds$data[[length(ds$data)]])
    } else if (grepl("Genes", ds$name)) {
      pull(ds$data[[length(ds$data)]], "gene_id")
    } else {
      pull(ds$data[[length(ds$data)]], "transcript_id")
    }
    
    # Get features present in all depths
    features_all_depths <- ds$data %>%
      map(function(df) {
        if (grepl("SC", ds$name)) {
          return(rownames(df))
        } else if (grepl("Genes", ds$name)) {
          return(pull(df, "gene_id"))
        } else {
          return(pull(df, "transcript_id"))
        }
      }) %>%
      reduce(intersect)
    
    # Get new features in the highest depth (not in any previous depth)
    prev_all_features <- ds$data[-length(ds$data)] %>%
      map(function(df) {
        if (grepl("SC", ds$name)) {
          return(rownames(df))
        } else if (grepl("Genes", ds$name)) {
          return(pull(df, "gene_id"))
        } else {
          return(pull(df, "transcript_id"))
        }
      }) %>%
      unlist() %>%
      unique()
    
    new_features_highest <- setdiff(last_features, prev_all_features)
    
    # Add to summary stats
    summary_stats <- summary_stats %>%
      add_row(
        dataset = ds$name,
        feature_type = ifelse(grepl("Genes", ds$name), "Genes", "Isoforms"),
        total_unique_features = length(all_features),
        features_preserved_all_depths = length(features_all_depths),
        percent_preserved = round(length(features_all_depths) / length(first_features) * 100, 1),
        new_features_highest_depth = length(new_features_highest)
      )
  }
  
  # Save summary stats
  write_csv(summary_stats, "sankey_plots/feature_flow_summary.csv")
  
  return(summary_stats)
}

# Generate and print summary statistics
flow_summary <- summarize_feature_flows()
print(flow_summary)

cat("Interactive Sankey diagrams created and saved to the 'sankey_plots' directory as HTML files.\n")
cat("Static versions saved as PNG files for reports.\n")
cat("Summary statistics saved to 'sankey_plots/feature_flow_summary.csv'.\n") 