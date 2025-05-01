# Script to put together a database of all gene/transcript metadata
library(tidyverse)
library(GenomicFeatures)
library(ensembldb)
library(AnnotationHub)

setwd("/data/gpfs/projects/punim2251/Aim1_LongBench/longbench-analysis-pipeline")

#######################################################
# Section to read and format count matrices if needed 
#######################################################

# Define directories for input and output matrices
SC_MATRICES_DIR <- "sc_processed_matrices"
BULK_MATRICES_DIR <- "bulk_processed_matrices"
SC_OUTPUT_DIR <- "sc_processed_matrices_with_metadata"
BULK_OUTPUT_DIR <- "bulk_processed_matrices_with_metadata"

# Create output directories if they don't exist
if (!dir.exists(SC_OUTPUT_DIR)) {
  dir.create(SC_OUTPUT_DIR)
}

if (!dir.exists(BULK_OUTPUT_DIR)) {
  dir.create(BULK_OUTPUT_DIR)
}

# Read single-cell matrices
message("Reading SC count matrices...")

# Function to read and format matrices
format_matrices <- function(file_path, is_sc = FALSE, is_gene = TRUE, column_mapping = NULL) {
  if (!file.exists(file_path)) {
    message(paste("File not found:", file_path))
    return(NULL)
  }
  
  # Read the RDS file
  df_list <- readRDS(file_path)
  
  # Apply column mapping to each dataframe in the list
  formatted_df_list <- lapply(df_list, function(df) {
    # For SC matrices, handle rownames as feature IDs
    if (is_sc) {
      if (is_gene) {
        # Add gene_id column based on rownames for gene matrices
        df$gene_id <- rownames(df)
      } else {
        # Add tx_id column based on rownames for transcript matrices
        df$tx_id <- rownames(df)
      }
    }
    
    # Apply column renaming if mappings are provided
    if (!is.null(column_mapping)) {
      for (old_name in names(column_mapping)) {
        if (old_name %in% colnames(df)) {
          colnames(df)[colnames(df) == old_name] <- column_mapping[[old_name]]
        }
      }
    }
    
    return(df)
  })
  
  return(formatted_df_list)
}

# Define column mappings if needed
sc_gene_column_mapping <- NULL
sc_tx_column_mapping <- NULL
bulk_gene_column_mapping <- NULL
bulk_tx_column_mapping <- c("transcript_id" = "tx_id")  # Example mapping

# Read and format single-cell matrices - specify is_sc=TRUE
sc_gene_ont <- format_matrices(
  file.path(SC_MATRICES_DIR, "filtered_genePseudobulkDfs_ont.rds"),
  is_sc = TRUE,
  is_gene = TRUE,
  column_mapping = sc_gene_column_mapping
)

sc_gene_pb <- format_matrices(
  file.path(SC_MATRICES_DIR, "filtered_genePseudobulkDfs_pb.rds"),
  is_sc = TRUE,
  is_gene = TRUE,
  column_mapping = sc_gene_column_mapping
)

sc_tx_ont <- format_matrices(
  file.path(SC_MATRICES_DIR, "filtered_isoPseudobulkDfs_ont.rds"),
  is_sc = TRUE,
  is_gene = FALSE,
  column_mapping = sc_tx_column_mapping
)

sc_tx_pb <- format_matrices(
  file.path(SC_MATRICES_DIR, "filtered_isoPseudobulkDfs_pb.rds"),
  is_sc = TRUE,
  is_gene = FALSE,
  column_mapping = sc_tx_column_mapping
)

# Read and format bulk matrices - keep is_sc=FALSE (default)
bulk_gene_ont <- format_matrices(
  file.path(BULK_MATRICES_DIR, "bulk_geneDfs_ont_filtered.rds"),
  column_mapping = bulk_gene_column_mapping
)

bulk_gene_pb <- format_matrices(
  file.path(BULK_MATRICES_DIR, "bulk_geneDfs_pb_filtered.rds"),
  column_mapping = bulk_gene_column_mapping
)

bulk_tx_ont <- format_matrices(
  file.path(BULK_MATRICES_DIR, "bulk_isoDfs_ont_filtered.rds"),
  column_mapping = bulk_tx_column_mapping
)

bulk_tx_pb <- format_matrices(
  file.path(BULK_MATRICES_DIR, "bulk_isoDfs_pb_filtered.rds"),
  column_mapping = bulk_tx_column_mapping
)

# Print diagnostic info
count_matrices_info <- function(matrix_list, name) {
  if (is.null(matrix_list)) {
    message(paste(name, "not found or empty"))
    return()
  }
  
  sample_df <- matrix_list[[1]]
  message(paste(name, "contains", length(matrix_list), "matrices"))
  message(paste("First matrix has", nrow(sample_df), "rows and", ncol(sample_df), "columns"))
  message(paste("Column names:", paste(colnames(sample_df), collapse=", ")))
  # For demonstrating rownames
  if (nrow(sample_df) > 0) {
    message(paste("First few rownames:", paste(head(rownames(sample_df), 3), collapse=", ")))
  }
  message("------------------------")
}

# Print diagnostics
message("\nCount matrices diagnostic information:")
count_matrices_info(sc_gene_ont, "SC ONT Gene")
count_matrices_info(sc_gene_pb, "SC PacBio Gene")
count_matrices_info(sc_tx_ont, "SC ONT Transcript")
count_matrices_info(sc_tx_pb, "SC PacBio Transcript")
count_matrices_info(bulk_gene_ont, "Bulk ONT Gene")
count_matrices_info(bulk_gene_pb, "Bulk PacBio Gene")
count_matrices_info(bulk_tx_ont, "Bulk ONT Transcript")
count_matrices_info(bulk_tx_pb, "Bulk PacBio Transcript")

#######################################################
# Section to generate gene/transcript metadata
#######################################################

## define gtf file path
gtfFilePath <- "/data/gpfs/projects/punim2251/LongBench_data/reference/gencode.v44.annotation.gtf"

# Make the transcript database from gtf
txs <- makeTxDbFromGFF(gtfFilePath, format="gtf")

#Use GTF to pull out information for each transcript/gene in the v46 reference------------
tx.gene.info.v46 <- transcriptLengths(txs, with.cds_len=TRUE, with.utr5_len=TRUE, with.utr3_len=TRUE)

# use annotationhub and ensembldb to get other information----------
ah <- AnnotationHub()
human_ens <- query(ah, c("Homo sapiens", "EnsDb"))
human_ens <- human_ens[["AH116291"]]

# Extract gene-level information
genes(human_ens, return.type = "data.frame")
# Extract transcript-level information
transcripts(human_ens, return.type = "data.frame")

# Create a gene-level dataframe 
annotations_geneMetadata <- genes(human_ens, return.type = "data.frame")  %>%
  dplyr::select(gene_id, gene_biotype, canonical_transcript, symbol) 

# Create a transcript-level dataframe 
annotations_txMetadata <- transcripts(human_ens, return.type = "data.frame")  %>%
  dplyr::select(tx_id, tx_biotype, gc_content, tx_is_canonical, tx_external_name) 

## Create final dataframes with metadata------
# Remove dots from gene and transcript ids
tx.gene.info.v46.formated <- tx.gene.info.v46 %>% 
  dplyr::select(-tx_id) %>% 
  dplyr::rename(tx_id = tx_name) %>% 
  separate(tx_id, c("tx_id"), "\\.", extra = "drop") %>% # remove dots after tx and gene ids
  separate(gene_id, c("gene_id"), "\\.", extra = "drop")

merged.tx.metadata.df <- merge(tx.gene.info.v46.formated, annotations_txMetadata, 
                               by = "tx_id") %>% 
  mutate(utr5_len = if_else(cds_len == 0, NA, utr5_len),
         utr3_len = if_else(cds_len == 0, NA, utr3_len),
         cds_len = if_else(cds_len == 0, NA, cds_len),)
  

# Make a new dataframe with average transcript information corresponding to each gene
gene.avg.len.info <- merged.tx.metadata.df %>% 
  group_by(gene_id) %>% 
  summarise(median.tx.len = median(tx_len),
            median.cds.len = median(cds_len, na.rm = T),
            median.utr5.len = median(utr5_len, na.rm = T),
            median.utr3.len = median(utr3_len, na.rm = T),
            median.gc.content = median(gc_content),
            num_txs = n()) 

merged.gene.metadata.df <- merge(gene.avg.len.info, annotations_geneMetadata, 
                                 by = "gene_id")

saveRDS(merged.gene.metadata.df, file = "geneMetadata.rds")
saveRDS(merged.tx.metadata.df, file = "transcriptMetadata.rds")

#######################################################
# New section to add metadata to filtered matrices
#######################################################

message("Adding metadata to filtered matrices...")

# Function to standardize and add metadata to gene count matrices
add_gene_metadata_to_matrix <- function(matrix_df, gene_metadata) {
  # Process gene IDs to match metadata format (removing version)
  if ("gene_id" %in% colnames(matrix_df)) {
    # Bulk format or SC format with added gene_id column
    matrix_df$gene_id_clean <- sub("\\.[0-9]+$", "", matrix_df$gene_id)
  } else {
    # SC format (if gene_id wasn't already added)
    matrix_df$gene_id <- rownames(matrix_df)
    matrix_df$gene_id_clean <- sub("\\.[0-9]+$", "", matrix_df$gene_id)
  }
  
  # Add gene metadata
  enhanced_df <- left_join(matrix_df, 
                          gene_metadata %>% 
                            dplyr::select(-canonical_transcript) %>%
                            dplyr::rename(gene_id_clean = gene_id, 
                                   gene_symbol = symbol),
                          by = "gene_id_clean")
  
  # Add gene_name if not present (for consistency with bulk format)
  if (!("gene_name" %in% colnames(enhanced_df))) {
    enhanced_df$gene_name <- enhanced_df$gene_symbol
  }
  
  # Standardize column order
  standard_cols <- c("gene_id", "gene_id_clean", "gene_name", "gene_symbol", 
                    "gene_biotype", "median.tx.len", "median.cds.len", 
                    "median.utr5.len", "median.utr3.len", "median.gc.content", 
                    "num_txs", "counts", "sampleID", "targetMedian", "rarefactionDesc")
  
  # Select only columns that exist in the dataframe
  existing_cols <- intersect(standard_cols, colnames(enhanced_df))
  enhanced_df <- enhanced_df[, existing_cols]
  
  # Set rownames for consistency
  rownames(enhanced_df) <- enhanced_df$gene_id
  
  return(enhanced_df)
}

# Function to standardize and add metadata to transcript count matrices
add_tx_metadata_to_matrix <- function(matrix_df, tx_metadata) {
  # Process transcript IDs to match metadata format (removing version)
  if ("tx_id" %in% colnames(matrix_df)) {
    # Bulk format or SC format with added tx_id column
    matrix_df$tx_id_clean <- sub("\\.[0-9]+$", "", matrix_df$tx_id)
  } else {
    # SC format (if tx_id wasn't already added)
    matrix_df$tx_id <- rownames(matrix_df)
    matrix_df$tx_id_clean <- sub("\\.[0-9]+$", "", matrix_df$tx_id)
  }
  
  # Prepare metadata for joining
  metadata_prepared <- tx_metadata %>%
    mutate(tx_id_clean = tx_id) %>%
    dplyr::rename(tx_name = tx_external_name) %>%
    dplyr::select(-tx_id)  # Remove original tx_id to avoid duplicates
  
  # Add transcript metadata
  enhanced_df <- left_join(matrix_df, metadata_prepared, by = "tx_id_clean")
  
  # Add tx_name if not present (for consistency)
  if (!("tx_name" %in% colnames(enhanced_df))) {
    enhanced_df$tx_name <- enhanced_df$tx_id
  }
  
  # Standardize column order
  standard_cols <- c("tx_id", "tx_id_clean", "tx_name", "gene_id", 
                    "tx_biotype", "tx_is_canonical", "tx_len", "cds_len", 
                    "utr5_len", "utr3_len", "gc_content", 
                    "counts", "sampleID", "targetMedian", "rarefactionDesc")
  
  # Select only columns that exist in the dataframe
  existing_cols <- intersect(standard_cols, colnames(enhanced_df))
  enhanced_df <- enhanced_df[, existing_cols]
  
  # Set rownames for consistency
  rownames(enhanced_df) <- enhanced_df$tx_id
  
  return(enhanced_df)
}

# Load gene and transcript metadata
gene_metadata <- readRDS("geneMetadata.rds")
tx_metadata <- readRDS("transcriptMetadata.rds")

#######################################################
# Process preformatted matrices
#######################################################

# Process single-cell matrices
# Gene matrices - ONT
if (!is.null(sc_gene_ont)) {
  message("Processing SC gene matrices for ONT...")
  enhanced_gene_df_list_ont <- lapply(sc_gene_ont, function(df) {
    add_gene_metadata_to_matrix(df, gene_metadata)
  })
  saveRDS(enhanced_gene_df_list_ont, 
          file = file.path(SC_OUTPUT_DIR, "filtered_genePseudobulkDfs_ont_with_metadata.rds"))
}

# Gene matrices - PacBio
if (!is.null(sc_gene_pb)) {
  message("Processing SC gene matrices for PacBio...")
  enhanced_gene_df_list_pb <- lapply(sc_gene_pb, function(df) {
    add_gene_metadata_to_matrix(df, gene_metadata)
  })
  saveRDS(enhanced_gene_df_list_pb, 
          file = file.path(SC_OUTPUT_DIR, "filtered_genePseudobulkDfs_pb_with_metadata.rds"))
}

# Transcript matrices - ONT
if (!is.null(sc_tx_ont)) {
  message("Processing SC transcript matrices for ONT...")
  enhanced_tx_df_list_ont <- lapply(sc_tx_ont, function(df) {
    add_tx_metadata_to_matrix(df, tx_metadata)
  })
  saveRDS(enhanced_tx_df_list_ont, 
          file = file.path(SC_OUTPUT_DIR, "filtered_isoPseudobulkDfs_ont_with_metadata.rds"))
}

# Transcript matrices - PacBio
if (!is.null(sc_tx_pb)) {
  message("Processing SC transcript matrices for PacBio...")
  enhanced_tx_df_list_pb <- lapply(sc_tx_pb, function(df) {
    add_tx_metadata_to_matrix(df, tx_metadata)
  })
  saveRDS(enhanced_tx_df_list_pb, 
          file = file.path(SC_OUTPUT_DIR, "filtered_isoPseudobulkDfs_pb_with_metadata.rds"))
}

# Process bulk matrices
# Gene matrices - ONT
if (!is.null(bulk_gene_ont)) {
  message("Processing bulk gene matrices for ONT...")
  enhanced_gene_df_list_ont <- lapply(bulk_gene_ont, function(df) {
    add_gene_metadata_to_matrix(df, gene_metadata)
  })
  saveRDS(enhanced_gene_df_list_ont, 
          file = file.path(BULK_OUTPUT_DIR, "bulk_geneDfs_ont_filtered_with_metadata.rds"))
}

# Gene matrices - PacBio
if (!is.null(bulk_gene_pb)) {
  message("Processing bulk gene matrices for PacBio...")
  enhanced_gene_df_list_pb <- lapply(bulk_gene_pb, function(df) {
    add_gene_metadata_to_matrix(df, gene_metadata)
  })
  saveRDS(enhanced_gene_df_list_pb, 
          file = file.path(BULK_OUTPUT_DIR, "bulk_geneDfs_pb_filtered_with_metadata.rds"))
}

# Transcript matrices - ONT
if (!is.null(bulk_tx_ont)) {
  message("Processing bulk transcript matrices for ONT...")
  enhanced_tx_df_list_ont <- lapply(bulk_tx_ont, function(df) {
    add_tx_metadata_to_matrix(df, tx_metadata)
  })
  saveRDS(enhanced_tx_df_list_ont, 
          file = file.path(BULK_OUTPUT_DIR, "bulk_isoDfs_ont_filtered_with_metadata.rds"))
}

# Transcript matrices - PacBio
if (!is.null(bulk_tx_pb)) {
  message("Processing bulk transcript matrices for PacBio...")
  enhanced_tx_df_list_pb <- lapply(bulk_tx_pb, function(df) {
    add_tx_metadata_to_matrix(df, tx_metadata)
  })
  saveRDS(enhanced_tx_df_list_pb, 
          file = file.path(BULK_OUTPUT_DIR, "bulk_isoDfs_pb_filtered_with_metadata.rds"))
}

message("Finished adding metadata to all matrices. Standardized matrices are saved in the *_with_metadata directories.")
