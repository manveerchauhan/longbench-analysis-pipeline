# Script to make rarefaction curve for matched bulk and sc PB and ONT data for all cell lines

library(Seurat)
library(tidyverse)
library(purrr)
library(readr)

setwd("/data/gpfs/projects/punim2251/Aim1_LongBench/ReadRarefaction_wFixedCells/scripts/rarefactionCurveAnalysis")
theme_set(ggmin::theme_min())

## Define functions used in the script---------
# Function to return Pseudobulk count matrix for a given count matrix file path
returnPseudoBulkMatrix <- function(countMatrix_Path, txMatrix = F, 
                                   removeNovelTxs = T, sample = "NA",
                                   targetMedianReads = "NA"){
  if(txMatrix){
    pseudobulkMatrix <- read.csv(countMatrix_Path, row.names = 1)[-1]
  } else{
    pseudobulkMatrix <- read.csv(countMatrix_Path, row.names = 1)
  }
  
  pseudobulkMatrix <- pseudobulkMatrix %>% 
    CreateSeuratObject(counts = .) %>% 
    AggregateExpression(object = .) %>% 
    as.data.frame() %>%
    dplyr::rename(counts = RNA) %>%
    arrange(desc(counts)) %>%
    mutate(sampleID = sample,
           targetMedian = targetMedianReads,
           rarefactionDesc = paste0(sample, "_", targetMedianReads, "M_sc"))
  
  if(removeNovelTxs){
    message("Removing bambu (novel) features")
    pseudobulkMatrix <- pseudobulkMatrix %>%
      filter(!grepl("Bambu", rownames(.)))
  }
  
  return(pseudobulkMatrix)
}
# Function to convert list of rarefied pseudobulk/bulk matrices into a df showing number of features in each matrix
generateFeatureNumSummary <- function(pseudobulk.df.list){
  num.Features.Summary <- data.frame()
  
  for(pseudobulk.matrix in names(pseudobulk.df.list)){
    current.matrix <- pseudobulk.df.list[[pseudobulk.matrix]]
    
    summary <- current.matrix %>%
      group_by(rarefactionDesc, sampleID, targetMedian) %>%
      summarise(featureNum = n(), .groups = 'drop')
    
    num.Features.Summary <- rbind(num.Features.Summary, summary)
  }
  
  num.Features.Summary <- num.Features.Summary %>%
    dplyr::mutate(targetMedian = as.numeric(targetMedian),
                  featureNum = as.numeric(featureNum))
  
  return(num.Features.Summary)
}

# Function to create a list of pseudobulk gene matrices for each rarified gene/tx count matrix file path
createPseudobulkDfList <- function(matrixFilePaths, targetMedians, sampleID, 
                                   isGene = T, filterNovelTxs = T){
  
  sampleID <- paste0(sampleID, "_sc")
  
  
  if(isGene){
    suppressWarnings({
      pseudobulk.dfs.list <- map2(matrixFilePaths, targetMedians, 
                                  ~returnPseudoBulkMatrix(.x, 
                                                          sample = sampleID, 
                                                          targetMedianReads = .y,
                                                          removeNovelTxs = filterNovelTxs))
    })
    names(pseudobulk.dfs.list) <- paste0(sampleID, "_gene_", targetMedians, "M")
  }else{
    suppressWarnings({
      pseudobulk.dfs.list <- map2(matrixFilePaths, targetMedians, 
                                  ~returnPseudoBulkMatrix(.x,
                                                          txMatrix = TRUE,
                                                          sample = sampleID, 
                                                          targetMedianReads = .y,
                                                          removeNovelTxs = filterNovelTxs))
    })
    names(pseudobulk.dfs.list) <- paste0(sampleID, "_tx_", targetMedians, "M")
  }
  
  return(pseudobulk.dfs.list)
}
# Function to read in bulk count data (from Isoquant) into the same format as sc pseudobulk matrices
readInBulkCountData <- function(countMatrix_Path,
                                sample = "NA",
                                targetMedianReads = "NA"){
  bulk_df <- read_tsv(countMatrix_Path) %>% 
    filter(!str_starts(`#feature_id`, "_")) %>%
    column_to_rownames(var = "#feature_id") %>% 
    dplyr::rename(counts = count) %>%
    mutate(sampleID = sample,
           targetMedian = targetMedianReads,
           rarefactionDesc = paste0(sample, "_", targetMedianReads, "M_sc"))
  
  return(bulk_df)
}

# Function that calls readInBulkCountData for a list of file paths and returns a formatted list
createBulkDfList <- function(matrixFilePaths, targetMedians, 
                             sampleID, 
                             isGene = T){
  
  sampleID <- paste0(sampleID, "_bulk")
  
  
  if(isGene){
    suppressWarnings({
      bulk.dfs.list <- map2(matrixFilePaths, targetMedians,
                            ~readInBulkCountData(.x, 
                                                 sample = sampleID,
                                                 targetMedianReads = .y))
    })
    names(bulk.dfs.list) <- paste0(sampleID, "_gene_", targetMedians, "M")
  }else{
    suppressWarnings({
      bulk.dfs.list <- map2(matrixFilePaths, targetMedians,
                            ~readInBulkCountData(.x, 
                                                 sample = sampleID,
                                                 targetMedianReads = .y))
    })
    names(bulk.dfs.list) <- paste0(sampleID, "_tx_", targetMedians, "M")
  }
  
  return(bulk.dfs.list)
}


# Step 1: Read in sc gene/tx count files as pseudobulk matrices-------
allSamples.ONT.sc.geneM.paths <- c(
  "/data/gpfs/projects/punim2251/Aim1_LongBench/ReadRarefaction_wFixedCells/data/LongBench_All/ont_sc/10Percent_FLAMES/gene_count.csv",
  "/data/gpfs/projects/punim2251/Aim1_LongBench/ReadRarefaction_wFixedCells/data/LongBench_All/ont_sc/25Percent_FLAMES/gene_count.csv",
  "/data/gpfs/projects/punim2251/Aim1_LongBench/ReadRarefaction_wFixedCells/data/LongBench_All/ont_sc/50Percent_FLAMES/gene_count.csv",
  "/data/gpfs/projects/punim2251/Aim1_LongBench/ReadRarefaction_wFixedCells/data/LongBench_All/ont_sc/75Percent_FLAMES/gene_count.csv",
  "/data/gpfs/projects/punim2251/LongBench_data/flames_out/ont_sc/gene_count.csv"
)

allSamples.ONT.sc.txM.paths <- c(
  "/data/gpfs/projects/punim2251/Aim1_LongBench/ReadRarefaction_wFixedCells/data/LongBench_All/ont_sc/10Percent_FLAMES/allSamples_ONT_sc_10Percent_iso_count.csv",
  "/data/gpfs/projects/punim2251/Aim1_LongBench/ReadRarefaction_wFixedCells/data/LongBench_All/ont_sc/25Percent_FLAMES/allSamples_ONT_sc_25Percent_iso_count.csv",
  "/data/gpfs/projects/punim2251/Aim1_LongBench/ReadRarefaction_wFixedCells/data/LongBench_All/ont_sc/50Percent_FLAMES/allSamples_ONT_sc_50Percent_iso_count.csv",
  "/data/gpfs/projects/punim2251/Aim1_LongBench/ReadRarefaction_wFixedCells/data/LongBench_All/ont_sc/75Percent_FLAMES/allSamples_ONT_sc_75Percent_iso_count.csv",
  "/data/gpfs/projects/punim2251/LongBench_data/flames_out/ont_sc/allSamples_ONT_sc_100Percent_iso_count.csv"
)
target_medians <- c("1569", "3923", "7846", "11769", "15692")

gene.pseudobulk.dfs <- createPseudobulkDfList(matrixFilePaths = allSamples.ONT.sc.geneM.paths,
                                              targetMedians = target_medians,
                                              sampleID = "All_Samples_ONT",
                                              isGene = T)

tx.pseudobulk.dfs <- createPseudobulkDfList(matrixFilePaths = allSamples.ONT.sc.txM.paths,
                                            targetMedians = target_medians,
                                            sampleID = "All_Samples_ONT",
                                            isGene = F)

# Step 2: Read in pacbio data
allSamples.PB.sc.geneM.paths <- c(
  "/data/gpfs/projects/punim2251/Aim1_LongBench/ReadRarefaction_wFixedCells/data/LongBench_All/pb_sc/10Percent_FLAMES/gene_count.csv",
  "/data/gpfs/projects/punim2251/Aim1_LongBench/ReadRarefaction_wFixedCells/data/LongBench_All/pb_sc/25Percent_FLAMES/gene_count.csv",
  "/data/gpfs/projects/punim2251/Aim1_LongBench/ReadRarefaction_wFixedCells/data/LongBench_All/pb_sc/50Percent_FLAMES/gene_count.csv",
  "/data/gpfs/projects/punim2251/Aim1_LongBench/ReadRarefaction_wFixedCells/data/LongBench_All/pb_sc/75Percent_FLAMES/gene_count.csv",
  "/data/gpfs/projects/punim2251/Aim1_LongBench/ReadRarefaction_wFixedCells/data/LongBench_All/pb_sc/100Percent_FLAMES/gene_count.csv"
)

allSamples.PB.sc.txM.paths <- c(
  "/data/gpfs/projects/punim2251/Aim1_LongBench/ReadRarefaction_wFixedCells/data/LongBench_All/pb_sc/10Percent_FLAMES/allSamples_PB_sc_10Percent_iso_count.csv",
  "/data/gpfs/projects/punim2251/Aim1_LongBench/ReadRarefaction_wFixedCells/data/LongBench_All/pb_sc/25Percent_FLAMES/allSamples_PB_sc_25Percent_iso_count.csv",
  "/data/gpfs/projects/punim2251/Aim1_LongBench/ReadRarefaction_wFixedCells/data/LongBench_All/pb_sc/50Percent_FLAMES/allSamples_PB_sc_50Percent_iso_count.csv",
  "/data/gpfs/projects/punim2251/Aim1_LongBench/ReadRarefaction_wFixedCells/data/LongBench_All/pb_sc/75Percent_FLAMES/allSamples_PB_sc_75Percent_iso_count.csv",
  "/data/gpfs/projects/punim2251/Aim1_LongBench/ReadRarefaction_wFixedCells/data/LongBench_All/pb_sc/100Percent_FLAMES/allSamples_ONT_sc_100Percent_iso_count.csv"
)

gene.pseudobulk.dfs.pb <- createPseudobulkDfList(matrixFilePaths = allSamples.PB.sc.geneM.paths,
                                                 targetMedians = target_medians,
                                                 sampleID = "All_Samples_PacBio",
                                                 isGene = T)

tx.pseudobulk.dfs.pb <- createPseudobulkDfList(matrixFilePaths = allSamples.PB.sc.txM.paths,
                                               targetMedians = target_medians,
                                               sampleID = "All_Samples_PacBio",
                                               isGene = F)
## Step 2a: Export pseudobulk and bulk dfs----
saveRDS(gene.pseudobulk.dfs.pb, 
        "allCellLines_genePseudobulkDfs_pb.rds")
saveRDS(tx.pseudobulk.dfs.pb, 
        "allCellLines_isoPseudobulkDfs_pb.rds")
saveRDS(gene.pseudobulk.dfs, 
        "allCellLines_genePseudobulkDfs_ont.rds")
saveRDS(tx.pseudobulk.dfs, 
        "allCellLines_isoPseudobulkDfs_ont.rds")

# Step 3: Prepare Rarefaction Curve Dataframe Inputs -------
gene.sc.featureNums <- generateFeatureNumSummary(gene.pseudobulk.dfs) %>% 
  dplyr::mutate(sampleID = "ONT SC Genes Identified and Quantified by FLAMES",
                type = "Gene Discovery",
                rarefactionDesc = NULL)

tx.sc.featureNums <- generateFeatureNumSummary(tx.pseudobulk.dfs) %>% 
  dplyr::mutate(sampleID = "ONT SC Txs Quantified by FLAMES (w/o novel)",
                type = "Isoform Discovery",
                rarefactionDesc = NULL)

gene.sc.featureNums.pb <- generateFeatureNumSummary(gene.pseudobulk.dfs.pb) %>% 
  dplyr::mutate(sampleID = "PB SC Genes Identified and Quantified by FLAMES",
                type = "Gene Discovery",
                rarefactionDesc = NULL)

tx.sc.featureNums.pb <- generateFeatureNumSummary(tx.pseudobulk.dfs.pb) %>% 
  dplyr::mutate(sampleID = "PB SC Txs Quantified by FLAMES (w/o novel)",
                type = "Isoform Discovery",
                rarefactionDesc = NULL)

## Step 4: Plot rarefaction curve [Efficient code]------
# Function to scale the axes of plot 2 to match plot 1
matchAxisScales <- function(plt1, plt2) {
  plt.metadata <- ggplot_build(plt1)
  
  x_scale <- plt.metadata$layout$panel_scales_x[[1]]
  y_scale <- plt.metadata$layout$panel_scales_y[[1]]
  x_limits <- x_scale$range$range
  y_limits <- y_scale$range$range
  
  # Function to update scale while preserving other properties
  update_scale <- function(scale, new_limits) {
    scale$limits <- new_limits
    scale
  }
  
  # Update plt2 scales
  plt2$scales$scales <- lapply(plt2$scales$scales, function(scale) {
    if (inherits(scale, "ScaleContinuousPosition")) {
      if (scale$aesthetics[1] == "x") {
        update_scale(scale, x_limits)
      } else if (scale$aesthetics[1] == "y") {
        update_scale(scale, y_limits)
      } else {
        scale
      }
    } else {
      scale
    }
  })
  
  return(plt2)
}
plotRarefactionCurve <- function(df.list.input,
                                 plt.title = "Rarefaction Curve",
                                 export.plt = F,
                                 file.prefix = "unnamed_plt",
                                 txt.size = 14,
                                 file.width = 10,
                                 file.height = 6,
                                 removeLegend = F){
  # Combine all data frames into one
  rarefactionCurveData <- do.call(rbind, df.list.input)
  
  rarefaction.plt <- ggplot(rarefactionCurveData, aes(x = targetMedian, y = featureNum, color = sampleID)) +
    geom_line() +
    geom_point() +
    labs(x = "Target Median",
         y = "# Unique Features Detected") +
    ggtitle(plt.title) +
    theme(text = element_text(size = txt.size),
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text(angle = 45, hjust = 1),
          legend.title = element_blank()) +
    scale_color_manual(values = color_palette)
  
  if(removeLegend){
    rarefaction.plt <- rarefaction.plt +
      theme(legend.position = "none") 
  }
  if(export.plt){
    filename <- paste0(file.prefix, ".png")
    ggsave(filename, plot = rarefaction.plt, 
           width = file.width, height = file.height, 
           dpi = 300)
  }
  
  return(rarefaction.plt)
}

### plts for masters talk---------
tx.sc.featureNums <- tx.sc.featureNums %>% 
  dplyr::mutate(sampleID = "Isoforms in ONT Single-Cell Data")

#removing the wonky ont 100% point temporarily
#tx.sc.featureNums <- tx.sc.featureNums %>% slice(-5)

gene.sc.featureNums <- gene.sc.featureNums %>% 
  dplyr::mutate(sampleID = "Genes in ONT Single-Cell Data")

tx.sc.featureNums.pb <- tx.sc.featureNums.pb %>% 
  dplyr::mutate(sampleID = "Isoforms in PB Single-Cell Data")
gene.sc.featureNums.pb <- gene.sc.featureNums.pb %>% 
  dplyr::mutate(sampleID = "Genes in PB Single-Cell Data")

# Define a color palette (add colors for each unique sampleID)
color_palette <- c("Isoforms in Bulk Data\n" = "maroon",
                   "Genes in Bulk Data\n" = "#00BFC4",
                   "Isoforms in ONT Single-Cell Data" = "#00BA38",
                   "Genes in ONT Single-Cell Data" = "#619CFF",
                   "Isoforms in PB Single-Cell Data" = "#F8766D",
                   "Genes in PB Single-Cell Data" = "#B79F00") ##F564E3

masters.talk.plt <- plotRarefactionCurve(df.list.input = list(gene.sc.featureNums,
                                                              tx.sc.featureNums,
                                                              gene.sc.featureNums.pb,
                                                              tx.sc.featureNums.pb),
                                         plt.title = "",
                                         export.plt = F,
                                         file.prefix = "allSamples_sc_rarefaction_v2",
                                         removeLegend = F,
                                         txt.size = 17,
                                         file.height = 8)

masters.talk.plt
