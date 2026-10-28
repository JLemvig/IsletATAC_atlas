# ==============================================================================
# 052_isletepigenetics_PairedInsertionCounting.R
# Purpose : Paired-insertion counting (PICsnATAC) of each sample's fragments in the
#           project peak set, written into the sample ArrowFile as matrix "PIC" (used
#           for the AnnData export in 053 and later analyses).
# Inputs  : data/public/scATACseq/ArchR/IsletProject;
#           data/public/scATACseq/<study>/Preprocessed/<sample>/frags.sort.bed.gz;
#           data/sensitive/snATACseq/FUSION/Preprocessed/<sample>/frags.sort.bed.gz
# Outputs : "PIC" matrix in each ArrowFile of IsletProject
# ==============================================================================

data_dir      <- normalizePath("data/public", mustWork = FALSE)     # public/processed data
sensitive_dir <- normalizePath("data/sensitive", mustWork = FALSE)  # controlled-access data (FUSION)

## Load dependencies
library(BSgenome.Hsapiens.UCSC.hg38)
library(ArchR)
library(tidyverse)
library(data.table)
library(GenomicRanges)
library(Matrix)
library(PICsnATAC)

## Set ArchR threads and reference genome
addArchRThreads(threads = round(parallelly::availableCores()/8))
addArchRGenome(genome = "hg38")

# Load ArchRProject
IsletProject <- loadArchRProject(path = file.path(data_dir, "scATACseq/ArchR/IsletProject/"))

# Fragment files of the samples in the project
Study <- c("HPAP","Chiou","FUSION","Wang","Benaglio","Thibodeau","Augsornworawat","Weng")
inputFiles <- c()
for(i in Study){
  folder <- ifelse(i == "FUSION",file.path(sensitive_dir,"snATACseq/"),file.path(data_dir,"scATACseq/"))
  dir <- paste0(folder,i,"/Preprocessed/")
  Samples <- list.files(dir)
  for(j in Samples){
    if(j %in% unique(IsletProject$Sample)) {
      input <- paste0(dir,j,"/frags.sort.bed.gz")
      names(input) = j
      inputFiles <- c(inputFiles, input)
    }
  }
}

pic_mat <- rep(list(), length = length(inputFiles))
for (s in seq_along(inputFiles)) {
  ## select one sample
  sample_sel <- names(inputFiles)[s]
  allCells <- rownames(getCellColData(IsletProject)) ## all cells after filtering
  cells_sample <- allCells[startsWith(allCells, paste0(sample_sel,"#"))]
  cells_sample_name_trunc <- sub(".*\\#", "", cells_sample)

  peak_sets <- IsletProject@peakSet ## the peak set, from MACS2 peak calling results
  names(peak_sets) <- NULL ## remove the peakset names
  fragment_tsv_gz_file_location <- inputFiles[s]

  pic_mat[[s]] <- PIC_counting(
    cells = cells_sample_name_trunc,
    fragment_tsv_gz_file_location = fragment_tsv_gz_file_location,
    peak_sets = peak_sets
  )

  ## change back the cell names to match with original input
  colnames(pic_mat[[s]]) <- cells_sample
  cat(paste0("Done with: ",sample_sel,", ",s," of ",length(inputFiles)))
}

## Add PIC matrix to arrowfiles (same features and cell order as the PeakMatrix), using
## ArchR internal functions
ArrowFiles <- paste0(file.path(data_dir,"scATACseq/ArchR/IsletProject/ArrowFiles/"),names(inputFiles),".arrow")
for (s in seq_along(inputFiles)) {
    ArrowFile <- ArrowFiles[s]
    sampleName <- ArchR:::.sampleName(ArrowFile)
    o <- h5closeAll()
    o <- ArchR:::.createArrowGroup(ArrowFile = ArrowFile, group = "PIC", force = TRUE, logFile = NULL)
    PM_tmp <- getMatrixFromArrow(ArrowFile, "PeakMatrix")
    orderedCells = colnames(PM_tmp)
    pic_mat[[s]] <- pic_mat[[s]][,match(orderedCells, colnames(pic_mat[[s]]))]
    cellNames <- colnames(pic_mat[[s]])
    dfParams <- data.frame(
      ceiling = 1000000,
      binarize = FALSE * 1,
      stringsAsFactors = FALSE)
    Class = "integer"
    Units = "Counts"
    featureDF = ArchR:::.getFeatureDF(ArrowFile, subGroup = "PeakMatrix")
    o <- ArchR:::.initializeMat(
      ArrowFile = ArrowFile,
      Group = "PIC",
      Class = Class,
      Units = Units,
      cellNames = cellNames,
      params = dfParams,
      featureDF = featureDF,
      force = TRUE
    )
    for (chr in unique(featureDF$seqnames)) {
      o <- ArchR:::.addMatToArrow(
        mat = pic_mat[[s]][ which(featureDF$seqnames %in% chr),],
        ArrowFile = ArrowFile,
        Group = paste0("PIC","/", chr),
        binarize = FALSE,
        addColSums = TRUE,
        addRowSums = TRUE
      )
  }
}
