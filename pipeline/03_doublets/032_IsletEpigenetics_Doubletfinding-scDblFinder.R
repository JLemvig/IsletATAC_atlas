# ==============================================================================
# 032_IsletEpigenetics_Doubletfinding-scDblFinder.R
# Purpose : Run scDblFinder 10 times per sample and add the number of doublet calls
#           per cell (0-10) to IsletProject3.
# Inputs  : data/public/scATACseq/ArchR/IsletProject3
# Outputs : data/public/scATACseq/ArchR/IsletProjectA-D (subsets, only to limit memory);
#           data/public/scATACseq/Signac/scDblFinder/<sample>.csv and scDblFinder.log;
#           IsletProject3$scDblFinder
# ==============================================================================

data_dir      <- normalizePath("data/public", mustWork = FALSE)     # public/processed data

## Load dependencies
library(Seurat)
library(scDblFinder)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ArchR)
library(Signac)

## Set ArchR threads and reference genome
addArchRThreads(threads = round(parallelly::availableCores()/8))
addArchRGenome(genome = "hg38")

# Load ArchRProject
IsletProject3 <- loadArchRProject(path = file.path(data_dir, "scATACseq/ArchR/IsletProject3"))

out_dir <- file.path(data_dir, "scATACseq/Signac/scDblFinder")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# The samples are split into four subsets only to limit memory use; scDblFinder is run
# per sample, so the split does not affect the result
Samples <- unique(IsletProject3$Sample)

{
  IsletProjectA <- subsetArchRProject(
    ArchRProj = IsletProject3,
    cells = IsletProject3$cellNames[IsletProject3$Sample %in% Samples[1:48]],
    outputDirectory = file.path(data_dir, "scATACseq/ArchR/IsletProjectA"),
    dropCells = TRUE,
    force = TRUE
  )

  IsletProjectB <- subsetArchRProject(
    ArchRProj = IsletProject3,
    cells = IsletProject3$cellNames[IsletProject3$Sample %in% Samples[49:67]],
    outputDirectory = file.path(data_dir, "scATACseq/ArchR/IsletProjectB"),
    dropCells = TRUE,
    force = TRUE
  )

  IsletProjectC <- subsetArchRProject(
    ArchRProj = IsletProject3,
    cells = IsletProject3$cellNames[IsletProject3$Sample %in% Samples[68:94]],
    outputDirectory = file.path(data_dir, "scATACseq/ArchR/IsletProjectC"),
    dropCells = TRUE,
    force = TRUE
  )

  IsletProjectD <- subsetArchRProject(
    ArchRProj = IsletProject3,
    cells = IsletProject3$cellNames[IsletProject3$Sample %in% Samples[95:122]],
    outputDirectory = file.path(data_dir, "scATACseq/ArchR/IsletProjectD"),
    dropCells = TRUE,
    force = TRUE
  )
}

## Create firstly Seurat object, then SingleCellExperiment, and then run scDblFinder.
set.seed(1)
{
  subsets <- c("A","B","C","D")
  logr::log_open(file.path(out_dir, "scDblFinder.log"), logdir = FALSE, compact = TRUE, show_notes = FALSE)
  for(i in subsets){
    IP <- loadArchRProject(paste0(file.path(data_dir, "scATACseq/ArchR/IsletProject"),i,"/"), showLogo = FALSE)
    PM <- getMatrixFromProject(ArchRProj = IP, useMatrix = "PeakMatrix", binarize = FALSE)
    rm(IP)
    chrom_assay <- CreateChromatinAssay(PM@assays@data$PeakMatrix, validate.fragments=FALSE, ranges = PM@rowRanges)
    metadata <- as.data.frame(PM@colData)
    seurat <- CreateSeuratObject(counts = chrom_assay, assay = "peaks", meta.data = metadata)
    rm(PM)
    rm(chrom_assay)
    rm(metadata)
    logr::log_print(paste0(Sys.time()," Seurat object ",i," created"))
    #Create singlecellexperiment
    sce <- as.SingleCellExperiment(seurat, assay = "peaks")
    rm(seurat)
    # Clusters given to scDblFinder: the merged clusters (ManuelCluster) from 031
    sce$clusters <- as.factor(as.numeric(sub("C","",sce$ManuelCluster)))
    sce.list <- list()
    #Subset each sample to an unique singlecellexperiment
    counter <- 1
    for (ds in unique(sce$Sample)) {
      sce.list[[counter]] <- sce[ ,sce$Sample == ds]
      counter <- counter + 1
    }
    rm(sce)
    rm(counter)
    #Run scDblfinder 10 times for each singlecellexperiement
    for (ds in 1:length(sce.list)) {
      # Setup to save stats
      doub.stats <- matrix(nrow = dim(sce.list[[ds]])[2], ncol = 10, 0)
      rownames(doub.stats) <- colnames(sce.list[[ds]])
      # Run scDblFinder; error tolerant, a failed run is repeated until 10 runs have
      # succeeded (at most 50 attempts per sample)
      rep <- 1
      attempts <- 0
      while (rep <= 10) {
        attempts <- attempts + 1
        if (attempts > 50) stop("scDblFinder failed too often for sample ", unique(sce.list[[ds]]$Sample))
        doub <- tryCatch({scDblFinder(sce = sce.list[[ds]],
                                      clusters = sce.list[[ds]]$clusters,
                                      aggregateFeatures=TRUE,
                                      nfeatures=25,
                                      processing="normFeatures",
                                      dbr.sd = 1)},
                         error=function(e){cat("ERROR")})
        if (class(doub$scDblFinder.score) == "numeric") {
          # The class factor is stored as its code: singlet = 1, doublet = 2
          doub.stats[,rep] <- doub$scDblFinder.class
          print(rep)
          rep <- rep + 1
          rm(doub)
        }
      }
      write.csv(
        x = doub.stats,
        file = file.path(out_dir, paste0(unique(sce.list[[ds]]$Sample), ".csv")),
        row.names = TRUE
      )
      logr::log_print(paste0(Sys.time()," Done with: ",unique(sce.list[[ds]]$Sample)," - ",ds,":",length(sce.list)))
    }
  }
  logr::log_close()
}

## Add the number of doublet calls per cell to the project
doub <- NULL
for(i in unique(IsletProject3$Sample)){
  csv <- read.csv(file.path(out_dir, paste0(i, ".csv")))
  doub <- dplyr::bind_rows(doub, csv)
}
colnames(doub)[1] <- "cellNames"
# Sum of the 10 codes minus 10 = number of runs that called the cell a doublet (0-10)
doub[,12] <- rowSums(doub[,2:11]) - 10
df <- dplyr::left_join(data.frame("cellNames" = IsletProject3$cellNames),doub[,c(1,12)], by = "cellNames")
IsletProject3$scDblFinder <- df[,2]

saveArchRProject(IsletProject3)
