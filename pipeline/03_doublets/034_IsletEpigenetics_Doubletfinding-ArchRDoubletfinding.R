# ==============================================================================
# 034_IsletEpigenetics_Doubletfinding-ArchRDoubletfinding.R
# Purpose : Run ArchR doublet scoring and filtering, and mark the cells ArchR would
#           remove as doublets in IsletProject3 (cells are not removed here).
# Inputs  : data/public/scATACseq/ArchR/IsletProject3
# Outputs : IsletProject3$ArchRDoublet (1 = doublet, 0 = not)
# ==============================================================================

data_dir      <- normalizePath("data/public", mustWork = FALSE)     # public/processed data

## Load dependencies
library(BSgenome.Hsapiens.UCSC.hg38)
library(ArchR)

## Set ArchR threads and reference genome
addArchRThreads(threads = round(parallelly::availableCores()/4))
addArchRGenome(genome = "hg38")

# Load ArchRProject
IsletProject3 <- loadArchRProject(path = file.path(data_dir, "scATACseq/ArchR/IsletProject3"))

# Adding ArchRDoublet scores to cells
IsletProject3 <- addDoubletScores(
  IsletProject3,
  useMatrix = "TileMatrix",
  dimsToUse = 1:50,
  nTrials = 10
)

# Filtering ArchR doublets (in memory only, nothing is saved for this object)
IsletProject3_filtered <- filterDoublets(
  ArchRProj = IsletProject3,
  filterRatio = 1
)

# Marking the cells filtered by ArchR doublet filtering on the original ArchR Project
all_cells <- data.frame("cellNames" = IsletProject3$cellNames, "ArchRDoublet" = 1)
filtered_cells <- data.frame("cellNames" = IsletProject3_filtered$cellNames, "ArchRDoublet" = 0)
archr_doublets <- dplyr::rows_update(all_cells,filtered_cells, by = "cellNames")
IsletProject3$ArchRDoublet <- archr_doublets$ArchRDoublet

saveArchRProject(IsletProject3)
