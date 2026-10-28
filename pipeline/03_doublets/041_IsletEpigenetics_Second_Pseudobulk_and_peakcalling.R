# ==============================================================================
# 041_IsletEpigenetics_Second_Pseudobulk_and_peakcalling.R
# Purpose : On the doublet-filtered project: LSI on the first peak set, Harmony and
#           clustering, second pseudobulk and MACS2 peak calling, then LSI, Harmony,
#           high-resolution clusters (resolution 20) and UMAP on the new peaks, and
#           merged clusters (ManuelCluster) for the second scDblFinder run in 042.
# Inputs  : data/public/scATACseq/ArchR/IsletProject4
# Outputs : IsletProject4 updated in place (new peak set and PeakMatrix, LSI, Harmony,
#           HarmonyCluster, HarmonyClusterR20, HarmonyUMAP, ManuelCluster)
# Requires MACS2 2.2.9.1 on the PATH (found with findMacs2()).
# ==============================================================================

data_dir      <- normalizePath("data/public", mustWork = FALSE)     # public/processed data

## Load dependencies
library(BSgenome.Hsapiens.UCSC.hg38)
library(ArchR)
library(tidyverse)

## Set ArchR threads and reference genome
addArchRThreads(threads = round(parallelly::availableCores()/4))
addArchRGenome(genome = "hg38")

# Load ArchRProject
IsletProject4 <- loadArchRProject(path = file.path(data_dir, "scATACseq/ArchR/IsletProject4"))

# Add Peak LSI (seed = 1 is the ArchR default, written out for reproducibility)
IsletProject4 <- addIterativeLSI(
  ArchRProj = IsletProject4,
  useMatrix = "PeakMatrix",
  name = "LSI",
  iterations = 4,
  clusterParams = list(
    resolution = c(0.1,0.15,0.2),
    sampleCells = 100000,
    n.start = 10,
    maxClusters = 8
  ),
  sampleCellsPre = 100000,
  nPlot = 100000,
  varFeatures = 50000,
  dimsToUse = 1:50,
  seed = 1,
  force = TRUE
)

# Add Harmony integration per sample (addHarmony has no seed argument)
set.seed(1)
IsletProject4 <- addHarmony(
  ArchRProj = IsletProject4,
  reducedDims = "LSI",
  name = "Harmony",
  groupBy = c("Sample"),
  dimsToUse = 2:50,
  force = TRUE,
  epsilon.harmony = 1e-04,
  max.iter.harmony = 10
)

# Clusters (resolution 0.2) used for the second pseudobulk and peak calling
IsletProject4 <- addClusters(
  input = IsletProject4,
  reducedDims = "Harmony",
  method = "Seurat",
  name = "HarmonyCluster",
  resolution = 0.2,
  force = TRUE,
  maxClusters = NULL,
  nOutlier = 100,
  sampleCells = 100000,
  seed = 1
)
saveArchRProject(IsletProject4)

## Second Peak calling, one pseudobulk per HarmonyCluster
IsletProject4 <- addGroupCoverages(
  ArchRProj = IsletProject4,
  groupBy = "HarmonyCluster",
  force = TRUE
)
# Call peaks with Macs2
IsletProject4 <- addReproduciblePeakSet(
  ArchRProj = IsletProject4,
  groupBy = "HarmonyCluster",
  peakMethod = "Macs2",
  pathToMacs2 = findMacs2(),
  force = TRUE
)
# Add peak matrix
IsletProject4 <- addPeakMatrix(
  ArchRProj = IsletProject4,
  force = TRUE
)

## LSI on new peaks
IsletProject4 <- addIterativeLSI(
  ArchRProj = IsletProject4,
  useMatrix = "PeakMatrix",
  name = "LSI",
  iterations = 4,
  clusterParams = list(
    resolution = c(0.1,0.2,0.3),
    sampleCells = 100000,
    n.start = 10,
    maxClusters = 9
  ),
  sampleCellsPre = 100000,
  nPlot = 100000,
  varFeatures = 50000,
  dimsToUse = 1:50,
  seed = 1,
  force = TRUE
)
# Add Harmony
set.seed(1)
IsletProject4 <- addHarmony(
  ArchRProj = IsletProject4,
  reducedDims = "LSI",
  name = "Harmony",
  groupBy = c("Sample"),
  dimsToUse = 2:50,
  force = TRUE,
  epsilon.harmony = 1e-04,
  max.iter.harmony = 10
)

# High-resolution clusters (resolution 20), merged into ManuelCluster below
IsletProject4 <- addClusters(
  input = IsletProject4,
  reducedDims = "Harmony",
  method = "Seurat",
  name = "HarmonyClusterR20",
  resolution = 20,
  force = TRUE,
  maxClusters = NULL,
  nOutlier = 50,
  seed = 1
)

IsletProject4 <- addUMAP(
  ArchRProj = IsletProject4,
  reducedDims = "Harmony",
  name = "HarmonyUMAP",
  nNeighbors = 30,
  minDist = 0.5,
  metric = "cosine",
  dimsToUse = 1:49,
  seed = 1,
  force = TRUE
)

## Merge the resolution 20 clusters into ManuelCluster, the clusters given to
## scDblFinder in 042 (ids valid for the clustering above with seed 1)
IsletProject4$ManuelCluster <- case_when(
  IsletProject4$HarmonyClusterR20 %in% paste0("C",c(1:13,15,16,174:177)) ~ "C1",
  IsletProject4$HarmonyClusterR20 %in% paste0("C",c(14)) ~ "C2",
  IsletProject4$HarmonyClusterR20 %in% paste0("C",c(17:22)) ~ "C3",
  IsletProject4$HarmonyClusterR20 %in% paste0("C",c(23:88,178:182)) ~"C4",
  IsletProject4$HarmonyClusterR20 %in% paste0("C",c(89:104)) ~ "C5",
  IsletProject4$HarmonyClusterR20 %in% paste0("C",c(105:173,183:206)) ~ "C6",
  IsletProject4$HarmonyClusterR20 %in% paste0("C",c(209:220,246:251,256)) ~ "C7",
  IsletProject4$HarmonyClusterR20 %in% paste0("C",c(221:232,234:245)) ~ "C8",
  IsletProject4$HarmonyClusterR20 %in% paste0("C",c(252:255)) ~ "C9",
  IsletProject4$HarmonyClusterR20 %in% paste0("C",c(257,258)) ~ "C10",
  IsletProject4$HarmonyClusterR20 %in% paste0("C",c(259,260)) ~ "C11",
  IsletProject4$HarmonyClusterR20 %in% paste0("C",c(261)) ~ "C12",
  IsletProject4$HarmonyClusterR20 %in% paste0("C",c(262:266)) ~ "C13",
  IsletProject4$HarmonyClusterR20 %in% paste0("C",c(267:270)) ~ "C14",
  IsletProject4$HarmonyClusterR20 %in% paste0("C",c(271)) ~ "C15",
  IsletProject4$HarmonyClusterR20 %in% paste0("C",c(272:276)) ~ "C16",
  IsletProject4$HarmonyClusterR20 %in% paste0("C",c(207,208)) ~ "C17",
  IsletProject4$HarmonyClusterR20 %in% paste0("C",c(233)) ~ "C18",
  .default = as.character(IsletProject4$HarmonyCluster))

saveArchRProject(IsletProject4)
