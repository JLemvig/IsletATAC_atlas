# ==============================================================================
# 031_IsletEpigenetics_CreateClustersforDoubletfinding.R
# Purpose : Iterative LSI on peaks and Harmony integration (used for the resolution
#           20 clusters in 035), and resolution 10 clusters merged into ManuelCluster
#           for the first scDblFinder run (032).
# Inputs  : data/public/scATACseq/ArchR/IsletProject3
# Outputs : IsletProject3 updated in place (LSI, Harmony, HarmonyCluster, ManuelCluster)
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

# Add Peak LSI (seed = 1 is the ArchR default, written out for reproducibility)
IsletProject3 <- addIterativeLSI(
  ArchRProj = IsletProject3,
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
IsletProject3 <- addHarmony(
  ArchRProj = IsletProject3,
  reducedDims = "LSI",
  name = "Harmony",
  groupBy = c("Sample"),
  dimsToUse = 2:50,
  force = TRUE,
  epsilon.harmony = 1e-04,
  max.iter.harmony = 10
)

# High-resolution clusters (resolution 10). This overwrites the resolution 0.2
# HarmonyCluster from 021 (kept as HarmonyClusterR02).
IsletProject3 <- addClusters(
  input = IsletProject3,
  reducedDims = "Harmony",
  method = "Seurat",
  name = "HarmonyCluster",
  resolution = 10,
  force = TRUE,
  maxClusters = NULL,
  nOutlier = 100,
  seed = 1
)

# Combine the resolution 10 clusters into ManuelCluster, the clusters given to
# scDblFinder in 032 (ids valid for the clustering above with seed 1)
IsletProject3$ManuelCluster <- dplyr::case_when(
  IsletProject3$HarmonyCluster %in% paste0("C",c(1:9)) ~ "C1",
  IsletProject3$HarmonyCluster %in% paste0("C",c(10:19,64,70:89)) ~ "C2",
  IsletProject3$HarmonyCluster %in% paste0("C",c(20:56,61)) ~ "C3",
  IsletProject3$HarmonyCluster %in% paste0("C",c(57:60,63,65)) ~ "C5",
  IsletProject3$HarmonyCluster %in% paste0("C",c(66:69)) ~ "C6",
  IsletProject3$HarmonyCluster %in% paste0("C",c(90:98)) ~ "C7",
  IsletProject3$HarmonyCluster %in% paste0("C",c(99)) ~ "C8",
  IsletProject3$HarmonyCluster %in% paste0("C",c(100:102)) ~ "C9",
  IsletProject3$HarmonyCluster %in% paste0("C",c(103:112)) ~ "C10",
  IsletProject3$HarmonyCluster %in% paste0("C",c(113:125)) ~ "C11",
  IsletProject3$HarmonyCluster %in% paste0("C",c(126,127)) ~ "C12",
  IsletProject3$HarmonyCluster %in% paste0("C",c(128,129)) ~ "C13",
  IsletProject3$HarmonyCluster %in% paste0("C",c(130:136)) ~ "C14",
  IsletProject3$HarmonyCluster %in% paste0("C",c(137:139)) ~ "C15",
  IsletProject3$HarmonyCluster %in% paste0("C",c(140:145)) ~ "C16"
)

saveArchRProject(IsletProject3)
