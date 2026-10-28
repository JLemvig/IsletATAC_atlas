# ==============================================================================
# 051_isletepigenetics_labeling.R
# Purpose : Rerun LSI, Harmony, clustering (resolution 0.2 HarmonyClusterR02 and
#           resolution 20 HarmonyClusterR20) and UMAP on the final project, and label
#           the clusters with cell types from gene score markers. Clusters C1, C2, C5,
#           C8 and C11 are subclustered (HarmonyClusterR05_subset) and labelled
#           separately, then merged back into IsletProject$labels.
# Inputs  : data/public/scATACseq/ArchR/IsletProject
# Outputs : IsletProject updated in place (LSI, Harmony, HarmonyClusterR02,
#           HarmonyClusterR20, HarmonyUMAP, labels);
#           data/public/scATACseq/ArchR/IsletProject_subset
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
IsletProject <- loadArchRProject(path = file.path(data_dir, "scATACseq/ArchR/IsletProject/"))

## LSI, Harmony, clusters and UMAP on the final project
IsletProject <- addIterativeLSI(
  ArchRProj = IsletProject,
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

# addHarmony has no seed argument
set.seed(1)
IsletProject <- addHarmony(
  ArchRProj = IsletProject,
  reducedDims = "LSI",
  name = "Harmony",
  groupBy = c("Sample"),
  dimsToUse = 2:50,
  force = TRUE,
  epsilon.harmony = 1e-04,
  max.iter.harmony = 10
)

# Clusters labelled below (replaces the HarmonyClusterR02 from 021)
IsletProject <- addClusters(
  input = IsletProject,
  reducedDims = "Harmony",
  method = "Seurat",
  name = "HarmonyClusterR02",
  resolution = 0.2,
  force = TRUE,
  maxClusters = NULL,
  nOutlier = 100,
  sampleCells = 100000,
  seed = 1
)

# High-resolution clusters (resolution 20), used in maping_clusters_to_labels.R
IsletProject <- addClusters(
  input = IsletProject,
  reducedDims = "Harmony",
  method = "Seurat",
  name = "HarmonyClusterR20",
  resolution = 20,
  force = TRUE,
  maxClusters = NULL,
  nOutlier = 50,
  seed = 1
)

IsletProject <- addUMAP(
  ArchRProj = IsletProject,
  reducedDims = "Harmony",
  name = "HarmonyUMAP",
  nNeighbors = 30,
  minDist = 0.5,
  metric = "cosine",
  dimsToUse = 1:49,
  seed = 1,
  force = TRUE
)

## Subcluster the clusters that contain several cell types
IsletProject_subset <- subsetArchRProject(IsletProject,
                                    cells = IsletProject$cellNames[IsletProject$HarmonyClusterR02 %in% c("C1","C2","C5","C8","C11")],
                                    dropCells = TRUE,
                                    outputDirectory = file.path(data_dir, "scATACseq/ArchR/IsletProject_subset/"),
                                    force = TRUE)

IsletProject_subset <- addIterativeLSI(
  ArchRProj = IsletProject_subset,
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

set.seed(1)
IsletProject_subset <- addHarmony(
  ArchRProj = IsletProject_subset,
  reducedDims = "LSI",
  name = "Harmony",
  groupBy = c("Sample"),
  dimsToUse = 2:50,
  force = TRUE,
  epsilon.harmony = 1e-04,
  max.iter.harmony = 10
)

IsletProject_subset <- addClusters(
  input = IsletProject_subset,
  reducedDims = "Harmony",
  method = "Seurat",
  name = "HarmonyClusterR05_subset",
  resolution = 0.5,
  force = TRUE,
  maxClusters = NULL,
  nOutlier = 50,
  seed = 1
)

## Labels, assigned from gene score markers of each cluster
# Subclusters
IsletProject_subset$labels <- case_when(
  IsletProject_subset$HarmonyClusterR05_subset %in% "C1" ~ "schwann",
  IsletProject_subset$HarmonyClusterR05_subset %in% "C2" ~ "quiescent_stellate",
  IsletProject_subset$HarmonyClusterR05_subset %in% "C3" ~ "activated_stellate",
  IsletProject_subset$HarmonyClusterR05_subset %in% "C4" ~ "epsilon",
  IsletProject_subset$HarmonyClusterR05_subset %in% "C5" ~ "beta",
  IsletProject_subset$HarmonyClusterR05_subset %in% "C6" ~ "pp",
  IsletProject_subset$HarmonyClusterR05_subset %in% "C7" ~ "alpha",
  IsletProject_subset$HarmonyClusterR05_subset %in% "C8" ~ "lymphoid",
  IsletProject_subset$HarmonyClusterR05_subset %in% "C9" ~ "myeloid",
  IsletProject_subset$HarmonyClusterR05_subset %in% "C10" ~ "mast")

# Main clusters; C1, C2, C5, C8 and C11 are replaced by the subcluster labels below
IsletProject$labels <- case_when(
  IsletProject$HarmonyClusterR02 %in% "C1" ~ "C1",
  IsletProject$HarmonyClusterR02 %in% "C2" ~ "C2",
  IsletProject$HarmonyClusterR02 %in% "C3" ~ "beta",
  IsletProject$HarmonyClusterR02 %in% "C4" ~ "delta",
  IsletProject$HarmonyClusterR02 %in% "C5" ~ "C5",
  IsletProject$HarmonyClusterR02 %in% "C6" ~ "cycling",
  IsletProject$HarmonyClusterR02 %in% "C7" ~ "alpha",
  IsletProject$HarmonyClusterR02 %in% "C8" ~ "C8",
  IsletProject$HarmonyClusterR02 %in% "C9" ~ "ductal",
  IsletProject$HarmonyClusterR02 %in% "C10" ~ "acinar",
  IsletProject$HarmonyClusterR02 %in% "C11" ~ "C11",
  IsletProject$HarmonyClusterR02 %in% "C12" ~ "endothelial"
)

IsletProject$labels <- dplyr::rows_update(x = data.frame("cellnames" = IsletProject$cellNames, "labels" = IsletProject$labels),
                   y = data.frame("cellnames" = IsletProject_subset$cellNames, "labels" = IsletProject_subset$labels), by = "cellnames")[,2]

saveArchRProject(IsletProject)
saveArchRProject(IsletProject_subset)
