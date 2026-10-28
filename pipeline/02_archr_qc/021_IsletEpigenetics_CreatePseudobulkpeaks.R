# ==============================================================================
# 021_IsletEpigenetics_CreatePseudobulkpeaks.R
# Purpose : Iterative LSI on tiles, Harmony integration, clustering and UMAP, then
#           pseudobulk coverages per cluster and the first MACS2 peak set.
# Inputs  : data/public/scATACseq/ArchR/IsletProject2
# Outputs : IsletProject2 updated in place (LSI, Harmony, HarmonyCluster,
#           HarmonyClusterR02, HarmonyUMAP, group coverages, peak set, PeakMatrix)
# Requires MACS2 2.2.9.1 on the PATH (found with findMacs2()).
# ==============================================================================

data_dir      <- normalizePath("data/public", mustWork = FALSE)     # public/processed data

## Load dependencies
library(BSgenome.Hsapiens.UCSC.hg38)
library(ArchR)

## Set ArchR threads and reference genome
addArchRThreads(threads = round(parallelly::availableCores()/4))
addArchRGenome(genome = "hg38")

# Load ArchRProject
IsletProject2 <- loadArchRProject(path = file.path(data_dir, "scATACseq/ArchR/IsletProject2"))

# Add Tile LSI (seed = 1 is the ArchR default, written out for reproducibility)
IsletProject2 <- addIterativeLSI(
  ArchRProj = IsletProject2,
  useMatrix = "TileMatrix",
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

# Add Harmony integration per sample; LSI component 1 is not used, so Harmony has
# 49 dimensions (addHarmony has no seed argument)
set.seed(1)
IsletProject2 <- addHarmony(
  ArchRProj = IsletProject2,
  reducedDims = "LSI",
  name = "Harmony",
  groupBy = c("Sample"),
  dimsToUse = 2:50,
  force = TRUE,
  epsilon.harmony = 1e-04,
  max.iter.harmony = 10
)
# Identify Clusters (resolution 0.2)
IsletProject2 <- addClusters(
  input = IsletProject2,
  reducedDims = "Harmony",
  method = "Seurat",
  name = "HarmonyCluster",
  resolution = 0.2,
  force = TRUE,
  maxClusters = NULL,
  nOutlier = 200,
  sampleCells = 100000,
  seed = 1
)
# Keep these resolution 0.2 clusters as HarmonyClusterR02: HarmonyCluster is overwritten
# in 031, and HarmonyClusterR02 is used in 043 and 04_annotation
IsletProject2$HarmonyClusterR02 <- IsletProject2$HarmonyCluster
# Add UMAP embedding on all 49 Harmony dimensions
IsletProject2 <- addUMAP(
  ArchRProj = IsletProject2,
  reducedDims = "Harmony",
  name = "HarmonyUMAP",
  nNeighbors = 30,
  minDist = 0.5,
  metric = "cosine",
  dimsToUse = 1:49,
  seed = 1,
  force = TRUE
)

# Save before the coverages and peak calling
saveArchRProject(IsletProject2)

# Add Group Coverages, one pseudobulk per Harmony cluster
IsletProject2 <- addGroupCoverages(
  ArchRProj = IsletProject2,
  groupBy = "HarmonyCluster"
)
# Call peaks with Macs2 per Harmony cluster
IsletProject2 <- addReproduciblePeakSet(
  ArchRProj = IsletProject2,
  groupBy = "HarmonyCluster",
  peakMethod = "Macs2",
  pathToMacs2 = findMacs2()
)
# Add peak matrix (also adds the FRIP column used in 022 and 023)
IsletProject2 <- addPeakMatrix(
  ArchRProj = IsletProject2
)

saveArchRProject(IsletProject2)
