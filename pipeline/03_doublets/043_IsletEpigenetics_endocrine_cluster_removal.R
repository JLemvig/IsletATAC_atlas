# ==============================================================================
# 043_IsletEpigenetics_endocrine_cluster_removal.R
# Purpose : Subset the endocrine cells (HarmonyClusterR02 C1-C6), recluster them
#           (resolution 10), and summarise doublet calls and QC metrics per cluster to
#           choose the clusters removed in 044.
# Inputs  : data/public/scATACseq/ArchR/IsletProject4 (with the doublet calls from 032-034
#           and 042)
# Outputs : data/public/scATACseq/ArchR/IsletProject4e;
#           data/project/doublets_endocrine.csv (doublet and QC overview per cluster)
# ==============================================================================

data_dir      <- normalizePath("data/public", mustWork = FALSE)     # public/processed data
project_dir   <- normalizePath("data/project", mustWork = FALSE)    # small project files

## Load dependencies
library(BSgenome.Hsapiens.UCSC.hg38)
library(dplyr)
library(ArchR)
library(tidyverse)

## Set ArchR threads and reference genome
addArchRThreads(threads = round(parallelly::availableCores()/4))
addArchRGenome(genome = "hg38")

# Load ArchRProject
IsletProject4 <- loadArchRProject(path = file.path(data_dir, "scATACseq/ArchR/IsletProject4"))

# Endocrine cells: broad clusters C1-C6 of the resolution 0.2 clustering from 021
endocrine_cells <- paste0("C",c(1:6))

IsletProject4e <- subsetArchRProject(IsletProject4,
                                     cells = IsletProject4$cellNames[IsletProject4$HarmonyClusterR02 %in% endocrine_cells],
                                     dropCells = TRUE,
                                     outputDirectory = file.path(data_dir, "scATACseq/ArchR/IsletProject4e/"),
                                     force = TRUE)

IsletProject4e <- loadArchRProject(path = file.path(data_dir, "scATACseq/ArchR/IsletProject4e"))

IsletProject4e <- addIterativeLSI(
  ArchRProj = IsletProject4e,
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

# addHarmony has no seed argument
set.seed(1)
IsletProject4e <- addHarmony(
  ArchRProj = IsletProject4e,
  reducedDims = "LSI",
  name = "Harmony",
  groupBy = c("Sample"),
  dimsToUse = 2:50,
  force = TRUE,
  epsilon.harmony = 1e-04,
  max.iter.harmony = 10
)

IsletProject4e <- addClusters(
  input = IsletProject4e,
  reducedDims = "Harmony",
  method = "Seurat",
  name = "HarmonyClusterR10e",
  resolution = 10,
  force = TRUE,
  maxClusters = NULL,
  nOutlier = 50,
  seed = 1
)

IsletProject4e <- addUMAP(
  ArchRProj = IsletProject4e,
  reducedDims = "Harmony",
  name = "HarmonyUMAP",
  nNeighbors = 30,
  minDist = 0.5,
  metric = "cosine",
  dimsToUse = 1:49,
  seed = 1,
  force = TRUE
)

saveArchRProject(IsletProject4e)

df <- data.frame("cellNames" = IsletProject4e$cellNames,
                 "Sample" = IsletProject4e$Sample,
                 "Cluster" = IsletProject4e$HarmonyClusterR10e,
                 "nFrags" = IsletProject4e$nFrags,
                 "tss_enrichment" = IsletProject4e$TSSEnrichment,
                 "blacklist_ratio" = IsletProject4e$BlacklistRatio,
                 "nucleosome_ratio" = IsletProject4e$NucleosomeRatio,
                 "frip" = IsletProject4e$FRIP)
df <- df %>%
  dplyr::group_by(Sample) %>%
  dplyr::mutate(norm_nFrags = nFrags/mean(nFrags)) %>%
  dplyr::mutate(norm_tss_enrichment = tss_enrichment/mean(tss_enrichment)) %>%
  dplyr::mutate(norm_blacklist_ratio = blacklist_ratio/mean(blacklist_ratio)) %>%
  dplyr::mutate(norm_nucleosome_ratio = nucleosome_ratio/mean(nucleosome_ratio)) %>%
  dplyr::mutate(norm_frip = frip/mean(frip))

## Creating overview of cells pr cluster identified as doublet, as well as QC metrics
## normalised to the sample mean.
# Doublet calls: ArchR doublet, AMULET q < 0.05, or scDblFinder doublet in 9 or 10 of
# 10 runs (first run 032, second run 042)
archr_doublets <- table("Cluster" = IsletProject4e$HarmonyClusterR10e,
                        "ArchR.Doublets" = ifelse(IsletProject4e$ArchRDoublet == 1,"ArchRDoublets",NA))
amulet_doublets <- table("Cluster" = IsletProject4e$HarmonyClusterR10e,
                         "Amulet.q5.Doublets" = ifelse(IsletProject4e$Amulet.q.value < 0.05,"AmuletDoublets",NA))
scdblfinder_doublets <- table("Cluster" = as.factor(IsletProject4e$HarmonyClusterR10e),
                              "scDblFinder.Doublets" = ifelse(IsletProject4e$scDblFinder %in% 9:10,"scDblFinderDoublets",NA))
scdblfinder2_doublets <- table("Cluster" = as.factor(IsletProject4e$HarmonyClusterR10e),
                               "scDblFinder.Doublets" = ifelse(IsletProject4e$scDblFinder2 %in% 9:10,"scDblFinderDoublets2nd",NA))
any_doublets <- table("Cluster" = IsletProject4e$HarmonyClusterR10e,
                      "Any.Doublets" = ifelse(IsletProject4e$ArchRDoublet == 1,"AnyDoublets",
                                              ifelse(IsletProject4e$Amulet.q.value < 0.05,"AnyDoublets",
                                                     ifelse(IsletProject4e$scDblFinder %in% 9:10,"AnyDoublets",
                                                            ifelse(IsletProject4e$scDblFinder2 %in% 9:10,"AnyDoublets",NA)))))

sum <- table("Cluster" = as.factor(IsletProject4e$HarmonyClusterR10e))
archr_ratio <- round(as.vector(archr_doublets) / as.vector(sum),2)
amulet_ratio <- round(as.vector(amulet_doublets) / as.vector(sum),2)
scdbl_ratio <- round(as.vector(scdblfinder_doublets) / as.vector(sum),2)
scdbl2_ratio <- round(as.vector(scdblfinder2_doublets) / as.vector(sum),2)
any_ratio <- round(as.vector(any_doublets) / as.vector(sum),2)

norm_nfrag <- df %>%
  group_by(Cluster) %>%
  summarise_at(vars(norm_nFrags), list(norm.nFrags = mean))

norm_tss_enrichment <- df %>%
  group_by(Cluster) %>%
  summarise_at(vars(norm_tss_enrichment), list(norm_tss_enrichment = mean))

norm_blacklist_ratio <- df %>%
  group_by(Cluster) %>%
  summarise_at(vars(norm_blacklist_ratio), list(norm_blacklist_ratio = mean))

norm_nucleosome_ratio <- df %>%
  group_by(Cluster) %>%
  summarise_at(vars(norm_nucleosome_ratio), list(norm_nucleosome_ratio = mean))

norm_frip <- df %>%
  group_by(Cluster) %>%
  summarise_at(vars(norm_frip), list(norm_frip = mean))

doublet_overview <- cbind(archr_doublets,
                          archr_ratio,
                          amulet_doublets,
                          amulet_ratio,
                          scdblfinder_doublets,
                          scdbl_ratio,
                          scdblfinder2_doublets,
                          scdbl2_ratio,
                          any_doublets,
                          any_ratio,
                          sum,
                          "norm.nfrags" = unlist(norm_nfrag[,2]),
                          "norm_tss_enrichment" = unlist(norm_tss_enrichment[,2]),
                          "norm_blacklist_ratio" = unlist(norm_blacklist_ratio[,2]),
                          "norm_nucleosome_ratio" = unlist(norm_nucleosome_ratio[,2]),
                          "norm_frip" = unlist(norm_frip[,2]))

doublet_overview <- doublet_overview[str_order(row.names(doublet_overview), numeric = TRUE), ]

write.csv(doublet_overview, file.path(project_dir, "doublets_endocrine.csv"))

# The clusters to remove were chosen from this overview; they are listed in 044
# (bad_clusters_endo).
