# ==============================================================================
# 035_IsletEpigenetics_Doubletfinding_RemoveDoubletClusters.R
# Purpose : Cluster at resolution 20, summarise the doublet calls of ArchR (034), AMULET
#           (033) and scDblFinder (032) per cluster and remove the doublet-enriched
#           clusters.
# Inputs  : data/public/scATACseq/ArchR/IsletProject3 (with the doublet calls from 032-034)
# Outputs : data/project/doublets.csv (doublet overview per HarmonyClusterR20 cluster);
#           data/public/scATACseq/ArchR/IsletProject4
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
IsletProject3 <- loadArchRProject(path = file.path(data_dir, "scATACseq/ArchR/IsletProject3"))

## Adding nFrags normalised to the sample mean to ArchRproject
df <- data.frame("cellNames" = IsletProject3$cellNames,
                 "Sample" = IsletProject3$Sample,
                 "nFrags" = IsletProject3$nFrags)
df <- df %>%
  dplyr::group_by(Sample) %>%
  dplyr::mutate(norm.nFrags = nFrags/mean(nFrags))

IsletProject3$norm.nFrags <- df$norm.nFrags

## High-resolution clusters (resolution 20) used for doublet cluster removal
## (seed = 1 is the ArchR default)
IsletProject3 <- addClusters(
  input = IsletProject3,
  reducedDims = "Harmony",
  method = "Seurat",
  name = "HarmonyClusterR20",
  resolution = 20,
  force = TRUE,
  maxClusters = NULL,
  nOutlier = 100,
  seed = 1
)
saveArchRProject(IsletProject3)

## Creating overview of cells pr HarmonyClusterR20 cluster identified as doublet, as well as
## normalized nFrags. The doublet clusters removed below were chosen from this overview.
# Doublet calls: ArchR doublet, AMULET q < 0.05, or scDblFinder doublet in 9 or 10 of 10 runs
archr_doublets <- table("Cluster" = IsletProject3$HarmonyClusterR20,
                        "ArchR.Doublets" = ifelse(IsletProject3$ArchRDoublet == 1,"ArchRDoublets",NA))
amulet_doublets <- table("Cluster" = IsletProject3$HarmonyClusterR20,
                         "Amulet.q5.Doublets" = ifelse(IsletProject3$Amulet.q.value < 0.05,"AmuletDoublets",NA))
scdblfinder_doublets <- table("Cluster" = as.factor(IsletProject3$HarmonyClusterR20),
                              "scDblFinder.Doublets" = ifelse(IsletProject3$scDblFinder %in% 9:10,"scDblFinderDoublets",NA))
any_doublets <- table("Cluster" = IsletProject3$HarmonyClusterR20,
                      "Any.Doublets" = ifelse(IsletProject3$ArchRDoublet == 1,"AnyDoublets",
                                       ifelse(IsletProject3$Amulet.q.value < 0.05,"AnyDoublets",
                                       ifelse(IsletProject3$scDblFinder %in% 9:10,"AnyDoublets",NA))))

# Fraction of cells per cluster called doublet by each method
sum <- table("Cluster" = as.factor(IsletProject3$HarmonyClusterR20))
archr_ratio <- round(as.vector(archr_doublets) / as.vector(sum),2)
amulet_ratio <- round(as.vector(amulet_doublets) / as.vector(sum),2)
scdbl_ratio <- round(as.vector(scdblfinder_doublets) / as.vector(sum),2)
any_ratio <- round(as.vector(any_doublets) / as.vector(sum),2)

nfrag <- data.frame("Cluster" = as.factor(IsletProject3$HarmonyClusterR20),
                    "norm.nFrags" = IsletProject3$norm.nFrags)
nfrag <- nfrag %>%
  group_by(Cluster) %>%
  summarise_at(vars(norm.nFrags), list(name = mean))

doublet_overview <- cbind(archr_doublets, archr_ratio, amulet_doublets, amulet_ratio,
                          scdblfinder_doublets, scdbl_ratio, any_doublets, any_ratio,
                          sum, "nFrags" = unlist(nfrag[,2]))
doublet_overview <- doublet_overview[str_order(row.names(doublet_overview), numeric = TRUE), ]
write.csv(doublet_overview, file.path(project_dir, "doublets.csv"))


## Subset ArchR project without the doublet clusters (HarmonyClusterR20 ids, valid for the
## clustering above with seed 1)
doublet_clusters <- paste0("C",c(157,204:210,233:236,238,240:245,249,253:255,268:270))
IsletProject4 <- subsetArchRProject(
  ArchRProj = IsletProject3,
  cells = IsletProject3$cellNames[!IsletProject3$HarmonyClusterR20 %in% doublet_clusters],
  outputDirectory = file.path(data_dir, "scATACseq/ArchR/IsletProject4"),
  dropCells = TRUE,
  force = TRUE
)
