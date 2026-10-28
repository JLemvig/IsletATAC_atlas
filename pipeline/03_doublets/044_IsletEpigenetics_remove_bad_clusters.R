# ==============================================================================
# 044_IsletEpigenetics_remove_bad_clusters.R
# Purpose : Remove the doublet and low-quality clusters chosen from the 043 overviews,
#           and combine the remaining endocrine and non-endocrine cells into the final
#           project.
# Inputs  : data/public/scATACseq/ArchR/IsletProject4, IsletProject4e, IsletProject4i
# Outputs : data/public/scATACseq/ArchR/IsletProject (final project used by all later
#           scripts)
# ==============================================================================

data_dir      <- normalizePath("data/public", mustWork = FALSE)     # public/processed data

## Load dependencies
library(BSgenome.Hsapiens.UCSC.hg38)
library(ArchR)

## Set ArchR threads and reference genome
addArchRThreads(threads = round(parallelly::availableCores()/4))
addArchRGenome(genome = "hg38")

# Load ArchRProject
IsletProject4 <- loadArchRProject(path = file.path(data_dir, "scATACseq/ArchR/IsletProject4"))
IsletProject4i <- loadArchRProject(path = file.path(data_dir, "scATACseq/ArchR/IsletProject4i"))
IsletProject4e <- loadArchRProject(path = file.path(data_dir, "scATACseq/ArchR/IsletProject4e"))

## Clusters removed, chosen from data/project/doublets_endocrine.csv and
## doublets_non_endocrine.csv (043); ids valid for the 043 clusterings with seed 1
# Endocrine (HarmonyClusterR10e)
bad_clusters_endo <- paste0("C",c(10,16))
# Non-endocrine (HarmonyClusterR10i):
#   low quality (low normalised nFrags and high blacklist ratio, with high or low TSS
#   enrichment and a FRIP or nucleosome ratio outlier, average > 0.22): C1, C31, C60, C91, C98
#   doublets (fraction of doublet calls > 0.20 for the three methods, or > 0.40 for any
#   method combined): C10, C29, C36, C88, C94
bad_clusters_non_endo <- paste0("C",c(1,60,91,98,31,10,29,36,88,94))

# Final project: the endocrine and non-endocrine cells without the bad clusters. Cells
# of IsletProject4 outside HarmonyClusterR02 C1-C12 are in neither subset and are dropped.
IsletProject5 <- subsetArchRProject(IsletProject4,
                                     cells = c(IsletProject4e$cellNames[!IsletProject4e$HarmonyClusterR10e %in% bad_clusters_endo],
                                               IsletProject4i$cellNames[!IsletProject4i$HarmonyClusterR10i %in% bad_clusters_non_endo]),
                                     dropCells = TRUE,
                                     outputDirectory = file.path(data_dir, "scATACseq/ArchR/IsletProject/"),
                                     force = TRUE)
