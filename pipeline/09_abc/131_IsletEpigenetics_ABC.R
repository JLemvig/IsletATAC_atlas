# ==============================================================================
# 131_IsletEpigenetics_ABC.R
# Purpose : Prepare the inputs for ABC (Activity-by-Contact) peak-to-gene linking: the
#           full peak set, fragment files per sample and cell type, and custom peak
#           sets (top 150,000 peaks) for all cell types other than alpha and beta.
#           ABC itself was run with the standard Snakemake workflow of the ABC
#           repository (ABC-Enhancer-Gene-Prediction) on these inputs, and its output
#           was used as is.
# Inputs  : data/public/scATACseq/ArchR/IsletProject (with labels and samplelabels)
# Outputs : data/public/scATACseq/abc/peakset.bed;
#           fragment files per samplelabels group in IsletProject/GroupFragments;
#           data/public/scATACseq/ABC/config/<celltype>_peaks.bed
# ==============================================================================

data_dir      <- normalizePath("data/public", mustWork = FALSE)     # public/processed data

library(ArchR)
library(rtracklayer)

addArchRThreads(threads = round(parallelly::availableCores()/8))

proj <- loadArchRProject(file.path(data_dir, "scATACseq/ArchR/IsletProject/"))

# Full peak set (all cell types), used as the alpha and beta peak sets for ABC
peakSet <- getPeakSet(proj)

dir.create(file.path(data_dir, "scATACseq/abc"), recursive = TRUE, showWarnings = FALSE)
export.bed(peakSet, file.path(data_dir, "scATACseq/abc/peakset.bed"))

# Fragment files per sample and cell type (samplelabels), written to the project folder
frags <- getGroupFragments(proj, groupBy = "samplelabels")

### Create peaksets for ABC-custom-peaks
# All cell types other than alpha and beta: top 150,000 peaks by fraction of cells
# in which the peak is open (PIC matrix)
dir.create(file.path(data_dir, "scATACseq/ABC/config"), recursive = TRUE, showWarnings = FALSE)
for(celltype in unique(proj$labels)){
  if(!(celltype == "alpha" | celltype == "beta")){
    proj_subset <- subsetCells(proj, cellNames = proj$cellNames[proj$labels %in% celltype])
    pm <- getMatrixFromProject(ArchRProj = proj_subset,
                               useMatrix = "PIC",
                               verbose = FALSE,
                               binarize = FALSE)
    peakmatrix <- pm@assays@data@listData$PIC
    open_peak_fraction <- rowSums(peakmatrix > 0) / ncol(peakmatrix)
    rm(pm)
    rm(peakmatrix)
    top_peaks <- rank(open_peak_fraction, ties.method = "min") > (length(open_peak_fraction) - 150000)
    peakset <- getPeakSet(ArchRProj = proj_subset)
    peakset_top <- peakset[top_peaks, ]
    export(peakset_top, con = paste0(file.path(data_dir, "scATACseq/ABC/config/"),celltype,"_peaks.bed"), format = "BED")
  }
}
