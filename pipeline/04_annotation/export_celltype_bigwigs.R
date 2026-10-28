# ==============================================================================
# export_celltype_bigwigs.R
# Purpose : Export per cell type BigWig coverage tracks (ReadsInTSS normalised, tile size
#           100 and 200 bp, at most 1,000 cells per group) for the browser tracks.
# Inputs  : data/public/scATACseq/ArchR/IsletProject (with labels from 051)
# Outputs : data/public/scATACseq/ArchR/bigwig/tile_<size>_maxcells_1000_ceiling_4/<celltype>_*.bw
# ==============================================================================

data_dir      <- normalizePath("data/public", mustWork = FALSE)     # public/processed data

## Load dependencies
set.seed(1)
options(scipen=999)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ArchR)
library(tidyverse)

addArchRThreads(threads = round(parallelly::availableCores()/2))
addArchRGenome(genome = "hg38")

IsletProject <- loadArchRProject(path = file.path(data_dir, "scATACseq/ArchR/IsletProject/"))

exportGroupBW <- function(ArchRProj, groupBy = "labels",
                          normMethod = "ReadsInTSS",
                          tile_size = 100, max_cells = 1000, ceiling = 4,
                          in_dir, out_dir) {

  # Run getGroupBW
  bw_files <- getGroupBW(
    ArchRProj = ArchRProj,
    groupBy = groupBy,
    normMethod = normMethod,
    tileSize = tile_size,
    maxCells = max_cells,
    ceiling = ceiling
  )

  # Create new folder in out_dir with param info
  folder_name <- paste0("tile_", tile_size,
                        "_maxcells_", max_cells,
                        "_ceiling_", ceiling)
  target_dir <- file.path(out_dir, folder_name)
  dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)

  # Find BigWig files in input dir
  bw_list <- list.files(in_dir, pattern = "\\.bw$", full.names = TRUE)

  # Move & rename files
  for (f in bw_list) {
    base <- basename(f) # e.g. "alpha.bw"
    prefix <- stringr::str_split(base,pattern = "-")[[1]][1]
    new_name <- paste0(prefix,
                       "_tile_", tile_size,
                       "_maxcells_", max_cells,
                       "_ceiling_", ceiling,
                       ".bw")
    file.copy(f, file.path(target_dir, new_name), overwrite = TRUE)
  }

  message("BigWigs moved to: ", target_dir)
  invisible(bw_files)
}

exportGroupBW(ArchRProj = IsletProject,
              groupBy = "labels",
              normMethod = "ReadsInTSS",
              tile_size = 100,
              max_cells = 1000,
              ceiling = 4,
              in_dir = file.path(data_dir, "scATACseq/ArchR/IsletProject/GroupBigWigs/labels"),
              out_dir = file.path(data_dir, "scATACseq/ArchR/bigwig/")
)

exportGroupBW(ArchRProj = IsletProject,
              groupBy = "labels",
              normMethod = "ReadsInTSS",
              tile_size = 200,
              max_cells = 1000,
              ceiling = 4,
              in_dir = file.path(data_dir, "scATACseq/ArchR/IsletProject/GroupBigWigs/labels"),
              out_dir = file.path(data_dir, "scATACseq/ArchR/bigwig/")
)
