# ==============================================================================
# sex_validation2.R
# Purpose : Donor-level sex check: chrX/autosome fragment ratio per donor (TileMatrix),
#           two-group k-means, and the group with the higher XIST gene score called
#           female. The calls are printed for comparison with the donor metadata.
# Inputs  : data/public/scATACseq/ArchR/IsletProject (with donor_harmonized from 054)
# Outputs : sex_df (printed)
# ==============================================================================

data_dir      <- normalizePath("data/public", mustWork = FALSE)     # public/processed data

suppressPackageStartupMessages({
  library(ArchR)
  library(Matrix)
  library(SummarizedExperiment)
  library(GenomeInfoDb)
  library(dplyr)
})

proj <- loadArchRProject(file.path(data_dir, "scATACseq/ArchR/IsletProject/"))

stopifnot("donor_harmonized" %in% colnames(as.data.frame(ArchR::getCellColData(proj))))

# Aggregate the TileMatrix to tiles x donors (avoids loading the full tiles x cells matrix)
gse <- ArchR::getGroupSE(
  ArchRProj = proj,
  useMatrix = "TileMatrix",
  groupBy   = "donor_harmonized",
  threads   = 8
)

tile_mat <- SummarizedExperiment::assay(gse, "TileMatrix")
chr <- as.character(SummarizedExperiment::rowData(gse)$seqnames)

is_autosome <- chr %in% paste0("chr", 1:22)
is_X <- chr %in% c("chrX", "X")

aut_frags <- Matrix::colSums(tile_mat[is_autosome, , drop = FALSE])
x_frags   <- Matrix::colSums(tile_mat[is_X, , drop = FALSE])

tile_df <- tibble::tibble(
  donor_harmonized = colnames(tile_mat),
  autosome_frags = as.numeric(aut_frags),
  chrX_frags = as.numeric(x_frags),
  x_over_auto = chrX_frags / pmax(autosome_frags, 1)
)

# Cluster donors by X dosage (2 groups)
set.seed(1)
kmx <- stats::kmeans(scale(log1p(tile_df$x_over_auto)), centers = 2, nstart = 50)
tile_df$cluster_x <- kmx$cluster

# --- XIST pseudobulk (donor-level), used to label clusters ---
if (!"GeneScoreMatrix" %in% ArchR::getAvailableMatrices(proj)) {
  proj <- ArchR::addGeneScoreMatrix(proj, force = FALSE)
}

gse_gs <- ArchR::getGroupSE(
  ArchRProj = proj,
  useMatrix = "GeneScoreMatrix",
  groupBy   = "donor_harmonized",
  threads   = 8
)

# GeneScoreMatrix assay name can vary; take the first if needed
gs_assays <- names(SummarizedExperiment::assays(gse_gs))
gs_mat <- if ("GeneScoreMatrix" %in% gs_assays) {
  SummarizedExperiment::assay(gse_gs, "GeneScoreMatrix")
} else {
  SummarizedExperiment::assay(gse_gs, gs_assays[1])
}
rownames(gs_mat) <- gse_gs@elementMetadata$name

if (!"XIST" %in% rownames(gs_mat)) stop("XIST not found in GeneScoreMatrix rownames().")

xist <- as.numeric(gs_mat["XIST", ])

sex_df <- tile_df %>%
  dplyr::left_join(
    tibble::tibble(donor_harmonized = colnames(gs_mat), xist_score = xist),
    by = "donor_harmonized"
  )

# Label which X-dosage cluster is female by higher XIST
xist_by_cluster <- tapply(sex_df$xist_score, sex_df$cluster_x, median, na.rm = TRUE)
female_cluster <- as.integer(names(which.max(xist_by_cluster)))

sex_df <- sex_df %>%
  dplyr::mutate(
    sex_call = dplyr::if_else(cluster_x == female_cluster, "female", "male")
  ) %>%
  dplyr::select(donor_harmonized, sex_call, xist_score, x_over_auto, chrX_frags, autosome_frags) %>%
  dplyr::arrange(sex_call, dplyr::desc(xist_score))

# Sex call per donor
sex_df
