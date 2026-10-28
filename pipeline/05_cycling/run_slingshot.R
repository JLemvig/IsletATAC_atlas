# ==============================================================================
# run_slingshot.R
# Purpose : Slingshot trajectories on Harmony from the cycling cluster to each endocrine
#           cell type (alpha, beta, delta, pp), via the type's cycling cluster.
# Inputs  : data/public/scATACseq/ArchR/cycling_plus_5k_abdp (with cycling_manuel2 from
#           cycling_clusters.R)
# Outputs : cycling_plus_5k_abdp updated in place (slingPseudotime_alpha, _beta, _delta,
#           _pp)
# ==============================================================================

data_dir      <- normalizePath("data/public", mustWork = FALSE)     # public/processed data

suppressPackageStartupMessages({
  library(ArchR)
  library(Matrix)
  library(SingleCellExperiment)
  library(S4Vectors)
  library(slingshot)
})

proj <- loadArchRProject(path = file.path(data_dir, "scATACseq/ArchR/cycling_plus_5k_abdp/"))

# Slingshot per endpoint, restricted to the start cluster, the endpoint's cycling cluster
# and the endpoint (cycling -> <type>_cycling -> <type>); diffuse cells are excluded
run_slingshot_supervised <- function(
    proj,
    group_col   = "cycling_manuel2",
    reduced_use = "Harmony",
    embed_use   = "HarmonyUMAP",
    start       = "cycling",
    endpoints   = c("alpha","beta","delta","pp"),
    exclude     = c("diffuse"),
    shared      = NULL,   # optional shared intermediates, e.g. "beta_cycling"
    stretch     = 0
){
  # --- reduced dims + embedding (aligned) ---
  rd <- getReducedDims(proj, reducedDims = reduced_use)

  um <- as.data.frame(getEmbedding(proj, embedding = embed_use))
  um <- um[match(rownames(rd), rownames(um)), , drop = FALSE]
  um_mat <- as.matrix(um[, 1:2, drop = FALSE])
  colnames(um_mat) <- c("UMAP1","UMAP2")

  # --- group labels (aligned) ---
  cl <- as.character(getCellColData(proj, select = group_col)[, 1])
  names(cl) <- proj$cellNames
  cl <- cl[match(rownames(rd), names(cl))]

  # base filtering
  keep_base <- !is.na(cl) & !(cl %in% exclude)

  # output containers
  out <- vector("list", length(endpoints))
  names(out) <- endpoints

  pt_mat <- matrix(
    NA_real_,
    nrow = length(proj$cellNames),
    ncol = length(endpoints),
    dimnames = list(proj$cellNames, endpoints)
  )

  for (ct in endpoints) {
    allowed <- unique(c(start, shared, paste0(ct, "_cycling"), ct))

    keep <- keep_base & (cl %in% allowed)

    rd2 <- rd[keep, , drop = FALSE]
    um2 <- um_mat[keep, , drop = FALSE]
    cl2 <- droplevels(factor(cl[keep]))

    # skip if start/end absent after filtering
    if (!(start %in% levels(cl2)) || !(ct %in% levels(cl2))) {
      out[[ct]] <- list(
        endpoint = ct,
        allowed  = allowed,
        error    = sprintf("Missing start (%s) or end (%s) after filtering.", start, ct)
      )
      next
    }

    # build SCE with proper cell names
    sce <- SingleCellExperiment(
      assays      = list(dummy = Matrix(0, nrow = 1, ncol = nrow(rd2), sparse = TRUE)),
      colData     = DataFrame(cluster = cl2),
      reducedDims = SimpleList(HARMONY = rd2, UMAP = um2)
    )
    colnames(sce) <- rownames(rd2)

    # run slingshot (single endpoint)
    sce <- slingshot(
      sce,
      clusterLabels = "cluster",
      reducedDim    = "HARMONY",
      start.clus    = start,
      end.clus      = ct,
      stretch       = stretch
    )

    # extract pseudotime (usually 1 lineage here)
    pt <- slingPseudotime(sce)[, 1]
    names(pt) <- colnames(sce)

    # store in matrix using ArchR cellNames
    common <- intersect(names(pt), rownames(pt_mat))
    pt_mat[common, ct] <- pt[common]

    out[[ct]] <- list(
      endpoint = ct,
      allowed  = allowed,
      lineage  = slingLineages(sce)[[1]],
      curves   = slingCurves(sce),
      pseudotime = pt,
      sce      = sce
    )
  }

  list(
    pseudotime_matrix = pt_mat,
    per_endpoint      = out
  )
}

set.seed(1)
res <- run_slingshot_supervised(
  proj       = proj,
  group_col  = "cycling_manuel2",
  start      = "cycling",
  endpoints  = c("alpha","beta","delta","pp"),
  exclude    = c("diffuse"),
  shared     = NULL,
  stretch    = 0
)

# Pseudotime per cell (rows) and endpoint (columns); added to the project as
# slingPseudotime_<endpoint>
pt_mat <- res$pseudotime_matrix

stopifnot(all(proj$cellNames %in% rownames(pt_mat)))

for (ct in colnames(pt_mat)) {
  nm  <- paste0("slingPseudotime_", ct)
  vec <- pt_mat[proj$cellNames, ct]
  names(vec) <- proj$cellNames

  proj <- addCellColData(
    ArchRProj = proj,
    data      = vec,
    cells     = proj$cellNames,
    name      = nm,
    force     = TRUE
  )
}

saveArchRProject(proj)
