# ==============================================================================
# peak_pc_association.R
# Purpose : Per cell type, a per-sample pseudobulk peak matrix and a DESeq2 likelihood
#           ratio test of each donor PC (adjusted for PC1-3, PC5-8 and PC18) against
#           peak accessibility.
# Inputs  : data/public/scATACseq/ArchR/IsletProject (with samplelabels from
#           04_annotation/054); data/public/scATACseq/correlations/pca.csv (donor PCs,
#           06_donor_model/donor_pca_covariate_tests.R)
# Outputs : per cell type in data/public/scATACseq/correlations/<celltype>/:
#           peakset_<celltype>.rds, se_peakmatrix_<celltype>.rds,
#           deseq_lrt_correlation_<celltype>.rds; for alpha, beta, delta and pp also
#           peakset_<celltype>_peaks.csv (peak set table)
# ==============================================================================

data_dir      <- normalizePath("data/public", mustWork = FALSE)     # public/processed data

library(BSgenome.Hsapiens.UCSC.hg38)
library(ArchR)
library(DESeq2)
addArchRThreads(threads = round(parallelly::availableCores()/4))

# Per cell type: maximum and minimum number of cells per sample, and the minimum
# fraction of cells in which a peak must be open
celltype_params <- data.frame(
  cell_type = c("beta", "alpha", "acinar", "ductal",
    "delta", "pp", "cycling", "activated_stellate", "quiescent_stellate",
    "endothelial", "myeloid", "lymphoid", "mast", "epsilon", "schwann"),
  max_cells = c(2000,2000,800,600,500,200,400,200,200,150,80,40,40,50,40),
  min_cells = c(100,100,50,50,50,25,25,10,10,5,0,0,0,0,0),
  peak_cutoff = c(0.01,0.01,0.02,0.02,0.02,0.05,0.05,0.05,0.05,0.05,0.10,0.10,
    0.10,0.10,0.10),
  stringsAsFactors = FALSE
)
pca <- read.csv(file.path(data_dir, "scATACseq/correlations/pca.csv"), row.names = 1)
IsletProject <- loadArchRProject(path = file.path(data_dir, "scATACseq/ArchR/IsletProject/"))

## Pseudobulk peak matrix per cell type
for(i in seq_along(celltype_params$cell_type)) {
  celltype <- celltype_params$cell_type[i]
  path <- paste0(file.path(data_dir, "scATACseq/correlations/"),celltype,"/")
  dir.create(path, recursive = TRUE, showWarnings = FALSE)

  # Samples with fewer than min_cells cells are dropped; samples with more than
  # max_cells cells are downsampled to max_cells
  set.seed(1)
  ccd <- IsletProject@cellColData
  ccd <- ccd[ ccd$labels == celltype,]
  barcodes_keep <- c()
  for (sample_name in unique(ccd$Sample)) {
    if (sum(ccd$Sample == sample_name) <= celltype_params[i,2]) {
      if (sum(ccd$Sample == sample_name) >= celltype_params[i,3]) {
        barcodes_keep = c(barcodes_keep, rownames(ccd[ ccd$Sample == sample_name,]))
      }
    } else {
      barcodes_keep = c(barcodes_keep, sample(rownames(ccd[ ccd$Sample == sample_name,]), celltype_params[i,2]))
    }
  }
  ccd <- ccd[ rownames(ccd) %in% barcodes_keep,]

  #Subset project
  subset_project <- subsetCells(IsletProject, cellNames = barcodes_keep)

  # Keep peaks open in at least peak_cutoff of the cells
  pm <- getMatrixFromProject(ArchRProj = subset_project,
                             useMatrix = "PeakMatrix",
                             verbose = FALSE,
                             binarize = FALSE)
  peakmatrix <- pm@assays@data@listData$PeakMatrix
  open_peak_fraction <- rowSums(peakmatrix > 0) / ncol(peakmatrix)
  # Summed counts per sample (columns <Sample>_<celltype>)
  se <- getGroupSE(ArchRProj = subset_project,
                   useMatrix = "PeakMatrix",
                   groupBy = "samplelabels",
                   divideN = FALSE,
                   threads = getArchRThreads(),
                   verbose = TRUE,
                   logFile = createLogFile("getGroupSE")
  )
  se_peakmatrix <- se@assays@data@listData$PeakMatrix
  se_peakmatrix_top <- se_peakmatrix[open_peak_fraction >= celltype_params[i,4], ]
  peakset <- getPeakSet(ArchRProj = subset_project)
  peakset$rownames <- rownames(se_peakmatrix)
  peakset$fraction_open <- open_peak_fraction
  peakset_top <- peakset[open_peak_fraction >= celltype_params[i,4], ]
  saveRDS(peakset_top, paste0(path,"peakset_",celltype,".rds"))
  saveRDS(se_peakmatrix_top,paste0(path,"se_peakmatrix_",celltype,".rds"))
}

align_pca_and_peakmatrix <- function(pca_df, peak_mat) {

  # remove last "_celltype" segment (works with multiple "_" earlier)
  pm_ids <- sub("_[^_]+$", "", colnames(peak_mat))

  # identify overlapping samples (cell types with "_" in the name need a second strip)
  common_ids <- sort(intersect(rownames(pca_df), pm_ids))
  if (length(common_ids) == 0L){
    pm_ids <- sub("_[^_]+$", "", pm_ids)
  }
  common_ids <- sort(intersect(rownames(pca_df), pm_ids))
  if (length(common_ids) == 0L)
    stop("No overlapping sample IDs between PCA and peakmatrix.")

  # subset PCA
  pca_sub <- pca_df[common_ids, , drop = FALSE]

  # subset peak matrix columns
  keep_cols <- pm_ids %in% common_ids
  peak_mat_sub <- peak_mat[, keep_cols, drop = FALSE]
  pm_ids_sub <- pm_ids[keep_cols]

  # reorder peak matrix to match PCA sample order
  ord <- match(common_ids, pm_ids_sub)
  peak_mat_aligned <- peak_mat_sub[, ord, drop = FALSE]

  # replace peak matrix column names with trimmed IDs
  colnames(peak_mat_aligned) <- pm_ids_sub[ord]

  # z-scale PCA after subsetting
  pca_scaled <- as.data.frame(scale(pca_sub))
  rownames(pca_scaled) <- rownames(pca_sub)
  colnames(pca_scaled) <- colnames(pca_sub)

  list(
    pca        = pca_scaled,
    peakmatrix = peak_mat_aligned
  )
}

## DESeq2 likelihood ratio test per donor PC
for(i in seq_along(celltype_params$cell_type)) {
  celltype <- celltype_params$cell_type[i]
  message("[", Sys.time(), "] ### Starting cell type: ", celltype)
  path <- paste0(file.path(data_dir, "scATACseq/correlations/"),celltype,"/")
  peakmatrix <- readRDS(paste0(path,"se_peakmatrix_",celltype,".rds"))
  peakset <- readRDS(paste0(path,"peakset_",celltype,".rds"))
  aligned <- align_pca_and_peakmatrix(pca,peakmatrix)
  pca_trimmed <- aligned$pca
  peakmatrix_aligned <- aligned$peakmatrix
  # Drop peaks with zero counts in more than 80% of the samples
  low_coverage_peaks <- rowSums(peakmatrix_aligned == 0) > round(ncol(peakmatrix_aligned)/100*80,digits = 0)
  peakmatrix_aligned_trimmed <- peakmatrix_aligned[!low_coverage_peaks,]
  # PCs to adjust for
  adjust_pcs <- c("PC1","PC2","PC3","PC5","PC6","PC7","PC8","PC18")
  # PCs to test
  test_pcs <- c("PC4","PC9","PC10","PC11","PC12","PC13",
                "PC14","PC15","PC16","PC17","PC19","PC20")
  lrt_results <- list()
  for (pc in test_pcs) {
    message("[", Sys.time(), "] Testing ", pc, " ...")
    full_terms    <- c(adjust_pcs, pc)
    full_formula  <- as.formula(paste("~", paste(full_terms, collapse = " + ")))
    reduced_formula <- as.formula(paste("~", paste(adjust_pcs, collapse = " + ")))

    # build DESeqDataSet for this PC
    dds_pc <- DESeqDataSetFromMatrix(
      countData = peakmatrix_aligned_trimmed,
      colData   = pca_trimmed,
      design    = full_formula
    )
    dds_pc <- estimateSizeFactors(dds_pc, type = "poscounts")
    dds_pc <- estimateDispersions(dds_pc)

    # LRT: test contribution of this PC
    dds_pc <- nbinomLRT(dds_pc, full = full_formula, reduced = reduced_formula)

    # results + coefficient for this PC
    res_pc <- as.data.frame(results(dds_pc, name = pc))
    coef_table <- coef(dds_pc)
    res_pc$coef <- coef_table[, pc]
    res_pc$feature <- rownames(res_pc)

    lrt_results[[pc]] <- res_pc
  }

  saveRDS(lrt_results,paste0(path,"deseq_lrt_correlation_",celltype,".rds"))
  message("[", Sys.time(), "] ### Finished cell type: ", celltype, "\n")
}

## Peak set tables (start is the 1-based GRanges start; PeakID is the peak matrix row name)
for (celltype in c("alpha", "beta", "delta", "pp")) {
  path <- paste0(file.path(data_dir, "scATACseq/correlations/"),celltype,"/")
  peakset <- as.data.frame(readRDS(paste0(path,"peakset_",celltype,".rds")))
  peakset_table <- data.frame(chromosome = peakset$seqnames,
                              start = peakset$start,
                              end = peakset$end,
                              nearestTSSDistance = peakset$nearestGene,
                              nearestTSS = peakset$distToTSS,
                              peakType = peakset$peakType,
                              GC = peakset$GC,
                              fraction = peakset$fraction_open,
                              PeakID = peakset$rownames)
  write.csv(peakset_table, paste0(path,"peakset_",celltype,"_peaks.csv"), quote = FALSE, row.names = FALSE)
}
