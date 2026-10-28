# ==============================================================================
# significant_peaks.R
# Purpose : Collect the significant peak-PC associations (adjusted p < 0.05) per cell
#           type, with their direction (sign of the PC coefficient).
# Inputs  : data/public/scATACseq/correlations/<celltype>/deseq_lrt_correlation_<celltype>.rds
#           (peak_pc_association.R)
# Outputs : data/public/scATACseq/correlations/<celltype>/deseq_lrt_joint_signif_<celltype>.rds
# ==============================================================================

data_dir      <- normalizePath("data/public", mustWork = FALSE)     # public/processed data

celltypes <- c("beta", "alpha", "acinar", "ductal",
               "delta", "pp", "cycling", "activated_stellate", "quiescent_stellate",
               "endothelial", "myeloid", "lymphoid", "mast", "epsilon", "schwann")

for (celltype in celltypes) {
  path <- paste0(file.path(data_dir, "scATACseq/correlations/"),celltype,"/")
  lrt_results <- readRDS(paste0(path,"deseq_lrt_correlation_",celltype,".rds"))
  joint_list <- lapply(names(lrt_results), function(pc) {
    df <- lrt_results[[pc]]

    # add PC label
    df$PC <- pc

    # keep only non-NA and significant
    df <- df[!is.na(df$padj) & df$padj < 0.05, , drop = FALSE]

    if (nrow(df) == 0L) return(NULL)

    # direction relative to coef
    df$direction <- ifelse(df$coef > 0, "positive", "negative")

    df
  })
  joint_res <- do.call(rbind, joint_list)
  saveRDS(joint_res,paste0(path, "deseq_lrt_joint_signif_", celltype, ".rds"))
}
