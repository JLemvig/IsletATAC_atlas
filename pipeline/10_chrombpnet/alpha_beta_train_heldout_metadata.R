# ==============================================================================
# alpha_beta_train_heldout_metadata.R
# Purpose : Split the alpha and beta cells into a training set and a held-out set
#           (pre-diabetic donors and the Weng study are held out). Donors with fewer
#           than 200 cells of the cell type are excluded. The cell lists with their
#           disease status are used to make the ND and T2D pseudo-bulks for the
#           fine-tuning of ChromBPNet.
# Inputs  : data/public/scATACseq/ArchR/IsletProject (final ArchR project with labels
#           and donor metadata)
# Outputs : data/project/chrombpnet/<celltype>_{train,heldout}_metadata.csv
#           (celltype = alpha, beta; one row per cell, column Cell = cell barcode,
#           all cellColData columns, e.g. disease_harmonized)
# ==============================================================================

data_dir      <- normalizePath("data/public", mustWork = FALSE)     # public/processed data
project_dir   <- normalizePath("data/project", mustWork = FALSE)    # project outputs

library(ArchR)

dir.create(file.path(project_dir, "chrombpnet"), showWarnings = FALSE, recursive = TRUE)

## Load data
proj <- loadArchRProject(file.path(data_dir, "scATACseq/ArchR/IsletProject"))
all_metadata <- as.data.frame(proj@cellColData)

for (celltype in c("beta", "alpha")) {
  ## Subset to the cell type
  metadata <- all_metadata[ all_metadata$labels == celltype,]

  ## Define donors to use
  donor_count <- as.matrix(table(metadata$Sample))
  donors_keep <- rownames(donor_count)[which(donor_count >= 200)]

  ## Cells of all donors to keep, in donor order
  cells <- unlist(lapply(donors_keep, function(donor) rownames(metadata[ metadata$Sample %in% donor,])))

  ## Split the data
  metadata <- metadata[ metadata$Sample %in% donors_keep, ]
  metadata <- metadata[ match(cells, rownames(metadata)),]
  cells_heldout <- rownames(metadata[ metadata$disease_harmonized_hba1c == "PRE-T2D" | metadata$study == "Weng",])

  ## Setup metadata
  metadata_heldout <- metadata[ rownames(metadata) %in% cells_heldout,]
  metadata_heldout$Cell <- rownames(metadata_heldout)
  metadata_train <- metadata[ !(rownames(metadata) %in% cells_heldout),]
  metadata_train$Cell <- rownames(metadata_train)

  ## Save outputs
  write.csv(metadata_heldout, file.path(project_dir, "chrombpnet", paste0(celltype, "_heldout_metadata.csv")), row.names = FALSE)
  write.csv(metadata_train, file.path(project_dir, "chrombpnet", paste0(celltype, "_train_metadata.csv")), row.names = FALSE)
}
