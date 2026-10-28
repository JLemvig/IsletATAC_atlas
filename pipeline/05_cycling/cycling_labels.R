# ==============================================================================
# cycling_labels.R
# Purpose : Add cycling_labels to the final project: the cell type label for all
#           cells, except cycling cells, which get their HarmonyClusterR20 cluster id.
#           R20 ids with fewer than 30 cells are combined into "Combined".
# Inputs  : data/public/scATACseq/ArchR/IsletProject (with labels and HarmonyClusterR20
#           from 04_annotation/051)
# Outputs : IsletProject updated in place (cycling_labels)
# ==============================================================================

data_dir      <- normalizePath("data/public", mustWork = FALSE)     # public/processed data

library(ArchR)

IsletProject <- loadArchRProject(path = file.path(data_dir, "scATACseq/ArchR/IsletProject"))

# Cycling cells are split by their resolution 20 cluster
IsletProject$cycling_labels <- IsletProject$labels
is_cycling <- IsletProject$labels == "cycling"
IsletProject$cycling_labels[is_cycling] <- IsletProject$HarmonyClusterR20[is_cycling]

# Clusters with fewer than 30 cells are combined
n_cells <- table(IsletProject$cycling_labels)
IsletProject$cycling_labels[IsletProject$cycling_labels %in% names(n_cells)[n_cells < 30]] <- "Combined"

saveArchRProject(IsletProject)
