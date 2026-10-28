# ==============================================================================
# maping_clusters_to_labels.R
# Purpose : Map each HarmonyClusterR20 cluster to the cell type label held by most of its
#           cells (majority vote), and report how many cells match the mapped label.
# Inputs  : data/public/scATACseq/ArchR/IsletProject (with labels from 051)
# Outputs : data/project/HarmonyClusterR20_to_labels.csv; eval_R20 and overall_R20 (printed)
# ==============================================================================

data_dir      <- normalizePath("data/public", mustWork = FALSE)     # public/processed data
project_dir   <- normalizePath("data/project", mustWork = FALSE)    # small project files

suppressPackageStartupMessages({
  library(ArchR)
  library(dplyr)
  library(readr)
})

proj <- loadArchRProject(file.path(data_dir, "scATACseq/ArchR/IsletProject/"))

# Pull needed columns
ccd <- as.data.frame(ArchR::getCellColData(
  proj,
  select = c("HarmonyClusterR20", "labels")
)) %>%
  dplyr::mutate(
    HarmonyClusterR20 = factor(HarmonyClusterR20),
    labels = factor(labels)
  )

# Majority-vote label per cluster (best match)
majority_map <- function(df, cluster_col, label_col = "labels") {
  cluster_col <- rlang::ensym(cluster_col)
  label_col   <- rlang::ensym(label_col)

  df %>%
    dplyr::count(!!cluster_col, !!label_col, name = "n") %>%
    dplyr::group_by(!!cluster_col) %>%
    dplyr::slice_max(order_by = n, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::rename(cluster = !!cluster_col, mapped_label = !!label_col) %>%
    dplyr::select(cluster, mapped_label)
}

# R20 -> labels
map_R20 <- majority_map(ccd, HarmonyClusterR20, labels)
readr::write_csv(map_R20, file.path(project_dir, "HarmonyClusterR20_to_labels.csv"))

# Evaluate R20 mapping: for each HarmonyClusterR20, how many cells match mapped_label
eval_R20 <- ccd %>%
  dplyr::left_join(map_R20, by = c("HarmonyClusterR20" = "cluster")) %>%
  dplyr::group_by(HarmonyClusterR20, mapped_label) %>%
  dplyr::summarise(
    n_total   = dplyr::n(),
    n_correct = sum(labels == mapped_label, na.rm = TRUE),
    prop_correct = n_correct / n_total,
    .groups = "drop"
  ) %>%
  dplyr::arrange(dplyr::desc(prop_correct), dplyr::desc(n_total))

eval_R20

overall_R20 <- eval_R20 %>%
  dplyr::summarise(
    n_total   = sum(n_total),
    n_correct = sum(n_correct),
    prop_correct = n_correct / n_total
  )
overall_R20
