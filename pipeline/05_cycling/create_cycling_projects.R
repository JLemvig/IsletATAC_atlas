# ==============================================================================
# create_cycling_projects.R
# Purpose : Create the trajectory project: all cycling cells plus up to 5,000 randomly
#           drawn alpha, beta, delta and pp cells each.
# Inputs  : data/public/scATACseq/ArchR/IsletProject (with cycling_labels from
#           cycling_labels.R)
# Outputs : data/public/scATACseq/ArchR/cycling_plus_5k_abdp
# ==============================================================================

data_dir      <- normalizePath("data/public", mustWork = FALSE)     # public/processed data

library(BSgenome.Hsapiens.UCSC.hg38)
library(ArchR)
library(tidyverse)
addArchRThreads(threads = round(parallelly::availableCores()/2))

IsletProject <- loadArchRProject(path = file.path(data_dir, "scATACseq/ArchR/IsletProject"))

cycling_cells <- IsletProject$cellNames[IsletProject$labels == "cycling"]

# Sample 5k from each endocrine label (or all available if < 5k)
set.seed(1)
sample_n_or_all <- function(label, n = 5000) {
  cells <- IsletProject$cellNames[IsletProject$labels == label]
  if (length(cells) <= n) {
    return(cells)
  } else {
    return(sample(cells, n))
  }
}

alpha_cells <- sample_n_or_all("alpha", 5000)
beta_cells  <- sample_n_or_all("beta",  5000)
delta_cells <- sample_n_or_all("delta", 5000)
pp_cells    <- sample_n_or_all("pp",    5000)

# Combine and de-duplicate
cells_keep <- unique(c(cycling_cells, alpha_cells, beta_cells, delta_cells, pp_cells))

IsletProject_cycling5k <- subsetArchRProject(
  ArchRProj = IsletProject,
  cells = cells_keep,
  outputDirectory = file.path(data_dir, "scATACseq/ArchR/cycling_plus_5k_abdp/"),
  dropCells = TRUE
)
