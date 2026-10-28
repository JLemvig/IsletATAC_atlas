# ==============================================================================
# cycling_clusters.R
# Purpose : Group the cycling R20 clusters (cycling_labels) of the trajectory project
#           into cycling, alpha_cycling, beta_cycling, delta_cycling, pp_cycling and
#           diffuse, chosen from gene score markers. Stored as cycling_manuel2 (used by
#           run_slingshot.R) and cycling_clusters (used by the figures).
# Inputs  : data/public/scATACseq/ArchR/cycling_plus_5k_abdp (with cycling_labels)
# Outputs : cycling_plus_5k_abdp updated in place (cycling_manuel2, cycling_clusters)
# ==============================================================================

data_dir      <- normalizePath("data/public", mustWork = FALSE)     # public/processed data

library(ArchR)

proj <- loadArchRProject(path = file.path(data_dir, "scATACseq/ArchR/cycling_plus_5k_abdp/"))

# cycling_labels ids (HarmonyClusterR20 of the final project)
alpha_cyc <- c("C6","C5","C4","C76","C15","C79","C174","C177","C178","C81","C50",
               "C175","C176","C180","C181","C182","C49","C12","C13","C9","C43","C69","C75","C62","C7")

beta_cyc  <- c("C167","C203","C110","C115","C116","C128","C202","C170","C114","C205",
               "C10","C2","C206","C3")

delta_cyc <- c("C11","C94","C95","C97","C102")

pp_cyc    <- c("C19","C8")

cycling   <- c("C1")

diffuse   <- c("C42","C16","Combined","C255","C256")

# map old -> new
map <- c(
  setNames(rep("alpha_cycling", length(alpha_cyc)), alpha_cyc),
  setNames(rep("beta_cycling",  length(beta_cyc)),  beta_cyc),
  setNames(rep("delta_cycling", length(delta_cyc)), delta_cyc),
  setNames(rep("pp_cycling",    length(pp_cyc)),    pp_cyc),
  setNames(rep("cycling",       length(cycling)),   cycling),
  setNames(rep("diffuse",       length(diffuse)),   diffuse)
)

# alpha, beta, delta and pp cells keep their label
old <- as.character(proj$cycling_labels)
new <- ifelse(old %in% names(map), unname(map[old]), old)

proj$cycling_manuel2 <- new
proj$cycling_clusters <- new

saveArchRProject(proj)
