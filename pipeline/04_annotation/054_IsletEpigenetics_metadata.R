# ==============================================================================
# 054_IsletEpigenetics_metadata.R
# Purpose : Add the harmonised donor metadata (donor, disease, HbA1c, BMI, age, sex,
#           ethnicity, cell type fractions) to the cells of IsletProject.
# Inputs  : data/public/scATACseq/ArchR/IsletProject (with labels from 051);
#           data/public/scATACseq/Project/meta.csv (donor and sample metadata, one row
#           per sample, column ID; same content as ESM Table 1)
# Outputs : IsletProject updated in place (metadata columns and samplelabels)
# ==============================================================================

data_dir      <- normalizePath("data/public", mustWork = FALSE)     # public/processed data

library(tidyverse)
library(ArchR)

IsletProject <- loadArchRProject(file.path(data_dir, "scATACseq/ArchR/IsletProject/"))
meta <- read.csv(file.path(data_dir, "scATACseq/Project/meta.csv"))
df <- data.frame("ID" = IsletProject$Sample)
df <- df %>%
  left_join(meta, by = "ID")

IsletProject$donor_harmonized <- df$donor_harmonized
IsletProject$disease_harmonized <- df$harmonized_disease
IsletProject$disease_harmonized_hba1c <- df$harmonized_disease_hba1c
IsletProject$metformin_harmonized <- df$metformin_harmonized
IsletProject$hba1c <- df$hba1c
IsletProject$bmi <- df$BMI
IsletProject$ethnicity_harmonized <- df$ethnicity_harmonized
IsletProject$gender <- df$Gender
IsletProject$age <- df$Age
IsletProject$total_cells <- df$total_cells
IsletProject$alpha_fraction <- df$alpha_fraction
IsletProject$beta_fraction <- df$beta_fraction
IsletProject$cycling_fraction <- df$cycling_fraction
IsletProject$delta_fraction <- df$delta_fraction
IsletProject$pp_fraction <- df$pp_fraction
IsletProject$total_endocrine_cells <- df$total_endocrine_cells
IsletProject$endocrine_fraction <- df$fraction_endocrine
IsletProject$alpha_fraction_endocrine <- df$alpha_fraction_endocrine
IsletProject$beta_fraction_endocrine <- df$beta_fraction_endocrine
IsletProject$delta_fraction_endocrine <- df$delta_fraction_endocrine
IsletProject$pp_fraction_endocrine <- df$pp_fraction_endocrine
IsletProject$epsilon_fraction_endocrine <- df$epsilon_fraction_endocrine
IsletProject$cycling_fraction_endocrine <- df$cycling_fraction_endocrine
IsletProject$total_stellate <- df$total_stellate
IsletProject$activated_stellate_fraction <- df$activated_stellate_fraction

# Sample and cell type, used to make per-sample pseudobulks of one cell type
IsletProject$samplelabels <- paste(IsletProject$Sample, IsletProject$labels, sep = "_")

saveArchRProject(IsletProject)

