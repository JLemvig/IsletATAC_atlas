# ==============================================================================
# 013_IsletEpigenetics_FilterQCandExclusions.R
# Purpose : Remove low-quality cells using the per-sample cutoffs in cutoffvalues.csv
#           (log10 fragments, TSS enrichment, blacklist ratio, nucleosome ratio).
# Inputs  : data/public/scATACseq/ArchR/IsletProject1;
#           pipeline/02_archr_qc/cutoffvalues.csv (columns documented in 012)
# Outputs : data/public/scATACseq/ArchR/IsletProject2
# ==============================================================================

data_dir      <- normalizePath("data/public", mustWork = FALSE)     # public/processed data

## Load dependencies
library(ArchR)

## Set ArchR threads and reference genome
addArchRThreads(threads = round(parallelly::availableCores()/4))
addArchRGenome(genome = "hg38")

# Load ArchRProject
IsletProject <- loadArchRProject(path = file.path(data_dir, "scATACseq/ArchR/IsletProject1"))

## Removing Cells with bad QC
# Loading QC cutoff values; every sample needs a row in the csv
cutoffvalues <- read.csv("pipeline/02_archr_qc/cutoffvalues.csv")
stopifnot(all(unique(IsletProject$Sample) %in% cutoffvalues$Name))
# QC metrics per cell, computed by createArrowFiles in 011
df <- getCellColData(IsletProject,
                     select = c("Sample", "log10(nFrags)", "TSSEnrichment", "BlacklistRatio", "NucleosomeRatio"))

# Identifying low quality cells: a cell passes (1) if all metrics are within the
# cutoffs of its sample (strict comparisons), otherwise 0. df rows follow the
# cell order of the project.
cv <- cutoffvalues[match(df$Sample, cutoffvalues$Name), ]
PassQC <- ifelse(df$`log10(nFrags)` > cv$log10 &
                   df$TSSEnrichment > cv$Min.TSS.Enrichment &
                   df$TSSEnrichment < cv$Max.TSS.Enrichment &
                   df$BlacklistRatio < cv$Blacklist.Ratio &
                   df$NucleosomeRatio < cv$Nucleosome.ratio, 1, 0)
IsletProject$PassQC <- PassQC

# Subsetting project (overwrites IsletProject2 on rerun)
IsletProject2 <- subsetArchRProject(
  ArchRProj = IsletProject,
  cells = IsletProject$cellNames[IsletProject$PassQC == 1],
  outputDirectory = file.path(data_dir, "scATACseq/ArchR/IsletProject2"),
  dropCells = TRUE,
  force = TRUE
)
