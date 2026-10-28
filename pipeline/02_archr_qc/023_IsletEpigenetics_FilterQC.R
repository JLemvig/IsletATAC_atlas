# ==============================================================================
# 023_IsletEpigenetics_FilterQC.R
# Purpose : Remove cells with FRIP below the per-sample cutoff in cutoffvalues.csv.
# Inputs  : data/public/scATACseq/ArchR/IsletProject2 (with the PeakMatrix from 021);
#           pipeline/02_archr_qc/cutoffvalues.csv (columns documented in 012)
# Outputs : data/public/scATACseq/ArchR/IsletProject3
# ==============================================================================

data_dir      <- normalizePath("data/public", mustWork = FALSE)     # public/processed data

## Load dependencies
library(ArchR)

## Set ArchR threads and reference genome
addArchRThreads(threads = round(parallelly::availableCores()/4))
addArchRGenome(genome = "hg38")

# Load ArchRProject
IsletProject2 <- loadArchRProject(path = file.path(data_dir, "scATACseq/ArchR/IsletProject2"))

## Removing Cells with bad QC
# Loading QC cutoff values (only the Name and FRIP columns are used); every sample
# needs a row in the csv
cutoffvalues <- read.csv("pipeline/02_archr_qc/cutoffvalues.csv")
stopifnot(all(unique(IsletProject2$Sample) %in% cutoffvalues$Name))

# Identify cells with low FRIP: a cell passes (1) if FRIP is above the cutoff of its
# sample, otherwise 0. Overwrites the PassQC column from 013. df rows follow the
# cell order of the project.
df <- getCellColData(IsletProject2, select = c("FRIP","Sample"))
PassQC <- ifelse(df$FRIP > cutoffvalues$FRIP[match(df$Sample, cutoffvalues$Name)], 1, 0)
IsletProject2$PassQC <- PassQC

# Subsetting project (overwrites IsletProject3 on rerun)
IsletProject3 <- subsetArchRProject(
  ArchRProj = IsletProject2,
  cells = IsletProject2$cellNames[IsletProject2$PassQC == 1],
  outputDirectory = file.path(data_dir, "scATACseq/ArchR/IsletProject3"),
  dropCells = TRUE,
  force = TRUE
)
