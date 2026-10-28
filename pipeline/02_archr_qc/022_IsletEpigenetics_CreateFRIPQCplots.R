# ==============================================================================
# 022_IsletEpigenetics_CreateFRIPQCplots.R
# Purpose : Per-sample plots of TSS enrichment against FRIP with the FRIP cutoff
#           from cutoffvalues.csv overlaid. The plots were used to choose the
#           per-sample FRIP cutoffs.
# Inputs  : data/public/scATACseq/ArchR/IsletProject2 (with the PeakMatrix from 021);
#           pipeline/02_archr_qc/cutoffvalues.csv (columns documented in 012)
# Outputs : data/public/scATACseq/ArchR/IsletProject2/FRIP/<sample>.pdf
# ==============================================================================

data_dir      <- normalizePath("data/public", mustWork = FALSE)     # public/processed data

## Load dependencies
library(ArchR)
library(ggplot2)

## Set ArchR threads and reference genome
addArchRThreads(threads = round(parallelly::availableCores()/4))
addArchRGenome(genome = "hg38")

# Load ArchRProject
IsletProject2 <- loadArchRProject(path = file.path(data_dir, "scATACseq/ArchR/IsletProject2"))

# Read cut off values from .csv (only the Name and FRIP columns are used)
cutoffvalues <- read.csv("pipeline/02_archr_qc/cutoffvalues.csv")

frip_dir <- file.path(data_dir, "scATACseq/ArchR/IsletProject2/FRIP")
dir.create(frip_dir, recursive = TRUE, showWarnings = FALSE)

df <- getCellColData(IsletProject2, select = c("FRIP", "TSSEnrichment","Sample"))
# One page per sample; the red number is the number of cells below the FRIP cutoff
invisible(lapply(unique(df[,"Sample"]),function(i) {
  # Read cut off values, set default value if not specified in .csv
  if(i %in% cutoffvalues[,1]){
    fripcut <- cutoffvalues[cutoffvalues[,"Name"] %in% i,"FRIP"]
  } else {
    fripcut <- 0.5
  }
  dfi <- df[df[,"Sample"] %in% i, ]
  # TSS / FRIP plot
  var <- ggPoint(x = dfi$FRIP,
                 y = dfi$TSSEnrichment,
                 colorDensity = TRUE,
                 continuousSet = "sambaNight",
                 xlabel = "Fraction of Reads in Peak regions (FRIP)",
                 ylabel = "TSS Enrichment",
                 xlim = c(0, max(dfi$FRIP)+0.1),
                 ylim = c(0, max(dfi$TSSEnrichment)+3)) +
    ggtitle(paste0(i," nCells: ",nrow(dfi))) +
    geom_vline(xintercept = fripcut, lty = "dashed") +
    annotate("text", x = fripcut - 0.05, y = 2, col = "red",
             label = sum(dfi$FRIP < fripcut))

  # Save plot
  file <- file.path(frip_dir, paste0(i, ".pdf"))
  cowplot::save_plot(file, var, base_height = 8, base_width = 7)
}))
