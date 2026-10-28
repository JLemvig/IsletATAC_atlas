# ==============================================================================
# 012_IsletEpigenetics_CreateQCplots.R
# Purpose : Per-sample QC plots of TSS enrichment against log10(nFrags), blacklist
#           ratio and nucleosome ratio, with the cutoffs from cutoffvalues.csv
#           overlaid. The plots were used to choose the per-sample cutoffs.
# Inputs  : data/public/scATACseq/ArchR/IsletProject1;
#           pipeline/02_archr_qc/cutoffvalues.csv
# Outputs : data/public/scATACseq/ArchR/IsletProject1/QC/<sample>.pdf
#
# cutoffvalues.csv, one row per sample. Columns:
#   Name                     ArchR sample name
#   Unique fragments cutoff  minimum unique fragments (the code uses the log10 column)
#   log10                    log10 of the previous column; cells need log10(nFrags) > log10
#   Min TSS Enrichment       cells need TSSEnrichment > this
#   Max TSS Enrichment       cells need TSSEnrichment < this
#   Blacklist Ratio          cells need BlacklistRatio < this
#   Nucleosome ratio         cells need NucleosomeRatio < this
#   FRIP                     cells need FRIP > this (used in 022 and 023)
# ==============================================================================

data_dir      <- normalizePath("data/public", mustWork = FALSE)     # public/processed data

## Load dependencies
library(ArchR)
library(ggplot2)
library(cowplot)

## Set ArchR threads and reference genome
addArchRThreads(threads = round(parallelly::availableCores()/4))
addArchRGenome(genome = "hg38")

## Load ArchRProject
IsletProject <- loadArchRProject(path = file.path(data_dir, "scATACseq/ArchR/IsletProject1"))

# Read cut off values from .csv
cutoffvalues <- read.csv("pipeline/02_archr_qc/cutoffvalues.csv")

qc_dir <- file.path(data_dir, "scATACseq/ArchR/IsletProject1/QC")
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

# QC metrics per cell, computed by createArrowFiles in 011
df <- getCellColData(IsletProject,
                     select = c("Sample", "log10(nFrags)", "TSSEnrichment", "BlacklistRatio", "NucleosomeRatio"))
# One page per sample; the red numbers are the number of cells in each region
# delimited by the dashed cutoff lines
invisible(lapply(unique(IsletProject$Sample),function(i) {
  # Read cut off values, set default values if not specified in .csv
  if(i %in% cutoffvalues[,1]){
    cv <- cutoffvalues[cutoffvalues[,1] %in% i, ]
  } else {
    cv <- data.frame("Nucleosome.ratio" = 6, "Blacklist.Ratio" = 0.02, "Min.TSS.Enrichment" = 4,
                     "Max.TSS.Enrichment" = 30, "log10" = 3)
  }
  dfi <- df[df[,1] %in% i, ]
  # TSS / nFrags plot
  p1 <- ggPoint(
    x = dfi$`log10(nFrags)`,
    y = dfi$TSSEnrichment,
    colorDensity = TRUE,
    continuousSet = "sambaNight",
    xlabel = "log10(nFrags)",
    ylabel = "TSS Enrichment",
    xlim = c(log10(500), max(dfi$`log10(nFrags)`)),
    ylim = c(0, max(dfi$TSSEnrichment))
  ) +
    geom_hline(yintercept = cv$Min.TSS.Enrichment, lty = "dashed") +
    geom_vline(xintercept = cv$log10, lty = "dashed") +
    geom_hline(yintercept = cv$Max.TSS.Enrichment, lty = "dashed") +
    annotate("text", x = cv$log10 + 0.2, y = 2, col = "red",
             label = sum(dfi$TSSEnrichment < cv$Min.TSS.Enrichment & dfi$`log10(nFrags)` > cv$log10)) +
    annotate("text", x = 2.9, y = 2, col = "red",
             label = sum(dfi$TSSEnrichment < cv$Min.TSS.Enrichment & dfi$`log10(nFrags)` < cv$log10)) +
    annotate("text", x = 2.9, y = mean(c(cv$Max.TSS.Enrichment,cv$Min.TSS.Enrichment)), col = "red",
             label = sum(dfi$`log10(nFrags)` < cv$log10 &
                           dfi$TSSEnrichment > cv$Min.TSS.Enrichment &
                           dfi$TSSEnrichment < cv$Max.TSS.Enrichment)) +
    annotate("text", x = 2.9, y = cv$Max.TSS.Enrichment + 2, col = "red",
             label = sum(dfi$`log10(nFrags)` < cv$log10 & dfi$TSSEnrichment > cv$Max.TSS.Enrichment)) +
    annotate("text", x = cv$log10 + 0.2, y = cv$Max.TSS.Enrichment + 2, col = "red",
             label = sum(dfi$`log10(nFrags)` > cv$log10 & dfi$TSSEnrichment > cv$Max.TSS.Enrichment)) +
    annotate("text", x = cv$log10 + 0.2, y = mean(c(cv$Max.TSS.Enrichment,cv$Min.TSS.Enrichment)), col = "red",
             label = sum(dfi$`log10(nFrags)` > cv$log10 &
                           dfi$TSSEnrichment > cv$Min.TSS.Enrichment &
                           dfi$TSSEnrichment < cv$Max.TSS.Enrichment))
  #TSS / Blacklist ratio plot
  p2 <- ggPoint(
    x = dfi$BlacklistRatio,
    y = dfi$TSSEnrichment,
    colorDensity = TRUE,
    continuousSet = "sambaNight",
    xlabel = "BlacklistRatio",
    ylabel = "TSS Enrichment",
    xlim = c(0, cv$Blacklist.Ratio + 0.005),
    ylim = c(0, max(dfi$TSSEnrichment))
  ) +
    geom_hline(yintercept = cv$Min.TSS.Enrichment, lty = "dashed") +
    geom_vline(xintercept = cv$Blacklist.Ratio, lty = "dashed") +
    geom_hline(yintercept = cv$Max.TSS.Enrichment, lty = "dashed") +
    annotate("text", x = cv$Blacklist.Ratio + 0.001, y = 2, col = "red",
             label = sum(dfi$BlacklistRatio > cv$Blacklist.Ratio & dfi$TSSEnrichment < cv$Min.TSS.Enrichment)) +
    annotate("text", x = cv$Blacklist.Ratio + 0.001, y = cv$Max.TSS.Enrichment + 2, col = "red",
             label = sum(dfi$BlacklistRatio > cv$Blacklist.Ratio & dfi$TSSEnrichment > cv$Max.TSS.Enrichment)) +
    annotate("text", x = cv$Blacklist.Ratio + 0.001, y = mean(c(cv$Max.TSS.Enrichment,cv$Min.TSS.Enrichment)),
             col = "red",
             label = sum(dfi$BlacklistRatio > cv$Blacklist.Ratio &
                           dfi$TSSEnrichment < cv$Max.TSS.Enrichment &
                           dfi$TSSEnrichment > cv$Min.TSS.Enrichment))
  #TSS / Nucleosome ratio plot
  p3 <- ggPoint(
    x = dfi$NucleosomeRatio,
    y = dfi$TSSEnrichment,
    colorDensity = TRUE,
    continuousSet = "sambaNight",
    xlabel = "NucleosomeRatio",
    ylabel = "TSS Enrichment",
    xlim = c(0, cv$Nucleosome.ratio + 0.5),
    ylim = c(0, max(dfi$TSSEnrichment))
  ) +
    geom_hline(yintercept = cv$Min.TSS.Enrichment, lty = "dashed") +
    geom_vline(xintercept = cv$Nucleosome.ratio, lty = "dashed") +
    geom_hline(yintercept = cv$Max.TSS.Enrichment, lty = "dashed") +
    annotate("text", x = cv$Nucleosome.ratio + 0.2, y = 2, col = "red",
             label = sum(dfi$NucleosomeRatio > cv$Nucleosome.ratio & dfi$TSSEnrichment < cv$Min.TSS.Enrichment)) +
    annotate("text", x = cv$Nucleosome.ratio + 0.2, y = cv$Max.TSS.Enrichment + 2, col = "red",
             label = sum(dfi$NucleosomeRatio > cv$Nucleosome.ratio & dfi$TSSEnrichment > cv$Max.TSS.Enrichment)) +
    annotate("text", x = cv$Nucleosome.ratio + 0.2, y = mean(c(cv$Max.TSS.Enrichment,cv$Min.TSS.Enrichment)),
             col = "red",
             label = sum(dfi$NucleosomeRatio > cv$Nucleosome.ratio &
                           dfi$TSSEnrichment < cv$Max.TSS.Enrichment &
                           dfi$TSSEnrichment > cv$Min.TSS.Enrichment))
  # Title
  title <- cowplot::ggdraw() + cowplot::draw_label(paste0(i, ": ", nrow(dfi), " cells"), fontface = "bold")

  # Finalize layout
  p4 <- cowplot::plot_grid(p1, p2, p3, nrow = 1)
  p5 <- cowplot::plot_grid(title, p4, ncol = 1, rel_heights = c(0.1, 1))

  # Save plot
  file <- file.path(qc_dir, paste0(i, ".pdf"))
  cowplot::save_plot(file, p5, base_height = 8, base_width = 19)
}))
