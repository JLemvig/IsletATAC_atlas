# ==============================================================================
# module_scores.R
# Purpose : Add gene score module scores (GeneScoreMatrix) for endocrine marker genes
#           and cell cycle genes: alpha, beta, delta and gamma (pp) markers and Seurat
#           S and G2M genes on the trajectory project, and the cell cycle score on the
#           final project.
# Inputs  : data/public/scATACseq/ArchR/cycling_plus_5k_abdp;
#           data/public/scATACseq/ArchR/IsletProject
# Outputs : both projects updated in place (ModuleScore.alpha_top, ModuleScore.beta_top,
#           ModuleScore.delta_top, ModuleScore.gamma_top, ModuleScore.seurat.cycling;
#           only ModuleScore.seurat.cycling on IsletProject)
# ==============================================================================

data_dir      <- normalizePath("data/public", mustWork = FALSE)     # public/processed data

suppressPackageStartupMessages({
  library(ArchR)
})

IsletProject_cycling5k <- loadArchRProject(path = file.path(data_dir, "scATACseq/ArchR/cycling_plus_5k_abdp/"))
IsletProject <- loadArchRProject(path = file.path(data_dir, "scATACseq/ArchR/IsletProject"))

# Top 10 marker genes per endocrine cell type
alpha_top <- c("GCG", "TTR", "IRX2", "HIGD1A", "GLS", "TM4SF4", "FAP", "GPX3", "SLC7A2", "GC")
beta_top  <- c("INS", "IAPP", "G6PC2", "ADCYAP1", "ERO1B", "DLK1", "NPTX2", "GSN", "INS-IGF2", "HADH")
gamma_top <- c("PPY", "ID2", "GCNT3", "FXYD2", "STMN2", "THSD7A", "SLITRK6", "SERTM1", "TM4SF4", "ETV1")
delta_top <- c("SST", "RBP4", "LEPR", "RGS2", "SEC11C", "PRG4", "BCHE", "ADGRL2", "HHEX", "SLC38A1")

# Seurat cell cycle genes (S and G2M), without s.genes[15] and g2m.genes[16] and [30]
seurat.s.genes   <- Seurat::cc.genes$s.genes[c(-15)]
seurat.g2m.genes <- Seurat::cc.genes$g2m.genes[c(-16,-30)]

modules <- list(
  alpha_top = alpha_top,
  beta_top  = beta_top,
  gamma_top = gamma_top,
  delta_top = delta_top,
  seurat.cycling = c(seurat.s.genes, seurat.g2m.genes)
)

# One module per call; each adds the column ModuleScore.<module name>
for (i in seq_along(modules)) {
  IsletProject_cycling5k <- addModuleScore(
    ArchRProj = IsletProject_cycling5k,
    useMatrix = "GeneScoreMatrix",
    name      = "ModuleScore",
    features  = modules[i],
    nBin      = 25,
    nBgd      = 100,
    seed      = 1
  )
}

IsletProject <- addModuleScore(
  ArchRProj = IsletProject,
  useMatrix = "GeneScoreMatrix",
  name      = "ModuleScore",
  features  = modules["seurat.cycling"],
  nBin      = 25,
  nBgd      = 100,
  seed      = 1
)

saveArchRProject(IsletProject_cycling5k)
saveArchRProject(IsletProject)
