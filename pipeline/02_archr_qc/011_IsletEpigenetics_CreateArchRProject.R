# ==============================================================================
# 011_IsletEpigenetics_CreateArchRProject.R
# Purpose : Create ArrowFiles from the sorted fragment files (per sample, minTSS 2,
#           minFrags 500) and build the initial ArchRProject (IsletProject1),
#           excluding samples from donors under 18 years.
# Inputs  : data/public/scATACseq/<study>/Preprocessed/<sample>/frags.sort.bed.gz;
#           data/sensitive/snATACseq/FUSION/Preprocessed/<sample>/frags.sort.bed.gz
# Outputs : <sample>.arrow in each sample folder;
#           data/public/scATACseq/ArchR/IsletProject1
# ==============================================================================

data_dir      <- normalizePath("data/public", mustWork = FALSE)     # public/processed data
sensitive_dir <- normalizePath("data/sensitive", mustWork = FALSE)  # controlled-access data (FUSION)

library(ArchR)

## Set ArchR threads and reference genome
addArchRThreads(threads = round(parallelly::availableCores()/4))
addArchRGenome(genome = "hg38")

## Create Arrowfiles from Fragment files.
# FUSION is controlled access and is read from data/sensitive
study <- c("HPAP","Chiou","FUSION","Wang","Benaglio","Thibodeau","Augsornworawat","Weng")
for(i in study){
  folder <- ifelse(i == "FUSION", file.path(sensitive_dir, "snATACseq/"), file.path(data_dir, "scATACseq/"))
  dir <- paste0(folder,i,"/Preprocessed/")
  Samples <- list.files(dir)
  for(j in Samples){
    # createArrowFiles writes <sample>.arrow to the working directory, so work inside
    # the sample folder and return afterwards
    old_wd <- getwd()
    setwd(paste0(dir,j))
    inputFiles <- paste0(dir,j,"/frags.sort.bed.gz")
    # minTSS and minFrags are permissive here; the per-sample cutoffs from
    # cutoffvalues.csv are applied in 013 and 023
    ArrowFiles <- createArrowFiles(
      inputFiles = inputFiles,
      sampleNames = j,
      minTSS = 2,
      minFrags = 500,
      addTileMat = TRUE,
      addGeneScoreMat = TRUE,
      force = TRUE   # overwrite existing ArrowFiles on rerun
    )
    setwd(old_wd)
  }
}

## Create ArchRproject from Arrowfiles
# Locating ArrowFiles
ArrowFiles <- c()
for(i in study){
  folder <- ifelse(i == "FUSION", file.path(sensitive_dir, "snATACseq/"), file.path(data_dir, "scATACseq/"))
  dir <- paste0(folder,i,"/Preprocessed/")
  Samples <- list.files(dir)
  for(j in Samples){
    ArrowFiles <- append(ArrowFiles, paste0(dir,j,"/",j,".arrow"))
  }
}

# Removing samples from donors under 18 years (age from the donor metadata)
samples_under_18 <- c("HPAP-039", "HPAP-042", "HPAP-044", "HPAP-047", "HPAP-084", "HPAP-104", "HPAP-107")
stopifnot(all(paste0(samples_under_18, ".arrow") %in% basename(ArrowFiles)))
ArrowFiles <- ArrowFiles[!basename(ArrowFiles) %in% paste0(samples_under_18, ".arrow")]

# Creating ArchRProject
IsletProject <- ArchRProject(
  ArrowFiles = ArrowFiles,
  outputDirectory = file.path(data_dir, "scATACseq/ArchR/IsletProject1"),
  showLogo = FALSE,
  copyArrows = TRUE
)

# Labeling samples by study (Wang et al. samples are named MM_* and JYH_*)
IsletProject$study <- dplyr::case_when(
  grepl("HPAP",IsletProject$Sample) ~ "HPAP",
  grepl("FUSION",IsletProject$Sample) ~ "FUSION",
  grepl("Chiou",IsletProject$Sample) ~ "Chiou",
  grepl("Benaglio",IsletProject$Sample) ~ "Benaglio",
  grepl("Thibodeau",IsletProject$Sample) ~ "Thibodeau",
  grepl("Augsornworawat",IsletProject$Sample) ~ "Augsornworawat",
  grepl("MM",IsletProject$Sample) ~ "Wang",
  grepl("JYH",IsletProject$Sample) ~ "Wang",
  grepl("Weng",IsletProject$Sample) ~ "Weng"
)
stopifnot(!anyNA(IsletProject$study))

saveArchRProject(IsletProject)
