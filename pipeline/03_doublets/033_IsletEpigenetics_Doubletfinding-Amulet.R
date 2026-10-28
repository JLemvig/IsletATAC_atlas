# ==============================================================================
# 033_IsletEpigenetics_Doubletfinding-Amulet.R
# Purpose : Extract the fragments of each sample from its ArrowFile, run AMULET
#           (scDblFinder::amulet) and add the p and q values to IsletProject3.
# Inputs  : data/public/scATACseq/ArchR/IsletProject3 (ArrowFiles)
# Outputs : data/public/scATACseq/ArchR/amuletdoublets.rds;
#           IsletProject3$Amulet.p.value, IsletProject3$Amulet.q.value
# ==============================================================================

data_dir      <- normalizePath("data/public", mustWork = FALSE)     # public/processed data

## Load dependencies
library(scDblFinder)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ArchR)

## Set ArchR threads and reference genome
addArchRThreads(threads = round(parallelly::availableCores()/4))
addArchRGenome(genome = "hg38")

# Load ArchRProject
IsletProject3 <- loadArchRProject(path = file.path(data_dir, "scATACseq/ArchR/IsletProject3"))

# Create Fragment files and run Amulet (one ArrowFile per sample)
amuletdoublets <- NULL
for(i in list.files(file.path(data_dir, "scATACseq/ArchR/IsletProject3/ArrowFiles/"))){
  fragment <- getFragmentsFromArrow(ArrowFile = paste0(file.path(data_dir, "scATACseq/ArchR/IsletProject3/ArrowFiles/"),i))
  # amulet() takes the cell barcode from the name column
  fragment$name <- fragment$RG
  score <- amulet(fragment)
  amuletdoublets <- dplyr::bind_rows(amuletdoublets,score)
  print(i)
}

saveRDS(amuletdoublets, file = file.path(data_dir, "scATACseq/ArchR/amuletdoublets.rds"))

# Order results as the cells of the project and add to the project
nameorder <- data.frame("cellNames" = IsletProject3$cellNames)
ad <- data.frame("cellNames" = row.names(amuletdoublets),
                 "Amulet.p.value" = amuletdoublets$p.value,
                 "Amulet.q.value" = amuletdoublets$q.value)
df <- dplyr::left_join(nameorder,ad, by = "cellNames")
IsletProject3$Amulet.p.value <- df$Amulet.p.value
IsletProject3$Amulet.q.value <- df$Amulet.q.value

# AMULET results for sample Chiou-04 are not used: the sample was mis-marked by AMULET,
# so its p and q values are set to 1 (no doublets)
IsletProject3$Amulet.q.value[IsletProject3$Sample %in% "Chiou-04"] <- 1
IsletProject3$Amulet.p.value[IsletProject3$Sample %in% "Chiou-04"] <- 1

IsletProject3 <- saveArchRProject(IsletProject3)
