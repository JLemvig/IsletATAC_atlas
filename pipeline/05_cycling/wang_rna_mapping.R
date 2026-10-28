# ==============================================================================
# wang_rna_mapping.R
# Purpose : Build the Wang et al. (GSE200044) multiome snRNA-seq Seurat object with RNA
#           profiles mapped to the ATAC barcodes of the ArchR project (ATAC and RNA
#           libraries matched per donor, RNA barcodes translated to ATAC barcodes with the
#           10x 737K-arc-v1 whitelists), and a subset of the cycling-trajectory cells with
#           their Slingshot pseudotimes.
# Inputs  : data/public/scATACseq/ArchR/IsletProject (cellColData with study, Sample,
#           cycling_labels; HarmonyUMAP);
#           data/public/scATACseq/ArchR/cycling_plus_5k_abdp (slingPseudotime_*, from
#           run_slingshot.R);
#           data/public/scRNAseq/Wang_Sander/Preprocessed/<donor>/Solo.out/GeneFull_Ex50pAS/raw/
#           (STARsolo output of pipeline/01_alignment/align_10x_snrna.sh);
#           data/whitelist/737K-arc-v1.txt and 737K-arc-v1-ATAC.txt (10x multiome gene
#           expression and ATAC barcodes, as in pipeline/01_alignment)
# Outputs : data/public/scATACseq/JM/Wang/RNA_mapped_all_cells_seurat.Rds,
#           data/public/scATACseq/JM/Wang/RNA_mapped_cycling_seurat.Rds
# ==============================================================================

data_dir      <- normalizePath("data/public", mustWork = FALSE)     # public/processed data
whitelist_dir <- normalizePath("data/whitelist", mustWork = FALSE)  # 10x barcode whitelists

suppressPackageStartupMessages({
  library(ArchR)
  library(Matrix)
  library(Seurat)
  library(Biostrings)
})

out_dir <- file.path(data_dir, "scATACseq/JM/Wang")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

proj <- ArchR::loadArchRProject(file.path(data_dir, "scATACseq/ArchR/IsletProject/"))

# Donor label transfer: ATAC library <-> RNA library per sample
RNA_library <- c("MM755", "MM725","MM749","MM740","MM729","MM739","MM738","MM754","MM732","MM752","MM750","MM756","MM730","MM731","MM737","MM751","MM753","MM726","MM727","MM508")
RNA_sample <- c("A0011", "C0019", "C0026","C0024","C0027","C0021","A0028","A0019","A0024","A0029","C0013","C0014","A0030","A0031","A0033","A0021","A0027","C0022","C0023","C0025")
RNA_id <- c(  "SAMP-1613_RNA",  "SAMP-1608_RNA",  "SAMP-1695_RNA",  "SAMP-1605_RNA",  "SAMP-1696_RNA",  "SAMP-1150_RNA",  "SAMP-1138_RNA",  "SAMP-1134_RNA",  "SAMP-1145_RNA",  "SAMP-1139_RNA",  "SAMP-1141_RNA",  "SAMP-1142_RNA",  "SAMP-1500_RNA",  "SAMP-1499_RNA",  "SAMP-1501_RNA",  "SAMP-1609_RNA",  "SAMP-1612_RNA",  "SAMP-1498_RNA",  "SAMP-1503_RNA",  "SAMP-1606_RNA")
ATAC_library <- c("MM747",  "MM717",  "MM741",  "MM736",  "MM721",  "MM735",  "MM734",  "MM746",  "MM724",  "MM744",  "MM742",  "MM748",  "MM722",  "MM723",  "MM733",  "MM743",  "MM745",  "MM718",  "MM719", "MM503")
ATAC_sample <- c(  "A0011",  "C0019",  "C0026",  "C0024",  "C0027",  "C0021",  "A0028",  "A0019",  "A0024",  "A0029",  "C0013",  "C0014",  "A0030",  "A0031",  "A0033",  "A0021",  "A0027",  "C0022",  "C0023", "C0025")
ATAC <- data.frame(Sample = ATAC_sample, ATAC_Library = ATAC_library)
RNA <- data.frame(Sample = RNA_sample, RNA_Library = RNA_library, ID = RNA_id)
transfer <- merge(ATAC, RNA, by = "Sample")
transfer$ATAC_Library <- gsub("MM", "MM_", transfer$ATAC_Library)

# Metadata of the Wang ATAC cells
md_atac <- as.data.frame(proj@cellColData)
md_atac <- md_atac[ md_atac$study == "Wang",]
md_atac$barcode <- rownames(md_atac)
md_atac$ATAC <- rownames(md_atac)
md_atac <- merge(md_atac, transfer, by.x = "Sample", by.y="ATAC_Library")

# Translate ATAC barcodes to RNA barcodes (ATAC whitelist is the reverse complement)
RNA_barcodes <- read.delim(file.path(whitelist_dir, "737K-arc-v1.txt"), header = FALSE)
ATAC_barcodes <- read.delim(file.path(whitelist_dir, "737K-arc-v1-ATAC.txt"), header = FALSE)
BC <- DNAStringSet(x = ATAC_barcodes[,1])
ATAC_barcodes[,"ATAC_barcode"] <- as.data.frame(Biostrings::reverseComplement(BC))
md_atac$barcode <- unlist(lapply(strsplit(md_atac$barcode, "#"), FUN = function(x) { x[2]}))
ATAC_barcodes[,"RNA_barcode"] <- RNA_barcodes[,1]
ATAC_barcodes <- ATAC_barcodes[,2:3]
md_atac <- merge(md_atac, ATAC_barcodes, by.x="barcode", by.y="ATAC_barcode")
md_atac$map <- paste(md_atac$Sample, md_atac$RNA_barcode, sep="#")

# Read STARsolo raw counts per donor, keeping barcodes mapped to ATAC cells
donors <- unique(md_atac$ID)
count_list <- list()
for (donor in donors) {
  solo_dir <- file.path(data_dir, "scRNAseq/Wang_Sander/Preprocessed", donor, "Solo.out/GeneFull_Ex50pAS/raw")
  counts <- readMM(file.path(solo_dir, "matrix.mtx"))
  bcs <- read.delim(file.path(solo_dir, "barcodes.tsv"), header = FALSE)
  feats <- read.delim(file.path(solo_dir, "features.tsv"), header = FALSE)
  counts <- counts[duplicated(feats$V2) == FALSE,]
  feats <- feats[duplicated(feats$V2) == FALSE,]
  rownames(counts) <- feats[,2]
  colnames(counts) <- bcs[,1]
  counts <- counts[,colnames(counts) %in% md_atac[ md_atac$ID %in% donor,"RNA_barcode"]]
  colnames(counts) <- paste(unique(md_atac[ md_atac$ID %in% donor,"Sample"]), colnames(counts), sep="#")
  count_list[[length(count_list) + 1]] <- counts
}

# Seurat object with the ATAC metadata
counts <- do.call("cbind", count_list)
rownames(md_atac) <- md_atac$map
md_atac <- md_atac[,2:80]
counts <-  counts[,colnames(counts) %in% rownames(md_atac)]
md_atac <- md_atac[rownames(md_atac) %in% colnames(counts),]
md_atac <- md_atac[match(colnames(counts), rownames(md_atac)),]
seu <- CreateSeuratObject(CreateAssayObject(counts), meta.data = md_atac)

# ATAC HarmonyUMAP as the embedding
UMAP <- ArchR::getEmbedding(proj, embedding = "HarmonyUMAP")
UMAP <- UMAP[ rownames(UMAP) %in% seu$ATAC,]
UMAP <- UMAP[ match(seu$ATAC, rownames(UMAP)),]
colnames(UMAP) <- c("UMAP_1", "UMAP_2")
rownames(UMAP) <- colnames(seu)
seu[["umap"]] <- CreateDimReducObject(as.matrix(UMAP), assay = "RNA")

seu <- NormalizeData(seu)
saveRDS(seu, file.path(out_dir, "RNA_mapped_all_cells_seurat.Rds"))

# Cycling-trajectory cells (without diffuse) with the Slingshot pseudotimes
cycling_proj <- loadArchRProject(file.path(data_dir, "scATACseq/ArchR/cycling_plus_5k_abdp/"))
cycling <- seu[,seu$ATAC %in% rownames(cycling_proj@cellColData)]
cycling <- subset(cycling, cycling_labels != "diffuse")
md <- as.data.frame(cycling_proj@cellColData)
md <- md[ rownames(md) %in% cycling$ATAC,]
md <- md[ match(cycling$ATAC, rownames(md)),]
cycling$slingPseudotime_beta <- md$slingPseudotime_beta
cycling$slingPseudotime_alpha <- md$slingPseudotime_alpha
cycling$slingPseudotime_delta <- md$slingPseudotime_delta
cycling$slingPseudotime_pp <- md$slingPseudotime_pp

saveRDS(cycling, file.path(out_dir, "RNA_mapped_cycling_seurat.Rds"))
