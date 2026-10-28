# ==============================================================================
# 053_isletepigenetics_export_to_anndata.R
# Purpose : Export the project (PIC peak matrix, LSI, Harmony, UMAP and cell metadata)
#           to one AnnData (.h5ad) file per sample via Seurat and SeuratDisk. Samples
#           with an existing .h5ad file are skipped.
# Inputs  : data/public/scATACseq/ArchR/IsletProject (with the PIC matrix from 052)
# Outputs : data/public/scATACseq/Project/AnnData/<sample>.h5ad;
#           ranges.bed (peak coordinates, written to the working directory)
# ==============================================================================

data_dir      <- normalizePath("data/public", mustWork = FALSE)     # public/processed data

## Load libraries
library(ArchR)
library(SeuratDisk)
library(Seurat)

## Load the project
proj <- loadArchRProject(file.path(data_dir, "scATACseq/ArchR/IsletProject/"))

anndata_dir <- file.path(data_dir, "scATACseq/Project/AnnData")
dir.create(anndata_dir, recursive = TRUE, showWarnings = FALSE)

## Loop across samples and export
counter <- 1
for (samp in unique(proj$Sample)) {
  if (!file.exists(file.path(anndata_dir, paste0(gsub("-", "", samp), ".h5ad")))) {
    # Subset to only cells from selected donor
    subs <- subsetCells(proj, cellNames = rownames(proj@cellColData)[ proj$Sample == samp])

    # Get and setup the PIC peak matrix; cell names lose '-' and '_', and '#' becomes '_'
    PM <- ArchR::getMatrixFromProject(subs, "PIC")
    counts <- assay(PM)
    colnames(counts) <- gsub("-", "", colnames(counts))
    colnames(counts) <- gsub("_", "", colnames(counts))
    colnames(counts) <- gsub("#", "_", colnames(counts))
    rownames(counts) <- paste(PM@rowRanges@seqnames,":",PM@rowRanges@ranges, sep="")

    # ranges.bed and the temporary tmp.h5Seurat and tmp.h5ad are written to the working directory
    # Write the bed coordinates to a csv file for the first sample
    if (counter == 1) {
      ranges <- as.data.frame(PM@rowRanges)
      colnames(ranges)[1] <- "chr"
      write.table(ranges[,c(1:3)], file = "ranges.bed", quote = FALSE, sep = ",", col.names = TRUE, row.names = FALSE)
    }

    # Get and setup metadata
    metadata <- as.data.frame(subs@cellColData)
    rownames(metadata) <- gsub("-", "", rownames(metadata))
    rownames(metadata) <- gsub("_", "", rownames(metadata))
    rownames(metadata) <- gsub("#", "_", rownames(metadata))
    metadata <- metadata[ match(colnames(counts), rownames(metadata)),]

    # Get and setup LSI
    LSI <- subs@reducedDims$LSI@listData$matSVD
    rownames(LSI) <- gsub("-", "", rownames(LSI))
    rownames(LSI) <- gsub("_", "", rownames(LSI))
    rownames(LSI) <- gsub("#", "_", rownames(LSI))
    LSI <- LSI[ rownames(LSI) %in% colnames(counts),]
    LSI <- LSI[ match(colnames(counts), rownames(LSI)),]
    colnames(LSI) <- paste("LSI_", 1:ncol(LSI), sep="")

    # Get and setup Harmony
    Harm <- subs@reducedDims$Harmony@listData$matDR
    rownames(Harm) <- gsub("-", "", rownames(Harm))
    rownames(Harm) <- gsub("_", "", rownames(Harm))
    rownames(Harm) <- gsub("#", "_", rownames(Harm))
    Harm <- Harm[ rownames(Harm) %in% colnames(counts),]
    Harm <- Harm[ match(colnames(counts), rownames(Harm)),]
    colnames(Harm) <- paste("Harmony_", 1:ncol(Harm), sep="")

    # Get and setup UMAP
    UMAP <- subs@embeddings$HarmonyUMAP@listData$df
    rownames(UMAP) <- gsub("-", "", rownames(UMAP))
    rownames(UMAP) <- gsub("_", "", rownames(UMAP))
    rownames(UMAP) <- gsub("#", "_", rownames(UMAP))
    UMAP <- UMAP[ rownames(UMAP) %in% colnames(counts),]
    UMAP <- UMAP[ match(colnames(counts), rownames(UMAP)),]
    colnames(UMAP) <- paste("UMAP_", 1:ncol(UMAP), sep="")

    # Create a Seurat object
    seu <- Seurat::CreateSeuratObject(counts = counts, project = gsub("-", "", samp), assay = "RNA", min.cells = 0, min.features = 0, meta.data = metadata)
    seu[["LSI"]] <- CreateDimReducObject(LSI, assay = "RNA", key = "LSI_")
    seu[["Harmony"]] <- CreateDimReducObject(Harm, assay = "RNA", key = "Harmony_")
    seu[["UMAP"]] <- CreateDimReducObject(as.matrix(UMAP), assay = "RNA", key = "UMAP_")

    # Create a h5Seurat file
    SaveH5Seurat(seu, "tmp.h5Seurat", overwrite = TRUE)

    # Convert to AnnData
    Convert("tmp.h5Seurat", dest = "h5ad")

    # Clean up
    file.remove("tmp.h5Seurat")
    file.copy(from = "tmp.h5ad", to = file.path(anndata_dir, paste0(gsub("-", "", samp), ".h5ad")), overwrite = TRUE)
    file.remove("tmp.h5ad")
    print(counter)
  }
  counter <- counter + 1
}

