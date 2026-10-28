# ==============================================================================
# hpap_labels.R
# Purpose : Compare our cell type labels with the published HPAP labels by
#           barcode overlap: match each ArchR sample to its best-overlapping
#           published donor library and write the barcode-level table.
# Inputs  : data/public/scATACseq/ArchR/IsletProject/, data/public/scATACseq/ArchR/hpap_atac.h5ad
# Outputs : data/public/scATACseq/ArchR/label_comparison/hpap_barcode_label_comparison.csv
# ==============================================================================

library(ArchR)

suppressPackageStartupMessages({
  library(zellkonverter)
  library(SingleCellExperiment)
})

data_dir      <- normalizePath("data/public", mustWork = FALSE)     # public/processed data

IsletProject <- loadArchRProject(file.path(data_dir, "scATACseq/ArchR/IsletProject/"))
path <- file.path(data_dir, "scATACseq/ArchR/hpap_atac.h5ad")

sce <- readH5AD(path)

## --- Helpers ---
revcomp <- function(x) {
  x <- toupper(x)
  rev_str <- vapply(strsplit(x, "", fixed = TRUE), function(s) paste(rev(s), collapse=""), character(1))
  chartr("ACGT", "TGCA", rev_str)
}
extract_after_hash <- function(x) sub("^.*#", "", x)    # "HPAP-035#BARCODE" -> "BARCODE"
extract_after_us   <- function(x) sub("^.*_", "", x)    # "HPAP-035_BARCODE-1" -> "BARCODE-1"
strip_10x_suffix   <- function(x) sub("-.*$", "", x)    # "BARCODE-1" -> "BARCODE"

map_hpap_celltype <- function(x) {
  m <- c(
    "Alpha"="alpha",
    "Beta"="beta",
    "Delta"="delta",
    "Gamma"="pp",
    "Acinar"="acinar",
    "Ductal"="ductal",
    "MUC5B_Ductal"="ductal",
    "Endothelial"="endothelial",
    "Schwann"="schwann",
    "A_Stellate"="activated_stellate",
    "Q_Stellate"="quiescent_stellate",
    "Immune"="immune"
  )
  out <- unname(m[as.character(x)])
  ifelse(is.na(x), NA_character_, ifelse(is.na(out), tolower(as.character(x)), out))
}

## --- 1) Build HPAP df from sce (h5ad) ---
sce_cells <- colnames(sce)  # e.g. "HPAP-035_AAACG...-1"
hpap_df <- data.frame(
  barcode     = strip_10x_suffix(extract_after_us(sce_cells)),  # "AAAC..."
  hpap_group  = as.character(sce$library),                      # "HPAP-035"
  celltype_raw = as.character(sce$Cell.Type),
  stringsAsFactors = FALSE
)
hpap_df$celltype_trim <- map_hpap_celltype(hpap_df$celltype_raw)

hpap_df <- hpap_df[!is.na(hpap_df$hpap_group) & hpap_df$hpap_group != "" &
                     !is.na(hpap_df$barcode) & hpap_df$barcode != "", , drop=FALSE]

valid_hpap_groups <- sort(unique(hpap_df$hpap_group))  # groups present in the h5ad

## --- 2) ArchR: HPAP study only, and only groups present in the h5ad ---
cells_hpap <- IsletProject$cellNames[IsletProject$study %in% "HPAP"]

archr_groups_all <- IsletProject$Sample[match(cells_hpap, IsletProject$cellNames)]
keep_cells <- archr_groups_all %in% valid_hpap_groups

cells_hpap <- cells_hpap[keep_cells]

archr_df <- data.frame(
  cellName    = cells_hpap,
  barcode     = extract_after_hash(cells_hpap),
  archr_group = IsletProject$Sample[match(cells_hpap, IsletProject$cellNames)],
  islet_label = IsletProject$labels[match(cells_hpap, IsletProject$cellNames)],
  stringsAsFactors = FALSE
)

archr_df <- archr_df[!is.na(archr_df$archr_group) & archr_df$archr_group != "" &
                       !is.na(archr_df$barcode) & archr_df$barcode != "", , drop=FALSE]

## --- 3) Choose forward or reverse-complement barcodes (whichever overlaps more) ---
overlap_forward <- sum(hpap_df$barcode %in% archr_df$barcode)
overlap_rc      <- sum(revcomp(hpap_df$barcode) %in% archr_df$barcode)

use_rc <- overlap_rc > overlap_forward
hpap_df$barcode_match <- if (use_rc) revcomp(hpap_df$barcode) else hpap_df$barcode

message("Barcode matching direction: ", if (use_rc) "reverse-complement" else "forward",
        " (forward=", overlap_forward, ", rc=", overlap_rc, ")")

## --- 4) Overlap df (barcode-level) ---
ov_all <- merge(
  archr_df[, c("barcode","archr_group","islet_label")],
  hpap_df[,  c("barcode_match","hpap_group","celltype_trim")],
  by.x = "barcode", by.y = "barcode_match",
  all = FALSE
)

## --- 5) Pairwise overlap counts + best match per ArchR sample ---
pair_counts <- as.data.frame(table(ov_all$archr_group, ov_all$hpap_group), stringsAsFactors = FALSE)
colnames(pair_counts) <- c("archr_group", "hpap_group", "n_overlap")
pair_counts$n_overlap <- as.integer(pair_counts$n_overlap)
pair_counts <- pair_counts[pair_counts$n_overlap > 0, , drop=FALSE]

archr_totals <- as.data.frame(table(archr_df$archr_group), stringsAsFactors = FALSE)
colnames(archr_totals) <- c("archr_group", "n_archr")

hpap_totals <- as.data.frame(table(hpap_df$hpap_group), stringsAsFactors = FALSE)
colnames(hpap_totals) <- c("hpap_group", "n_hpap")

pair_counts <- merge(pair_counts, archr_totals, by="archr_group", all.x=TRUE)
pair_counts <- merge(pair_counts, hpap_totals, by="hpap_group", all.x=TRUE)

pair_counts$frac_of_archr <- pair_counts$n_overlap / pair_counts$n_archr
pair_counts$frac_of_hpap  <- pair_counts$n_overlap / pair_counts$n_hpap
pair_counts <- pair_counts[order(pair_counts$archr_group, -pair_counts$n_overlap), ]

best_match <- do.call(rbind, lapply(split(pair_counts, pair_counts$archr_group), function(df) {
  df[which.max(df$n_overlap), , drop=FALSE]
}))
best_match <- best_match[order(best_match$archr_group), ]

## --- 6) Restrict to matched pairs ---
pairs_to_use <- best_match[, c("archr_group","hpap_group")]
pairs_to_use$key <- paste(pairs_to_use$archr_group, pairs_to_use$hpap_group, sep="||")

ov <- ov_all
ov$key <- paste(ov$archr_group, ov$hpap_group, sep="||")
ov_matched <- ov[ov$key %in% pairs_to_use$key, , drop=FALSE]

## --- 7) Save barcode-level comparison ---
comparison_df <- ov_matched[, c("barcode","archr_group","islet_label","hpap_group","celltype_trim")]

dir.create(file.path(data_dir, "scATACseq/ArchR/label_comparison"), recursive = TRUE, showWarnings = FALSE)
write.csv(
  comparison_df,
  file.path(data_dir, "scATACseq/ArchR/label_comparison/hpap_barcode_label_comparison.csv"),
  row.names = FALSE
)
