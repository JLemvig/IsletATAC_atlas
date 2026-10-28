# ==============================================================================
# wang_labels.R
# Purpose : Compare our cell type labels with the published Wang labels by
#           barcode overlap: match each ArchR sample to its best-overlapping
#           published donor and write the barcode-level table.
# Inputs  : data/public/scATACseq/ArchR/IsletProject/, data/public/scATACseq/ArchR/Wang_labels_leiden.csv.gz
# Outputs : data/public/scATACseq/ArchR/label_comparison/wang_barcode_label_comparison.csv
# ==============================================================================

library(ArchR)

data_dir      <- normalizePath("data/public", mustWork = FALSE)     # public/processed data

IsletProject <- loadArchRProject(file.path(data_dir, "scATACseq/ArchR/IsletProject/"))
wang_labels <- read.csv(file.path(data_dir, "scATACseq/ArchR/Wang_labels_leiden.csv.gz"), row.names = 1)
leiden_to_celltype <- c("0"="alpha","1"="beta","2"="beta","3"="beta","4"="alpha","5"="beta","6"="delta","7"="acinar","8"="alpha","9"="ductal","10"="gamma","11"="stellate","12"="immune","13"="endothelial")
wang_labels$celltype <- unname(leiden_to_celltype[as.character(wang_labels$leiden)])

## --- Helpers ---
extract_after_hash <- function(x) sub("^.*#", "", x)   # ArchR cellNames -> barcode
extract_after_us   <- function(x) sub("^.*_", "", x)   # Wang rownames -> barcode
rm_us              <- function(x) gsub("_", "", x)     # remove "_" for donor/sample matching

## --- 1) Restrict ArchR to Wang study only ---
cells_wang <- IsletProject$cellNames[IsletProject$study %in% "Wang"]

archr_df <- data.frame(
  cellName    = cells_wang,
  barcode     = extract_after_hash(cells_wang),
  archr_group = rm_us(IsletProject$Sample[match(cells_wang, IsletProject$cellNames)]),
  islet_label = IsletProject$labels[match(cells_wang, IsletProject$cellNames)],
  stringsAsFactors = FALSE
)

archr_df <- archr_df[!is.na(archr_df$archr_group) & archr_df$archr_group != "" &
                       !is.na(archr_df$barcode) & archr_df$barcode != "", , drop=FALSE]

## --- 2) Clean Wang barcodes + celltypes (NA-robust) ---
wang_df <- data.frame(
  barcode = extract_after_us(rownames(wang_labels)),
  wang_group = rm_us(wang_labels$donor),     # group = donor (cleaned)
  celltype_raw = wang_labels$celltype,
  stringsAsFactors = FALSE
)

wang_df$celltype_trim <- ifelse(
  is.na(wang_df$celltype_raw),
  NA_character_,
  tolower(sub("_.*$", "", wang_df$celltype_raw))
)

wang_df <- wang_df[!is.na(wang_df$wang_group) & wang_df$wang_group != "" &
                     !is.na(wang_df$barcode) & wang_df$barcode != "", , drop=FALSE]

## --- 3) Overlaps per (ArchR sample group) x (Wang donor group) ---
ov_all <- merge(
  archr_df[, c("barcode","archr_group","islet_label")],
  wang_df[,  c("barcode","wang_group","celltype_trim")],
  by = "barcode",
  all = FALSE
)

pair_counts <- as.data.frame(table(ov_all$archr_group, ov_all$wang_group), stringsAsFactors = FALSE)
colnames(pair_counts) <- c("archr_group", "wang_group", "n_overlap")
pair_counts$n_overlap <- as.integer(pair_counts$n_overlap)
pair_counts <- pair_counts[pair_counts$n_overlap > 0, , drop = FALSE]

archr_totals <- as.data.frame(table(archr_df$archr_group), stringsAsFactors = FALSE)
colnames(archr_totals) <- c("archr_group", "n_archr")

wang_totals <- as.data.frame(table(wang_df$wang_group), stringsAsFactors = FALSE)
colnames(wang_totals) <- c("wang_group", "n_wang")

pair_counts <- merge(pair_counts, archr_totals, by = "archr_group", all.x = TRUE)
pair_counts <- merge(pair_counts, wang_totals, by = "wang_group", all.x = TRUE)

pair_counts$frac_of_archr <- pair_counts$n_overlap / pair_counts$n_archr
pair_counts$frac_of_wang  <- pair_counts$n_overlap / pair_counts$n_wang

pair_counts <- pair_counts[order(pair_counts$archr_group, -pair_counts$n_overlap), ]

## --- 4) Best match per ArchR sample group ---
best_match <- do.call(rbind, lapply(split(pair_counts, pair_counts$archr_group), function(df) {
  df[which.max(df$n_overlap), , drop = FALSE]
}))
best_match <- best_match[order(best_match$archr_group), ]

## --- 5) Restrict overlaps to matched pairs ---
pairs_to_use <- best_match[, c("archr_group","wang_group")]
pairs_to_use$key <- paste(pairs_to_use$archr_group, pairs_to_use$wang_group, sep="||")

ov <- ov_all
ov$key <- paste(ov$archr_group, ov$wang_group, sep="||")
ov_matched <- ov[ov$key %in% pairs_to_use$key, , drop=FALSE]

## --- 6) Save barcode-level comparison ---
comparison_df <- ov_matched[, c("barcode","archr_group","islet_label","wang_group","celltype_trim")]

dir.create(file.path(data_dir, "scATACseq/ArchR/label_comparison"), recursive = TRUE, showWarnings = FALSE)
write.csv(
  comparison_df,
  file.path(data_dir, "scATACseq/ArchR/label_comparison/wang_barcode_label_comparison.csv"),
  row.names = FALSE
)
