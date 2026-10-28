# ==============================================================================
# chiou_labels.R
# Purpose : Compare our cell type labels with the published Chiou labels by
#           barcode overlap: match each ArchR sample to its best-overlapping
#           published donor and write the barcode-level table.
# Inputs  : data/public/scATACseq/ArchR/IsletProject/, data/public/scATACseq/ArchR/Chiou_labels.txt
# Outputs : data/public/scATACseq/ArchR/label_comparison/chiou_barcode_label_comparison.csv
# ==============================================================================

library(ArchR)

data_dir      <- normalizePath("data/public", mustWork = FALSE)     # public/processed data

IsletProject <- loadArchRProject(file.path(data_dir, "scATACseq/ArchR/IsletProject/"))
chiou_labels <- read.delim(file.path(data_dir, "scATACseq/ArchR/Chiou_labels.txt"))
chiou_labels <- chiou_labels[chiou_labels$cluster_name %in% c("acinar","activated_stellate","ductal","endothelial","gamma","GCGhigh_alpha","GCGlow_alpha","INShigh_beta","INSlow_beta","quiescent_stellate","SSThigh_delta","SSTlow_delta"),]

## --- Helpers ---
extract_after_hash <- function(x) sub("^.*#", "", x)            # ArchR cellNames -> barcode
extract_between_dash_us <- function(x) sub("^.*-([^_]+)_.*$", "\\1", x)  # between "-" and "_"
extract_after_us <- function(x) sub("^.*_", "", x)              # after "_"

## --- 1) Restrict ArchR to Chiou study only ---
cells_chiou <- IsletProject$cellNames[IsletProject$study %in% "Chiou"]

archr_df <- data.frame(
  cellName    = cells_chiou,
  barcode     = extract_after_hash(cells_chiou),
  archr_group = IsletProject$Sample[match(cells_chiou, IsletProject$cellNames)],
  islet_label = IsletProject$labels[match(cells_chiou, IsletProject$cellNames)],
  stringsAsFactors = FALSE
)

archr_df <- archr_df[!is.na(archr_df$archr_group) & archr_df$archr_group != "" &
                       !is.na(archr_df$barcode) & archr_df$barcode != "", , drop=FALSE]

## --- 2) Clean Chiou labels: donor + barcode from cellName ---
chiou_df <- data.frame(
  cellName = chiou_labels$cellName,
  chiou_group = extract_between_dash_us(chiou_labels$cellName),  # donor between "-" and "_"
  barcode = extract_after_us(chiou_labels$cellName),             # barcode after "_"
  cluster_raw = chiou_labels$cluster_name,
  stringsAsFactors = FALSE
)

# Keep text after "_" and lowercase (e.g. "GCGhigh_alpha" -> "alpha")
chiou_df$celltype_trim <- ifelse(
  is.na(chiou_df$cluster_raw),
  NA_character_,
  tolower(sub("^.*_", "", chiou_df$cluster_raw))
)

chiou_df <- chiou_df[!is.na(chiou_df$chiou_group) & chiou_df$chiou_group != "" &
                       !is.na(chiou_df$barcode) & chiou_df$barcode != "", , drop=FALSE]

## --- 3) Overlaps per (ArchR sample) x (Chiou donor) ---
ov_all <- merge(
  archr_df[, c("barcode","archr_group","islet_label")],
  chiou_df[, c("barcode","chiou_group","celltype_trim")],
  by = "barcode",
  all = FALSE
)

pair_counts <- as.data.frame(table(ov_all$archr_group, ov_all$chiou_group), stringsAsFactors = FALSE)
colnames(pair_counts) <- c("archr_group", "chiou_group", "n_overlap")
pair_counts$n_overlap <- as.integer(pair_counts$n_overlap)
pair_counts <- pair_counts[pair_counts$n_overlap > 0, , drop = FALSE]

archr_totals <- as.data.frame(table(archr_df$archr_group), stringsAsFactors = FALSE)
colnames(archr_totals) <- c("archr_group", "n_archr")

chiou_totals <- as.data.frame(table(chiou_df$chiou_group), stringsAsFactors = FALSE)
colnames(chiou_totals) <- c("chiou_group", "n_chiou")

pair_counts <- merge(pair_counts, archr_totals, by = "archr_group", all.x = TRUE)
pair_counts <- merge(pair_counts, chiou_totals, by = "chiou_group", all.x = TRUE)

pair_counts$frac_of_archr <- pair_counts$n_overlap / pair_counts$n_archr
pair_counts$frac_of_chiou <- pair_counts$n_overlap / pair_counts$n_chiou

pair_counts <- pair_counts[order(pair_counts$archr_group, -pair_counts$n_overlap), ]

## --- 4) Best match per ArchR sample ---
best_match <- do.call(rbind, lapply(split(pair_counts, pair_counts$archr_group), function(df) {
  df[which.max(df$n_overlap), , drop = FALSE]
}))
best_match <- best_match[order(best_match$archr_group), ]

## --- 5) Restrict overlaps to matched pairs ---
pairs_to_use <- best_match[, c("archr_group","chiou_group")]
pairs_to_use$key <- paste(pairs_to_use$archr_group, pairs_to_use$chiou_group, sep="||")

ov <- ov_all
ov$key <- paste(ov$archr_group, ov$chiou_group, sep="||")
ov_matched <- ov[ov$key %in% pairs_to_use$key, , drop=FALSE]

## --- 6) Save barcode-level comparison ---
comparison_df <- ov_matched[, c("barcode","archr_group","islet_label","chiou_group","celltype_trim")]

dir.create(file.path(data_dir, "scATACseq/ArchR/label_comparison"), recursive = TRUE, showWarnings = FALSE)
write.csv(
  comparison_df,
  file.path(data_dir, "scATACseq/ArchR/label_comparison/chiou_barcode_label_comparison.csv"),
  row.names = FALSE
)
