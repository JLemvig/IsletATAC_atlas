# ==============================================================================
# benaglio_labels.R
# Purpose : Compare our cell type labels with the published Benaglio labels by
#           barcode overlap: match each ArchR sample to its best-overlapping
#           published donor__treatment group and write the barcode-level table.
# Inputs  : data/public/scATACseq/ArchR/IsletProject/, data/public/scATACseq/ArchR/Benaglio_labels.csv.gz
# Outputs : data/public/scATACseq/ArchR/label_comparison/benaglio_barcode_label_comparison.csv
# ==============================================================================

library(ArchR)

data_dir      <- normalizePath("data/public", mustWork = FALSE)     # public/processed data

IsletProject <- loadArchRProject(file.path(data_dir, "scATACseq/ArchR/IsletProject/"))
benaglio_labels <- read.csv(file.path(data_dir, "scATACseq/ArchR/Benaglio_labels.csv.gz"), row.names = 1)

## --- Helpers ---
revcomp <- function(x) {
  x <- toupper(x)
  rev_str <- vapply(strsplit(x, "", fixed = TRUE), function(s) paste(rev(s), collapse=""), character(1))
  chartr("ACGT", "TGCA", rev_str)
}

# "CB:Z:AAAC...-1" -> "AAAC..."
extract_cb <- function(x) sub("^(?:[^:]*:){2}([^:-]+)-.*$", "\\1", x)

# "Benaglio-01#GCTA..." -> "GCTA..."
extract_after_hash <- function(x) sub("^.*#", "", x)

## --- 1) Clean Benaglio label barcodes ---
benaglio_labels$bc <- extract_cb(benaglio_labels$barcode)
benaglio_labels$bc_rc <- revcomp(benaglio_labels$bc)

## --- 2) Build ArchR Benaglio cell table (barcodes + metadata) ---
cc <- as.data.frame(ArchR::getCellColData(IsletProject))

cells_ben <- IsletProject$cellNames[IsletProject$study %in% "Benaglio"]
cells_ben_bc <- extract_after_hash(cells_ben)

cc_ben <- cc[match(cells_ben, rownames(cc)), , drop = FALSE]

archr_df <- data.frame(
  cellName    = cells_ben,
  barcode     = cells_ben_bc,
  islet_label = IsletProject$labels[match(cells_ben, IsletProject$cellNames)],
  stringsAsFactors = FALSE
)

## --- 3) ArchR sample (Benaglio-01 to Benaglio-05) ---
archr_df$archr_ben_group <- as.character(cc_ben[["Sample"]])

## --- 4) Choose forward or reverse-complement barcodes (whichever overlaps more) ---
overlap_forward <- sum(benaglio_labels$bc    %in% archr_df$barcode)
overlap_rc      <- sum(benaglio_labels$bc_rc %in% archr_df$barcode)

use_rc <- overlap_rc > overlap_forward

benaglio_labels$bc_match <- if (use_rc) benaglio_labels$bc_rc else benaglio_labels$bc

message("Barcode matching direction: ", if (use_rc) "reverse-complement" else "forward",
        " (forward=", overlap_forward, ", rc=", overlap_rc, ")")

## --- 5) Benaglio label groups = donor x treatment ---
benaglio_labels$ben_group <- paste(benaglio_labels$donor, benaglio_labels$treatment, sep = "__")

ben_df <- data.frame(
  barcode = benaglio_labels$bc_match,
  donor   = benaglio_labels$donor,
  treatment = benaglio_labels$treatment,
  ben_group = benaglio_labels$ben_group,
  assigned_celltype = benaglio_labels$assigned_celltype,
  stringsAsFactors = FALSE
)

## --- 6) Overlaps per (ArchR sample) x (Benaglio donor__treatment group) ---
ov <- merge(
  archr_df[, c("barcode", "archr_ben_group", "islet_label")],
  ben_df[, c("barcode", "ben_group", "assigned_celltype")],
  by = "barcode",
  all = FALSE
)

pair_counts <- as.data.frame(table(ov$archr_ben_group, ov$ben_group), stringsAsFactors = FALSE)
colnames(pair_counts) <- c("archr_group_1to5", "ben_group", "n_overlap")
pair_counts$n_overlap <- as.integer(pair_counts$n_overlap)
pair_counts <- pair_counts[pair_counts$n_overlap > 0, , drop = FALSE]

archr_totals <- as.data.frame(table(archr_df$archr_ben_group), stringsAsFactors = FALSE)
colnames(archr_totals) <- c("archr_group_1to5", "n_archr")

ben_totals <- as.data.frame(table(ben_df$ben_group), stringsAsFactors = FALSE)
colnames(ben_totals) <- c("ben_group", "n_ben")

pair_counts <- merge(pair_counts, archr_totals, by = "archr_group_1to5", all.x = TRUE)
pair_counts <- merge(pair_counts, ben_totals, by = "ben_group", all.x = TRUE)

pair_counts$frac_of_archr <- pair_counts$n_overlap / pair_counts$n_archr
pair_counts$frac_of_ben   <- pair_counts$n_overlap / pair_counts$n_ben

pair_counts <- pair_counts[order(pair_counts$archr_group_1to5, -pair_counts$n_overlap), ]

## --- 7) Best match per ArchR sample ---
best_match <- do.call(rbind, lapply(split(pair_counts, pair_counts$archr_group_1to5), function(df) {
  df[which.max(df$n_overlap), , drop = FALSE]
}))
best_match <- best_match[order(as.integer(best_match$archr_group_1to5)), ]

## --- 8) Flag non-unique mappings ---
dup_ben_groups <- best_match$ben_group[duplicated(best_match$ben_group)]
dup_ben_groups <- unique(dup_ben_groups)

if (length(dup_ben_groups) > 0) {
  message("WARNING: Some Benaglio donor__treatment groups are the best match for multiple ArchR groups: ",
          paste(dup_ben_groups, collapse = ", "))
}

## --- 9) Restrict overlaps to the matched pairs ---
pairs_to_use <- best_match[, c("archr_group_1to5", "ben_group")]
pairs_to_use$key <- paste(pairs_to_use$archr_group_1to5, pairs_to_use$ben_group, sep="||")

ov$key <- paste(ov$archr_ben_group, ov$ben_group, sep="||")

ov$assigned_celltype_trim <- tolower(sub("_.*$", "", ov$assigned_celltype))

ov_matched <- ov[ov$key %in% pairs_to_use$key, , drop=FALSE]

## --- 10) Export the matched barcode-level comparison ---
comparison_df <- ov_matched[, c("barcode","archr_ben_group","islet_label","ben_group","assigned_celltype_trim")]

dir.create(file.path(data_dir, "scATACseq/ArchR/label_comparison"), recursive = TRUE, showWarnings = FALSE)
write.csv(comparison_df, file.path(data_dir, "scATACseq/ArchR/label_comparison/benaglio_barcode_label_comparison.csv"), row.names = FALSE)
