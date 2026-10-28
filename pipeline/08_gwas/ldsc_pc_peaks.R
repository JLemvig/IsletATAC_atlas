# ==============================================================================
# ldsc_pc_peaks.R
# Purpose : LD score regression (partitioned heritability) of the peaks associated with
#           each donor PC: BED files of the positive, negative and combined
#           significant peaks per PC, lifted to hg19 and widened to 1 kb windows,
#           LDSC annotations and LD scores, and --h2 with baselineLD v2.3 for each
#           GWAS trait.
# Inputs  : data/public/scATACseq/correlations/<celltype>/deseq_lrt_joint_signif_<celltype>.rds
#           and peakset_<celltype>.rds (07_peak_associations);
#           data/project/gwas_environment (LDSC, 1000G Phase 3 reference files,
#           baselineLD v2.3, munged GWAS summary statistics)
# Outputs : data/public/scATACseq/correlations/<celltype>/bed/<PC>_{positive,negative,combined}.bed;
#           <celltype>/ldsc/annotation/; <celltype>/ldsc/results/<annotation>/<GWAS>.results;
#           <celltype>/ldsc/results/<celltype>_results.rds (all results of the cell type)
# Requires the LDSC conda environment (env/ldsc.yml) and GNU parallel.
# ==============================================================================

data_dir      <- normalizePath("data/public", mustWork = FALSE)     # public/processed data
tools_dir     <- normalizePath("data/project/gwas_environment", mustWork = FALSE)  # LDSC, reference files, sumstats

library(dplyr)
library(stringr)
library(GenomicRanges)
library(rtracklayer)
library(data.table)

alpha <- 0.05  # adjusted p-value cutoff

celltypes <- c("beta", "alpha", "acinar", "ductal",
               "delta", "pp", "cycling", "activated_stellate", "quiescent_stellate",
               "endothelial", "myeloid", "lymphoid", "mast", "epsilon", "schwann")


## 1) BED files of the significant peaks per PC and direction
for (celltype in celltypes) {

  path     <- file.path(file.path(data_dir, "scATACseq/correlations"), celltype)

  joint_file <- file.path(path, sprintf("deseq_lrt_joint_signif_%s.rds", celltype))
  peak_file  <- file.path(path, sprintf("peakset_%s.rds", celltype))

  if (!file.exists(joint_file)) {
    message("No joint results file for ", celltype, ", skipping.")
    next
  }
  if (!file.exists(peak_file)) {
    message("No peakset file for ", celltype, ", skipping.")
    next
  }

  joint_res <- readRDS(joint_file)
  peakset   <- readRDS(peak_file)

  # Ensure a feature column in peakset for matching
  if (!"feature" %in% colnames(mcols(peakset))) {
    mcols(peakset)$feature <- rownames(peakset)
  }

  # make sure joint_res has a feature column with same IDs
  if (!"feature" %in% colnames(joint_res)) {
    stop("joint_res is missing a 'feature' column for celltype: ", celltype)
  }

  # bed output folder
  bed_dir <- file.path(path, "bed")
  dir.create(bed_dir, showWarnings = FALSE, recursive = TRUE)

  pcs <- sort(unique(joint_res$PC))

  for (pc in pcs) {
    for (dirn in c("positive", "negative")) {

      sub <- joint_res %>%
        filter(
          PC == pc,
          direction == dirn,
          !is.na(padj),
          padj < alpha
        )

      if (nrow(sub) == 0L) {
        message("No sig peaks for ", celltype,
                ", PC", pc, " ", dirn, ", skipping.")
        next
      }

      # match peaks by feature id
      peak_idx <- peakset$rownames %in% sub$feature
      gr_sub   <- peakset[peak_idx]

      if (length(gr_sub) == 0L) {
        message("No matching peaks in peakset for ", celltype,
                ", ", pc, " ", dirn, ", skipping.")
        next
      }

      # build BED (0-based start)
      bed_df <- data.frame(
        chr   = as.character(seqnames(gr_sub <- gr_sub)), # ensure character
        start = start(gr_sub) - 1L,   # BED start is 0-based
        end   = end(gr_sub),
        name  = gr_sub$rownames
      )

      # sanity: drop any NA or weird entries
      bed_df <- bed_df %>%
        filter(!is.na(chr), !is.na(start), !is.na(end))

      out_bed <- file.path(
        bed_dir,
        sprintf("%s_%s.bed", pc, dirn)
      )

      write.table(
        bed_df,
        file      = out_bed,
        sep       = "\t",
        quote     = FALSE,
        row.names = FALSE,
        col.names = FALSE
      )

      message("Wrote ", nrow(bed_df),
              " peaks to ", out_bed)
    }
  }
}

## 2) Combined BED per PC: union of the positive and negative peaks
base_root <- file.path(data_dir, "scATACseq/correlations")

for (celltype in celltypes) {

  bed_dir <- file.path(base_root, celltype, "bed")
  if (!dir.exists(bed_dir)) {
    message("[", celltype, "] no bed directory, skipping.")
    next
  }

  bed_files <- list.files(bed_dir, pattern = "\\.bed$", full.names = TRUE)
  if (length(bed_files) == 0) {
    message("[", celltype, "] no bed files, skipping.")
    next
  }

  meta <- data.table(
    file = bed_files,
    base = basename(bed_files)
  )[
    ,
    `:=`(
      pc  = sub("_(positive|negative)\\.bed$", "", base),
      dir = fifelse(grepl("positive", base), "positive",
                    fifelse(grepl("negative", base), "negative", NA_character_))
    )
  ]

  pcs <- unique(meta$pc)

  message("=== ", celltype, " ===")
  for (pc_id in pcs) {
    sub_meta <- meta[pc == pc_id]

    # require both positive and negative
    if (!all(c("positive", "negative") %in% sub_meta$dir)) next

    pos_file <- sub_meta[dir == "positive", file][1]
    neg_file <- sub_meta[dir == "negative", file][1]

    dt_pos <- fread(pos_file, header = FALSE)
    dt_neg <- fread(neg_file, header = FALSE)

    dt_combined <- unique(rbind(dt_pos, dt_neg, use.names = FALSE))

    out_file <- file.path(bed_dir, sprintf("%s_combined.bed", pc_id))
    fwrite(dt_combined, out_file, sep = "\t", col.names = FALSE)

    message("  ", pc_id, " -> combined")
  }
}

## 3) LDSC annotations and LD scores
# BED windows are lifted hg38 -> hg19 because the 1000G Phase3 / baselineLD reference is hg19
path <- system.file(package="liftOver", "extdata", "hg38ToHg19.over.chain")
ch <- import.chain(path)

## parameters for window and minimum peaks per annotation
window_size <- 1000
min_peaks   <- 200L

# Lift a BED file to hg19 and resize to window_size around the centre; NULL if fewer
# than min_peaks remain
process_bed_to_hg19 <- function(bedfile, ch, window_size, min_peaks) {
  bed <- read.delim(bedfile, header = FALSE)
  gr  <- GRanges(seqnames = bed$V1, ranges = IRanges(start = bed$V2, end = bed$V3))

  # liftover
  lo   <- liftOver(gr, ch)
  keep <- elementNROWS(lo) == 1
  if (!any(keep)) return(NULL)
  gr_hg19 <- unlist(lo[keep], use.names = FALSE)

  # center + window
  gr_hg19 <- resize(gr_hg19, width = window_size, fix = "center")

  if (length(gr_hg19) < min_peaks) return(NULL)
  gr_hg19
}

make_annots_and_l2 <- function(
    sel_idx,                            # numeric indices into anno.list, e.g. c(1,3,6)
    anno.list,                          # named GRanges list
    out_dir,                            # where to write 1000G.EUR.QC.<chr>.annot.gz
    bfile_prefix = file.path(tools_dir, "1000G_EUR_Phase3_plink/1000G.EUR.QC."),
    hm3_list    = file.path(tools_dir, "hm3_no_MHC.list.txt"),
    ldsc_py     = file.path(tools_dir, "ldsc/ldsc.py"),
    # conda env 'ldsc' (python2, env/ldsc.yml) in tools_dir/miniconda3, with GNU parallel
    conda_sh    = file.path(tools_dir, "miniconda3/etc/profile.d/conda.sh"),
    env_name    = "ldsc",
    n_jobs      = 8,
    ld_wind_cm  = 1
) {
  sel_idx <- as.integer(sel_idx)
  if (anyNA(sel_idx)) stop("sel_idx contains non-numeric / NA values")
  if (!length(sel_idx)) stop("sel_idx is empty")
  if (any(sel_idx < 1 | sel_idx > length(anno.list))) stop("sel_idx out of range")

  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  sel_names <- names(anno.list)[sel_idx]
  message("Selected annotations (", length(sel_names), "): ",
          paste(sel_names, collapse = ", "))

  # 1) Build per-chr .annot.gz with selected columns
  for (chr in 1:22) {
    # Read BIM and prep base columns CHR,BP,SNP,CM
    bim <- fread(sprintf("%s%d.bim", bfile_prefix, chr),
                 col.names = c("CHR","SNP","CM","BP","A1","A2"))[, .(CHR,BP,SNP,CM)]

    # Genomic positions for overlap
    bim_gr <- GRanges(seqnames = paste0("chr", bim$CHR),
                      ranges   = IRanges(bim$BP, width = 1),
                      rsid     = bim$SNP)

    # Build annotation matrix
    mat <- matrix(0L, nrow = nrow(bim), ncol = length(sel_idx),
                  dimnames = list(NULL, sel_names))
    for (j in seq_along(sel_idx)) {
      hits <- unique(queryHits(findOverlaps(bim_gr, anno.list[[sel_idx[j]]], ignore.strand = TRUE)))
      if (length(hits)) mat[hits, j] <- 1L
    }

    out <- cbind(bim, as.data.table(mat))
    fn  <- file.path(out_dir, sprintf("1000G.EUR.QC.%d.annot", chr))
    fwrite(out, fn, sep = "\t", quote = FALSE)
    system2("gzip", c("-f", fn))
  }

  # 2) Run LDSC --l2 on the files we just wrote (uses Conda env)
  cmd <- sprintf(
    paste(
      "bash -c 'source %s && conda activate %s &&",
      "parallel -j %d python2 %s",
      "--l2",
      "--bfile %s{1}",
      "--ld-wind-cm %d",
      "--annot %s/1000G.EUR.QC.{1}.annot.gz",
      "--out %s/1000G.EUR.QC.{1}",
      "--print-snps %s ::: {1..22}'"
    ),
    conda_sh, env_name, n_jobs, ldsc_py,
    bfile_prefix, ld_wind_cm,
    out_dir, out_dir, hm3_list
  )

  message("Running LDSC --l2 via Conda...")
  status <- system(cmd)
  if (status != 0) stop("LDSC --l2 failed with status ", status)

  invisible(list(
    annot_dir = out_dir,
    selected  = sel_names,
    l2_prefix = file.path(out_dir, "1000G.EUR.QC.")
  ))
}

for (celltype in celltypes) {
  message("=== Processing cell type: ", celltype, " ===")

  base_path <- file.path(file.path(data_dir, "scATACseq/correlations"), celltype)
  bed_dir   <- file.path(base_path, "bed")

  if (!dir.exists(bed_dir)) {
    message("No bed folder for ", celltype, ", skipping.")
    next
  }

  bed_files <- list.files(bed_dir, pattern = "\\.bed$", full.names = TRUE)
  if (!length(bed_files)) {
    message("No .bed files for ", celltype, ", skipping.")
    next
  }

  # Parse PC and direction from file names like "PC4_positive.bed"
  bed_info <- tibble(
    file = bed_files,
    base = basename(bed_files)
  ) %>%
    mutate(
      base = str_remove(base, "\\.bed$"),
      PC   = str_extract(base, "PC[0-9]+|[0-9]+"),
      dirn = case_when(
        str_detect(base, "positive") ~ "positive",
        str_detect(base, "negative") ~ "negative",
        str_detect(base, "combined") ~ "combined",
        TRUE ~ NA_character_
      )
    )

  if (all(is.na(bed_info$dirn))) {
    message("Could not detect pos/neg in filenames for ", celltype, ", skipping.")
    next
  }

  anno.list <- list()

  # One annotation per BED file (PC x positive/negative/combined)
  for (k in seq_len(nrow(bed_info))) {
    this_file <- bed_info$file[k]
    this_pc   <- bed_info$PC[k]
    this_dir  <- bed_info$dirn[k]

    if (is.na(this_pc) || is.na(this_dir)) next

    gr <- process_bed_to_hg19(
      bedfile     = this_file,
      ch          = ch,
      window_size = window_size,
      min_peaks   = min_peaks
    )

    if (is.null(gr)) {
      message("  ", basename(this_file),
              " had too few peaks after liftover/window; skipping.")
      next
    }

    nm <- sprintf("%s_%s_%s", celltype, this_pc, this_dir)
    anno.list[[nm]] <- gr
    message("  Added annotation: ", nm, " (", length(gr), " windows)")
  }

  if (!length(anno.list)) {
    message("No valid annotations (pos/neg) for ", celltype, ", skipping LDSC.")
    next
  }

  # Save anno.list and run LDSC --l2 for all annotations of this cell type
  annotation_folder <- file.path(base_path, "ldsc", "annotation")
  dir.create(annotation_folder, recursive = TRUE, showWarnings = FALSE)

  rds_path <- file.path(
    annotation_folder,
    sprintf("annotation_%d.rds", window_size)
  )
  saveRDS(anno.list, rds_path)
  message("Saved ", length(anno.list), " annotations to: ", rds_path)

  anno.list_loaded <- readRDS(rds_path)

  message("Launching LDSC --l2 for ", celltype, " ...")
  res <- make_annots_and_l2(
    sel_idx   = seq_along(anno.list_loaded),  # every PC pos, PC neg, PC combined
    anno.list = anno.list_loaded,
    out_dir   = annotation_folder
  )

  message("Done with ", celltype, ". LDSC output prefix: ",
          res$l2_prefix, " (", length(res$selected), " annotations)")
}

## 4) Partitioned heritability (--h2) per annotation and GWAS trait
ss <- c(
  "Metabolic_UKB-HbA1c" = file.path(tools_dir, "sumstats/UKB_460K.biochemistry_HbA1c.sumstats.gz"),
  "Metabolic_UKB-T2D" = file.path(tools_dir, "sumstats/UKB_460K.disease_T2D.sumstats.gz"),
  "Metabolic_Suzuki-T2D" = file.path(tools_dir, "sumstats/Suzuki.ALL-T2D.sumstats.gz"),
  "Metabolic_UKB-BMI" = file.path(tools_dir, "sumstats/UKB_460K.body_BMIz.sumstats.gz"),
  "Metabolic_UKB-WHR" = file.path(tools_dir, "sumstats/UKB_460K.body_WHRadjBMIz.sumstats.gz"),
  "Metabolic_UKB-CRP" = file.path(tools_dir, "sumstats/UKB_460K.biochemistry_CreactiveProtein.sumstats.gz"),
  "Metabolic_UKB-Glucose" = file.path(tools_dir, "sumstats/UKB_460K.biochemistry_Glucose.sumstats.gz"),
  "Metabolic_UKB-IGF1" = file.path(tools_dir, "sumstats/UKB_460K.biochemistry_IGF1.sumstats.gz"),
  "Metabolic_UKB-HDL" = file.path(tools_dir, "sumstats/UKB_460K.biochemistry_HDLcholesterol.sumstats.gz"),
  "Metabolic_UKB-LDL" = file.path(tools_dir, "sumstats/UKB_460K.biochemistry_LDLdirect.sumstats.gz"),
  "Metabolic_UKB-TG" = file.path(tools_dir, "sumstats/UKB_460K.biochemistry_Triglycerides.sumstats.gz"),
  "Control_UKB-Asthma" = file.path(tools_dir, "sumstats/UKB_460K.disease_ASTHMA_DIAGNOSED.sumstats.gz"),
  "Control_UKB-Lefthanded" = file.path(tools_dir, "sumstats/UKB_460K.body_LEFT_HANDED.sumstats.gz"),
  "Control_UKB-Balding" = file.path(tools_dir, "sumstats/UKB_460K.body_BALDING4.sumstats.gz"),
  "Control_UKB-Hair" = file.path(tools_dir, "sumstats/UKB_460K.pigment_HAIR.sumstats.gz"),
  "Control_UKB-Skin" = file.path(tools_dir, "sumstats/UKB_460K.pigment_SKIN.sumstats.gz")
)

.ensure_dot <- function(x) {
  x <- as.character(x)
  sub("(?<!\\.)$", ".", x, perl = TRUE)
}

# List which annotations are available in a set of .annot.gz files
list_available_annots <- function(annot_dir,
                                  prefix      = "1000G.EUR.QC.",
                                  chr         = 1,
                                  base_cols   = c("CHR","BP","SNP","CM")) {
  f <- file.path(annot_dir, sprintf("%s%d.annot.gz", .ensure_dot(prefix), chr))
  if (!file.exists(f)) stop("Annotation file not found: ", f)
  dt <- fread(cmd = sprintf("zcat -f %s", shQuote(f)), nrows = 1)
  setdiff(names(dt), base_cols)
}

# Subset an existing annotation set to only specific columns
subset_annots_gz <- function(in_dir,
                             out_dir,
                             keep_annots,
                             prefix     = "1000G.EUR.QC.",
                             chrs       = 1:22,
                             base_cols  = c("CHR","BP","SNP","CM"),
                             ldscore_suffix = ".l2.ldscore.gz") {

  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  # validate availability on one chromosome (assume consistent across chrs)
  avail <- list_available_annots(in_dir, prefix, chr = chrs[1], base_cols = base_cols)
  if (!all(keep_annots %in% avail)) {
    missing <- setdiff(keep_annots, avail)
    stop("These requested annotations are not present: ", paste(missing, collapse = ", "))
  }

  # alias helper (ldscore sometimes appends L2/.L2/_L2)
  ann_aliases <- function(a) unique(c(a, paste0(a, "L2"), paste0(a, ".L2"), paste0(a, "_L2")))
  all_ann_aliases <- unique(unlist(lapply(avail, ann_aliases)))

  map_one_col <- function(a, cols) {
    hit <- intersect(ann_aliases(a), cols)
    if (length(hit)) hit[1] else NA_character_
  }

  # read helper for .M / .M_5_50 (supports optional .gz)
  read_M_vec <- function(path) {
    if (file.exists(path)) {
      scan(path, what = double(), quiet = TRUE)
    } else if (file.exists(paste0(path, ".gz"))) {
      con <- pipe(sprintf("zcat -f %s", shQuote(paste0(path, ".gz"))))
      on.exit(close(con), add = TRUE)
      scan(con, what = double(), quiet = TRUE)
    } else {
      NULL
    }
  }

  for (chr in chrs) {
    # --- subset .annot.gz ---
    fin_annot  <- file.path(in_dir,  sprintf("%s%d.annot.gz", .ensure_dot(prefix), chr))
    fout_annot <- file.path(out_dir, sprintf("%s%d.annot",    .ensure_dot(prefix), chr))
    if (!file.exists(fin_annot)) stop("Missing input: ", fin_annot)

    dt_annot  <- fread(cmd = sprintf("zcat -f %s", shQuote(fin_annot)))
    sub_annot <- dt_annot[, c(base_cols, keep_annots), with = FALSE]
    fwrite(sub_annot, fout_annot, sep = "\t", quote = FALSE)
    system2("gzip", c("-f", fout_annot))

    # --- subset .l2.ldscore.gz ---
    fin_l2  <- file.path(in_dir,  sprintf("%s%d%s", .ensure_dot(prefix), chr, ldscore_suffix))
    fout_l2 <- file.path(out_dir, sprintf("%s%d.l2.ldscore", .ensure_dot(prefix), chr))
    if (!file.exists(fin_l2)) stop("Missing input: ", fin_l2)

    dt_l2 <- fread(cmd = sprintf("zcat -f %s", shQuote(fin_l2)))

    base_l2_cols   <- c("CHR","SNP","BP")
    annot_l2_cols  <- setdiff(names(dt_l2), base_l2_cols)

    if (length(avail) == 1L &&
        length(annot_l2_cols) == 1L &&
        identical(annot_l2_cols, "L2")) {

      # Special case: one annotation, LDSC wrote generic column 'L2'
      ann_cols_full  <- "L2"
      keep_ld_cols   <- "L2"
      non_annot_cols <- setdiff(names(dt_l2), "L2")

    } else {

    # order of ALL annotation columns in ldscore (aligned to 'avail' order)
    ann_cols_full <- vapply(avail, map_one_col, character(1), cols = names(dt_l2))
    if (anyNA(ann_cols_full)) {
      bad <- avail[is.na(ann_cols_full)]
      stop("Could not find ldscore columns for: ", paste(bad, collapse = ", "))
    }

    # columns to keep in ldscore
    keep_ld_cols <- vapply(keep_annots, map_one_col, character(1), cols = names(dt_l2))
    if (anyNA(keep_ld_cols)) {
      bad <- keep_annots[is.na(keep_ld_cols)]
      stop("Requested annotations not present in ldscore columns: ", paste(bad, collapse = ", "))
    }

    non_annot_cols <- setdiff(names(dt_l2), all_ann_aliases)
    }

    sub_l2 <- dt_l2[, c(non_annot_cols, keep_ld_cols), with = FALSE]
    fwrite(sub_l2, fout_l2, sep = "\t", quote = FALSE)
    system2("gzip", c("-f", fout_l2))

    # --- subset .l2.M and .l2.M_5_50 (by annotation indices) ---
    keep_idx <- match(keep_ld_cols, ann_cols_full)  # positions in M vectors

    fin_M      <- file.path(in_dir,  sprintf("%s%d.l2.M",       .ensure_dot(prefix), chr))
    fin_M55    <- file.path(in_dir,  sprintf("%s%d.l2.M_5_50",  .ensure_dot(prefix), chr))
    fout_M     <- file.path(out_dir, sprintf("%s%d.l2.M",       .ensure_dot(prefix), chr))
    fout_M55   <- file.path(out_dir, sprintf("%s%d.l2.M_5_50",  .ensure_dot(prefix), chr))

    mv  <- read_M_vec(fin_M)
    if (!is.null(mv)) {
      if (length(mv) < max(keep_idx)) stop("M file has fewer entries than annotations for chr ", chr)
      cat(paste(mv[keep_idx], collapse = "\t"), "\n", file = fout_M)
    }

    mv55 <- read_M_vec(fin_M55)
    if (!is.null(mv55)) {
      if (length(mv55) < max(keep_idx)) stop("M_5_50 file has fewer entries than annotations for chr ", chr)
      cat(paste(mv55[keep_idx], collapse = "\t"), "\n", file = fout_M55)
    }
  }

  invisible(list(
    annot_dir   = out_dir,
    prefix_used = .ensure_dot(prefix),
    kept        = keep_annots
  ))
}

run_partitioned_h2 <- function(
    sumstats,                                     # various formats (see notes)
    ref_ld_prefixes,                              # character vector of prefixes (WITHOUT chr number; WITH or WITHOUT trailing '.')
    weights_prefix,                               # e.g. "path/to/weights."
    frq_prefix,                                   # e.g. "path/to/1000G.mac5eur."
    out_dir = NULL,                               # directory to write outputs when using multiple sumstats or named inputs
    out_prefix = NULL,                            # kept for backward-compat (single run)
    ldsc_py     = file.path(tools_dir, "ldsc/ldsc.py"),
    # conda env 'ldsc' (python2, env/ldsc.yml) in tools_dir/miniconda3
    conda_sh    = file.path(tools_dir, "miniconda3/etc/profile.d/conda.sh"),
    env_name    = "ldsc",
    print_coefficients = FALSE,
    overlap_annot      = TRUE
) {

  .ensure_dot <- function(x) sub("(?<!\\.)$", ".", as.character(x), perl = TRUE)

  # normalize sumstats input to a data.frame(name, path, out_prefix)
  normalize_sumstats <- function(ss, out_dir, out_prefix) {
    # case: single character path + out_prefix (legacy)
    if (is.character(ss) && length(ss) == 1L && !is.null(out_prefix)) {
      if (!file.exists(ss)) stop("sumstats not found: ", ss)
      return(data.table::data.table(name=basename(out_prefix), path=ss, out_prefix=out_prefix))
    }

    # named character vector
    if (is.character(ss) && length(ss) >= 1L && !is.null(names(ss)) && all(nzchar(names(ss)))) {
      dt <- data.table::data.table(name = names(ss), path = unname(ss))
      if (is.null(out_dir)) stop("Provide out_dir when passing a named vector of sumstats.")
      dt[, out_prefix := file.path(out_dir, name)]
      return(dt)
    }

    # data.frame/data.table with name, path
    if (is.data.frame(ss) && all(c("name","path") %in% names(ss))) {
      dt <- data.table::as.data.table(ss)[, .(name=as.character(name), path=as.character(path))]
      if (is.null(out_dir)) stop("Provide out_dir when passing a data.frame/list of sumstats.")
      dt[, out_prefix := file.path(out_dir, name)]
      return(dt)
    }

    # list of lists: each with name, path
    if (is.list(ss) && length(ss) >= 1L && all(vapply(ss, is.list, TRUE))) {
      dt <- data.table::rbindlist(lapply(ss, function(x) {
        if (!all(c("name","path") %in% names(x))) stop("Each list element must have 'name' and 'path'.")
        data.table::data.table(name = as.character(x$name), path = as.character(x$path))
      }))
      if (is.null(out_dir)) stop("Provide out_dir when passing a list of sumstats.")
      dt[, out_prefix := file.path(out_dir, name)]
      return(dt)
    }

    stop("Unsupported 'sumstats' input. Use a named character vector, a data.frame with columns (name, path), a list of lists, or single path + out_prefix.")
  }

  # sanity
  if (!(length(ref_ld_prefixes) >= 1)) stop("ref_ld_prefixes must have length >= 1")
  ref_ld_chr_arg <- paste(vapply(ref_ld_prefixes, .ensure_dot, ""), collapse = ",")

  # normalize
  ss_dt <- normalize_sumstats(sumstats, out_dir, out_prefix)

  # check files + prepare output dir
  ss_dt[, exists := file.exists(path)]
  if (any(!ss_dt$exists)) {
    bad <- ss_dt[!exists]
    stop("These sumstats paths do not exist:\n",
         paste(sprintf(" - %s (%s)", bad$name, bad$path), collapse = "\n"))
  }
  if (!is.null(out_dir)) dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  unique_dirs <- unique(dirname(ss_dt$out_prefix))
  for (d in unique_dirs) dir.create(d, showWarnings = FALSE, recursive = TRUE)

  # common flags (except --h2 and --out which vary per sumstats)
  common_flags <- c(
    sprintf("--ref-ld-chr %s",  shQuote(ref_ld_chr_arg)),
    sprintf("--w-ld-chr %s",    shQuote(.ensure_dot(weights_prefix))),
    sprintf("--frqfile-chr %s", shQuote(.ensure_dot(frq_prefix)))
  )
  if (isTRUE(overlap_annot))       common_flags <- c(common_flags, "--overlap-annot")
  if (isTRUE(print_coefficients))  common_flags <- c(common_flags, "--print-coefficients")

  run_one <- function(sumstats_path, out_pref) {
    flags <- c(
      sprintf("--h2 %s",  shQuote(sumstats_path)),
      common_flags,
      sprintf("--out %s", shQuote(out_pref))
    )
    cmd <- sprintf(
      "bash -c 'source %s && conda activate %s && python2 %s %s'",
      shQuote(conda_sh),
      shQuote(env_name),
      shQuote(ldsc_py),
      paste(flags, collapse = " ")
    )
    message("Running LDSC --h2 for: ", out_pref)
    status <- system(cmd)
    if (status != 0) stop("LDSC --h2 failed (", basename(out_pref), ") with status ", status)
    list(results = paste0(out_pref, ".results"), log = paste0(out_pref, ".log"))
  }

  # run all
  outs <- lapply(seq_len(nrow(ss_dt)), function(i) run_one(ss_dt$path[i], ss_dt$out_prefix[i]))
  ss_dt[, `:=`(out_results = vapply(outs, `[[`, "", "results"),
               out_log     = vapply(outs, `[[`, "", "log"))]

  invisible(ss_dt[])
}

force <- FALSE

for (celltype in celltypes) {
  message("=== Cell type: ", celltype, " ===")

  base_outdir      <- file.path(file.path(data_dir, "scATACseq/correlations"), celltype)
  full_annot_dir   <- file.path(base_outdir, "ldsc", "annotation")
  if (!dir.exists(full_annot_dir)) {
    message("  No annotation directory, skipping.")
    next
  }

  # All available annotation columns in this custom panel
  avail <- list_available_annots(
    annot_dir = full_annot_dir,
    prefix    = "1000G.EUR.QC."
  )

  # Restrict to annotations for this cell type (names like "beta_PC4_positive")
  cell_prefix <- paste0(celltype, "_")
  cell_annots <- grep(paste0("^", cell_prefix), avail, value = TRUE)

  if (!length(cell_annots)) {
    message("  No annotations for ", celltype, " in ", full_annot_dir, ", skipping.")
    next
  }

  annot_df <- data.table(annot = cell_annots)[
    # 1) split into celltype and the rest ("PC9_negative")
    , c("celltype", "rest") := tstrsplit(annot, "_(?=PC[0-9]+_)", perl = TRUE)
  ][
    # 2) split rest into pc_label and dirn
    , c("pc_label", "dirn") := tstrsplit(rest, "_", fixed = TRUE)
  ][
    # 3) drop helper column
    , rest := NULL
  ]

  # pc_label is like "PC4"; get numeric part for ordering
  annot_df[
    , pc_num := as.integer(str_extract(pc_label, "[0-9]+"))
  ]

  # enforce order: positive, negative, combined
  annot_df[
    , dirn := factor(dirn, levels = c("positive", "negative", "combined"))
  ]

  # Order: PC4 pos, PC4 neg, PC4 combined, PC9 pos, ...
  annot_df <- annot_df[order(pc_num, dirn)]

  message("  Found annotations (ordered):")
  message("  ", paste(annot_df$annot, collapse = ", "))

  # Baseline prefix and common LDSC inputs
  baseline_prefix <- file.path(tools_dir, "baselineLD_v2.3/baselineLD.")
  weights_prefix  <- file.path(tools_dir, "1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.")
  frq_prefix      <- file.path(tools_dir, "1000G_Phase3_frq/1000G.EUR.QC.")

  # Results base dir for this celltype
  results_base_dir <- file.path(base_outdir, "ldsc", "results")
  dir.create(results_base_dir, showWarnings = FALSE, recursive = TRUE)

  # For each annotation (PCxdirn), create a single-annot panel and run LDSC
  for (i in seq_len(nrow(annot_df))) {
    ann_name <- annot_df$annot[i]
    pc_lab   <- annot_df$pc_label[i]
    dirn     <- as.character(annot_df$dirn[i])

    message("  --- Running ", ann_name, " (", pc_lab, " ", dirn, ") ---")

    # Results directory for this particular annotation
    this_res_dir <- file.path(results_base_dir, ann_name)
    dir.create(this_res_dir, showWarnings = FALSE, recursive = TRUE)

    ## 1) Check which .results files are present for this annotation
    gwas_ids <- names(ss)  # e.g. "Metabolic_UKB-HbA1c", ...
    expected_results <- file.path(this_res_dir, paste0(gwas_ids, ".results"))

    done_vec <- file.exists(expected_results)

    if (!force && all(done_vec)) {
      # All GWAS done for this PCxdirection, skip to next
      message("  All .results present for ", ann_name, ", skipping.")
      next
    }

    # Subset to only missing GWAS unless force == TRUE
    if (force) {
      missing_gwas_ids <- gwas_ids
    } else {
      missing_gwas_ids <- gwas_ids[!done_vec]
    }
    ss_missing <- ss[missing_gwas_ids]

    message("  Missing results for ",
            length(missing_gwas_ids), " / ", length(gwas_ids), " GWAS: ",
            paste(missing_gwas_ids, collapse = ", "))

    ## 2) Only now prepare the single-annotation LDSC panel
    single_annot_dir <- file.path(full_annot_dir, ann_name)
    dir.create(single_annot_dir, showWarnings = FALSE, recursive = TRUE)

    subset_annots_gz(
      in_dir      = full_annot_dir,
      out_dir     = single_annot_dir,
      keep_annots = ann_name,
      prefix      = "1000G.EUR.QC.",
      chrs        = 1:22
    )

    custom_prefix <- file.path(single_annot_dir, "1000G.EUR.QC.")

    ## 3) Run LDSC only for the chosen GWAS
    run_partitioned_h2(
      sumstats           = ss_missing,
      ref_ld_prefixes    = c(baseline_prefix, custom_prefix),
      weights_prefix     = weights_prefix,
      frq_prefix         = frq_prefix,
      out_dir            = this_res_dir,
      print_coefficients = TRUE,
      overlap_annot      = TRUE
    )

    message("  Finished LDSC --h2 for ", ann_name,
            " (results in: ", this_res_dir, ")")
  }


  message("=== Done with ", celltype, " ===")
}

## 5) Collect the results: the custom annotation row (last row) of each .results file,
## saved per cell type as ldsc/results/<celltype>_results.rds (read by figures/fig3)
for (celltype in celltypes) {
  results_base_dir <- file.path(data_dir, "scATACseq/correlations", celltype, "ldsc", "results")
  results_files <- list.files(results_base_dir, pattern = "\\.results$", recursive = TRUE, full.names = TRUE)
  if (!length(results_files)) next

  collected <- rbindlist(lapply(results_files, function(f) {
    res_row <- tail(fread(f), 1)
    annotation <- basename(dirname(f))   # e.g. beta_PC4_combined
    res_row[, `:=`(
      celltype   = celltype,
      annotation = annotation,
      PC         = str_match(annotation, "_(PC[0-9]+)_(positive|negative|combined)$")[, 2],
      subset     = str_match(annotation, "_(PC[0-9]+)_(positive|negative|combined)$")[, 3],
      gwas       = sub("\\.results$", "", basename(f))
    )]
  }))
  saveRDS(as.data.frame(collected), file.path(results_base_dir, paste0(celltype, "_results.rds")))
}
