# ==============================================================================
# torus_pc_peaks.R
# Purpose : TORUS enrichment of GWAS signal in the peaks associated with each donor PC
#           (positive, negative and combined BED files, lifted to hg19 and widened to
#           1 kb windows), for metabolic and control GWAS traits (mapgen).
# Inputs  : data/public/scATACseq/correlations/<celltype>/bed/*.bed (ldsc_pc_peaks.R);
#           data/project/gwas_environment (1000G Phase 3 plink and frq files, GWAS
#           summary statistics, compiled TORUS)
# Outputs : data/public/scATACseq/correlations/<celltype>/torus/<celltype>_results.rds and
#           <bed>_torus_fits.rds
# TORUS: https://github.com/xqwen/torus, compiled in tools_dir/torus/torus-master/src/torus
# ==============================================================================

data_dir      <- normalizePath("data/public", mustWork = FALSE)     # public/processed data
tools_dir     <- normalizePath("data/project/gwas_environment", mustWork = FALSE)  # reference files, sumstats, TORUS

library(genio)
library(liftOver)
library(mapgen)
library(GenomicRanges)
library(data.table)
library(rtracklayer)
library(BiocParallel)
library(stringr)

register(MulticoreParam(workers = 4))
options(scipen = 999)

celltypes <- c("beta", "alpha", "acinar", "ductal",
               "delta", "pp", "cycling", "activated_stellate", "quiescent_stellate",
               "endothelial", "myeloid", "lymphoid", "mast", "epsilon", "schwann")

## Prepare for TORUS
LD_Blocks <- readRDS(system.file('extdata', 'LD.blocks.EUR.hg19.rds', package='mapgen'))

# Get SNP infos
for (chr in 1:22) {
  bim <- read_bim(paste(file.path(tools_dir, "1000G_EUR_Phase3_plink/1000G.EUR.QC."), chr, ".bim", sep=""))
  bim <- as.data.frame(bim)
  bim <- bim[,c(1,4,2)]
  colnames(bim) <- c("CHR", "BP","SNP")
  bim <- bim[ order(bim$SNP),]
  freq <- fread(paste(file.path(tools_dir, "1000G_Phase3_frq/1000G.EUR.QC."), chr, ".frq", sep=""), data.table = FALSE)
  freq <- freq[,c("SNP","MAF")]
  freq <- freq[ order(freq$SNP),]
  bim$MAF <- freq$MAF
  if (chr == 1) {
    snp_info <- bim
  } else {
    snp_info <- rbind(snp_info, bim)
  }
}

## Define a list of interesting GWAS
gwas.list <- c(
  file.path(tools_dir, "sumstats/Suzuki.ALL-T2D.sumstats.gz"),
  file.path(tools_dir, "sumstats/UKB_460K.biochemistry_HbA1c.sumstats.gz"),
  file.path(tools_dir, "sumstats/UKB_460K.body_BMIz.sumstats.gz"),
  file.path(tools_dir, "sumstats/UKB_460K.body_WHRadjBMIz.sumstats.gz"),
  file.path(tools_dir, "sumstats/UKB_460K.disease_T2D.sumstats.gz"),
  file.path(tools_dir, "sumstats/UKB_460K.biochemistry_CreactiveProtein.sumstats.gz"),
  file.path(tools_dir, "sumstats/UKB_460K.biochemistry_Glucose.sumstats.gz"),
  file.path(tools_dir, "sumstats/UKB_460K.biochemistry_IGF1.sumstats.gz"),
  file.path(tools_dir, "sumstats/UKB_460K.biochemistry_HDLcholesterol.sumstats.gz"),
  file.path(tools_dir, "sumstats/UKB_460K.biochemistry_LDLdirect.sumstats.gz"),
  file.path(tools_dir, "sumstats/UKB_460K.biochemistry_Triglycerides.sumstats.gz"),
  file.path(tools_dir, "sumstats/UKB_460K.disease_ASTHMA_DIAGNOSED.sumstats.gz"),
  file.path(tools_dir, "sumstats/UKB_460K.body_LEFT_HANDED.sumstats.gz"),
  file.path(tools_dir, "sumstats/UKB_460K.body_BALDING4.sumstats.gz"),
  file.path(tools_dir, "sumstats/UKB_460K.pigment_HAIR.sumstats.gz"),
  file.path(tools_dir, "sumstats/UKB_460K.pigment_SKIN.sumstats.gz"))

gwas.metadata <- data.frame(
  gwas = gsub(".sumstats.gz", "", gsub(file.path(tools_dir, "sumstats/"), "", gwas.list)),
  type = c(rep("Metabolic",11), rep("Control", 5)),
  trait = c("Suzuki-T2D","UKB-HbA1c","UKB-BMI","UKB-WHR","UKB-T2D","UKB-CRP","UKB-Glucose","UKB-IGF1","UKB-HDL","UKB-LDL","UKB-TG","UKB-Asthma", "UKB-Handedness","UKB-Balding","UKB-Hair","UKB-Skin")
)

# Prepare GWAS for Torus analysis
# Calculate SE and beta following this: https://pubmed.ncbi.nlm.nih.gov/27019110/
gwas.list.prepared <- list()
for (gwas in gwas.list[1:length(gwas.list)]) {
  # Import and setup
  dat <- fread(gwas, data.table = FALSE)
  if(startsWith(gwas,file.path(tools_dir, "sumstats/UKB_460K"))){
    dat$Pval <- 1 - pchisq(dat$CHISQ, 1)
    snps_tmp <- snp_info[ snp_info$SNP %in% dat$SNP,]
    snps_tmp <- snps_tmp[ match(dat$SNP, snps_tmp$SNP),]
    dat$MAF <- snps_tmp$MAF
    dat$CHR <- snps_tmp$CHR
    dat$POS <- snps_tmp$BP
    dat$SE <- sqrt((2 * dat$MAF * (1 - dat$MAF) * (dat$N + dat$Z^2)))
    dat$Beta <- dat$Z / dat$SE
    dat$SE <- 1 / dat$SE

    # Prepare
    gwas.list.prepared[[length(gwas.list.prepared)+1]] <- process_gwas_sumstats(
      dat,
      chr = "CHR",
      pos = "POS",
      beta = "Beta",
      se = "SE",
      a0 = "A1",
      a1 = "A2",
      snp = "SNP",
      pval = "Pval",
      remove_indels = TRUE,
      LD_Blocks = LD_Blocks)
    names(gwas.list.prepared)[length(gwas.list.prepared)] <- gsub(".sumstats.gz", "", gsub(file.path(tools_dir, "sumstats/"), "", gwas))
  } else {
    gwas.list.prepared[[length(gwas.list.prepared)+1]] <- process_gwas_sumstats(
      dat,
      chr = "CHR",
      pos = "BP",
      beta = "BETA",
      se = "SE",
      a0 = "A2",
      a1 = "A1",
      snp = "SNP",
      pval = "P",
      remove_indels = TRUE,
      LD_Blocks = LD_Blocks)
    names(gwas.list.prepared)[length(gwas.list.prepared)] <- gsub(".sumstats.gz", "", gsub(file.path(tools_dir, "sumstats/"), "", gwas))
  }
}

# Define function to calculate P-values
pval_from_ci <- function(mean, upper, ci){
  nsamp <- length(mean)
  pval.res <- rep(0, nsamp)
  for(i in 1:nsamp){
    alph <- (1-ci)/2
    zval <- qnorm(p = 1-alph)
    se <- (upper[i]-mean[i])/zval

    pval.res[i] <- 1 - pnorm(q = mean[i] / se)
  }
  return(pval.res)
}

torus_path <- file.path(tools_dir, "torus/torus-master/src/torus")
path <- system.file(package="liftOver", "extdata", "hg38ToHg19.over.chain")
ch <- import.chain(path)

## Run TORUS per cell type and BED file (BED files with fewer than min_peaks windows are
## skipped); the summary table and the full fits are saved
window_size = 1000
min_peaks   = 200

for (celltype in celltypes) {
  message("=== Cell type: ", celltype, " ===")

  base_outdir      <- file.path(file.path(data_dir, "scATACseq/correlations"), celltype)
  bed_dir          <- file.path(base_outdir, "bed")

  if (!dir.exists(bed_dir)) {
    message("  No bed directory, skipping.")
    next
  }

  peaklist <- list.files(bed_dir, full.names = TRUE)
  if (isEmpty(peaklist)) {
    message("  No bed files, skipping.")
    next
  }

  ## output dir for this cell type (used for both summary + per-bed fits)
  results_outdir <- file.path(base_outdir, "torus")
  dir.create(results_outdir, showWarnings = FALSE, recursive = TRUE)

  res <- list()  # summary table rows for this cell type

  for (peak in peaklist) {
    base  <- basename(peak)

    message(str_replace(base, ".bed", ""), " ",
            which(peaklist %in% peak), "/", length(peaklist))

    filename <- peak
    if (file.exists(filename) && file.size(filename) > 0) {

      ## Load peaks and make GRanges object
      bed   <- read.delim(filename, header = FALSE)
      bed_gr <- GRanges(bed$V1, ranges = IRanges(start = bed$V2, end = bed$V3))

      ## Lift from hg38 to hg19
      lo      <- liftOver(bed_gr, ch)
      one2one <- lo[elementNROWS(lo) == 1L]
      res_gr  <- unlist(one2one, use.names = FALSE)
      df_res  <- as.data.frame(res_gr)

      ## Resize peaks
      df_res$center <- as.integer((df_res$start + df_res$end) / 2)
      df_res$start  <- df_res$center - window_size / 2
      df_res$end    <- df_res$center + window_size / 2

      ## Prepare peaks BED for TORUS
      peaks   <- df_res[, 1:3]
      peaks[,1] <- gsub("chr", "", as.character(peaks[,1]))

      if (nrow(peaks) < min_peaks) {
        message("  Only ", nrow(peaks), " peaks (min_peaks = ", min_peaks,
                "), skipping TORUS for ", base)
        next  # skip this bed file
      }

      write.table(
        peaks,
        file.path(tools_dir, "torus.peaks.bed"),
        quote = FALSE, sep = "\t",
        row.names = FALSE, col.names = FALSE
      )

      ## Per-BED list of TORUS fits (one element per GWAS SS)
      torus_fits_for_bed <- list()

      ## Loop across GWAS and run TORUS
      for (idx in seq_along(gwas.list.prepared)) {
        gwas_name <- names(gwas.list.prepared)[idx]

        message(celltype,", ",str_replace(base, ".bed", ""), " ",
                which(peaklist %in% peak), "/", length(peaklist),
                " - ", gwas_name, " ", idx, "/", length(gwas.list.prepared))

        torus.files  <- prepare_torus_input_files(
          gwas.list.prepared[[idx]],
          file.path(tools_dir, "torus.peaks.bed")
        )

        torus.result <- run_torus(
          torus.files$torus_annot_file,
          torus.files$torus_zscore_file,
          option     = "est",
          torus_path = file.path(tools_dir, "torus/torus-master/src/torus")
        )

        ## Store full TORUS result object for this GWAS
        torus_fits_for_bed[[gwas_name]] <- torus.result

        ## Capture enrichment summary
        enrichment <- torus.result$enrich
        enrichment <- as.data.frame(enrichment)
        enrichment <- enrichment[2, ]   # second row = annotation term

        enrichment$term <- str_replace(base, ".bed", "")
        enrichment$Trait <- gwas.metadata[
          gwas.metadata$gwas %in% gwas_name, 3
        ]
        enrichment$Type <- gwas.metadata[
          gwas.metadata$gwas %in% gwas_name, 2
        ]
        enrichment$Celltpe <- basename(dirname(dirname(peak)))
        enrichment$Pval <- pval_from_ci(
          enrichment[1, 2],
          enrichment[1, 4],
          0.95
        )

        res[[length(res) + 1]] <- enrichment
      }

      ## Save per-BED list of TORUS fits
      bed_label <- str_replace(base, ".bed", "")
      saveRDS(
        torus_fits_for_bed,
        file = file.path(results_outdir,
                         paste0(bed_label, "_torus_fits.rds"))
      )
    }
  }

  ## Combine summary across all beds for this cell type
  if (length(res) > 0) {
    torus_results <- do.call(rbind, res)

    ## FDR etc.
    torus_results$`-log10_Pval` <- -log10(torus_results$Pval)
    torus_results$FDR           <- p.adjust(torus_results$Pval, "fdr")
    torus_results$`-log10_FDR`  <- -log10(torus_results$FDR)

    ## Save summary table for this cell type
    saveRDS(
      torus_results,
      file = file.path(results_outdir,
                       paste0(celltype, "_results.rds"))
    )
  } else {
    message("  No TORUS results for cell type ", celltype,
            " (all beds skipped or empty).")
  }
}

