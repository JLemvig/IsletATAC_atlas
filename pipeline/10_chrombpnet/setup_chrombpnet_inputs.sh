#!/bin/bash
# ==============================================================================
# setup_chrombpnet_inputs.sh
# Purpose : Download the reference files for ChromBPNet (chromosome folds, ENCODE
#           blacklist, GRCh38 genome, chromosome sizes, motif collection), build the
#           blacklist-filtered peak set with 250 bp summits, and copy the PC-associated
#           peak BEDs of alpha and beta cells into the ChromBPNet data folder.
# Inputs  : data/public/scATACseq/ArchR/macs2_peaks.narrowPeak (peak set of the final
#           ArchR project, 10-column narrowPeak format);
#           data/public/scATACseq/correlations/<celltype>/bed/PC*_{positive,negative}.bed
#           (alpha and beta; from pipeline/08_gwas/ldsc_pc_peaks.R)
# Outputs : data/project/chrombpnet/data/: fold_{0..4}.json, blacklist.bed,
#           hg38.chrom.sizes, hg38.fa, motifs.meme, motifs_metadata.tsv,
#           peaks_no_blacklist_summits.bed, peaks/correlated/<celltype>/bed/PC*.bed,
#           empty folders peaks/correlated/hits and peaks/correlated/old
# Usage:
#   bash pipeline/10_chrombpnet/setup_chrombpnet_inputs.sh
#
# Tools: wget, gunzip, bedtools, awk.
# ==============================================================================
set -euo pipefail

data_dir=$(realpath "data/public")
out_dir="data/project/chrombpnet/data"
mkdir -p "$out_dir"
out_dir=$(realpath "$out_dir")
cd "$out_dir"

# Peak set of the final ArchR project
cp "$data_dir/scATACseq/ArchR/macs2_peaks.narrowPeak" .

# Chromosome folds (Zenodo 7445373), ENCODE blacklist, chromosome sizes, genome
wget https://zenodo.org/records/7445373/files/fold_0.json
wget https://zenodo.org/records/7445373/files/fold_1.json
wget https://zenodo.org/records/7445373/files/fold_2.json
wget https://zenodo.org/records/7445373/files/fold_3.json
wget https://zenodo.org/records/7445373/files/fold_4.json
wget https://www.encodeproject.org/files/ENCFF356LFX/@@download/ENCFF356LFX.bed.gz -O blacklist.bed.gz
gunzip blacklist.bed.gz
wget https://www.encodeproject.org/files/GRCh38_EBV.chrom.sizes/@@download/GRCh38_EBV.chrom.sizes.tsv -O hg38.chrom.sizes
wget https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15/@@download/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz -O hg38.fa.gz
gunzip hg38.fa.gz

# Motif collection (Vierstra et al. consensus PWMs, v2.1beta) and its annotation
wget https://resources.altius.org/~jvierstra/projects/motif-clustering-v2.1beta/consensus_pwms.meme -O motifs.meme
wget https://resources.altius.org/~jvierstra/projects/motif-clustering-v2.1beta/metadata.tsv -O motifs_metadata.tsv

# Remove blacklisted peaks
bedtools slop -i blacklist.bed -g hg38.chrom.sizes -b 1057 > temp.bed
awk '{
  for(i=1; i<=NF; i++) {
    if(i==2 || i==3)
      printf "%.0f", $i
    else
      printf "%s", $i
    if(i<NF) printf "\t"
    else printf "\n"
  }
}' macs2_peaks.narrowPeak > peaks.narrowPeak
bedtools intersect -v -a peaks.narrowPeak -b temp.bed > peaks_no_blacklist.bed

# Adjust summits
awk '{ print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t250" }' peaks_no_blacklist.bed > peaks_no_blacklist_summits.bed

# PC-associated peaks of alpha and beta cells (positive and negative direction only)
mkdir -p peaks/correlated/hits peaks/correlated/old
for celltype in alpha beta; do
  mkdir -p "peaks/correlated/$celltype/bed"
  cp "$data_dir/scATACseq/correlations/$celltype/bed/"PC*_positive.bed "peaks/correlated/$celltype/bed/"
  cp "$data_dir/scATACseq/correlations/$celltype/bed/"PC*_negative.bed "peaks/correlated/$celltype/bed/"
done
