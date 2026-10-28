#!/bin/bash
# ==============================================================================
# align_combinatorial.sh
# Purpose : Align combinatorial-indexing snATAC-seq with STAR and create fragment
#           files with sinto. The cell barcode is part of the read name, so no
#           whitelist is used; sinto extracts the barcode from the read name.
# Inputs  : <raw_dir>: one subfolder per sample containing one .sra file;
#           data/annotation/GRCh38_index (STAR index, as in align_10x_atac.sh)
# Outputs : <out_dir>/<sample>/: Aligned.sortedByCoord.out.bam(.bai), Log.final.out,
#           frags.sort.bed.gz(.tbi)
# Usage:
#   bash pipeline/01_alignment/align_combinatorial.sh <raw_dir> <out_dir> [threads]
#   threads defaults to 24
#
# Note: the FASTQ files are written to (and left in) the raw sample folders.
#
# Tools: sra-tools (fasterq-dump), STAR 2.7.10a, samtools 1.17, sinto, htslib 1.21.
# ==============================================================================
set -euo pipefail

genome_dir="data/annotation/GRCh38_index"

[[ $# -ge 2 ]] || { echo "Usage: bash pipeline/01_alignment/align_combinatorial.sh <raw_dir> <out_dir> [threads]" >&2; exit 1; }

raw_dir=$(realpath "$1")
out_dir=$(realpath "$2")
THREADS="${3:-24}"
genome_dir=$(realpath "$genome_dir")

cd "$raw_dir" || exit 1
VAR=$(ls)
for i in $VAR
do
    mkdir -p "$out_dir/$i"
    cd "$raw_dir/$i" || exit 1
    # One .sra file per sample folder
    n_sra=$(ls *.sra 2>/dev/null | wc -l || true)
    if [[ $n_sra -ne 1 ]]; then
        echo "Sample $i: expected one .sra file, found $n_sra" >&2
        exit 1
    fi
    # SRA to FASTQ; '@$sn' sets the read name to the SRA spot name, which holds the
    # combinatorial cell barcode (barcode before the first colon)
    fasterq-dump *.sra --seq-defline '@$sn' --qual-defline '+'
    gzip *.fastq
    # One FASTQ pair expected, so leftovers of an earlier run are not aligned
    n_fq1=$(ls *1.fastq.gz 2>/dev/null | wc -l || true)
    n_fq2=$(ls *2.fastq.gz 2>/dev/null | wc -l || true)
    if [[ $n_fq1 -ne 1 ]] || [[ $n_fq2 -ne 1 ]]; then
        echo "Sample $i: expected one *1.fastq.gz and one *2.fastq.gz file, found $n_fq1 and $n_fq2" >&2
        exit 1
    fi
    # Align (plain STAR, no STARsolo: the barcode stays in the read name; options as in
    # align_10x_atac.sh, without clipping)
    #   --alignIntronMax 1, --alignSJDBoverhangMin 999: no spliced alignments (ATAC)
    #   --alignMatesGapMax 1000: maximum fragment gap between mates
    #   --outMultimapperOrder Random --outSAMmultNmax 1: one random alignment per
    #   multimapper; --runRNGseed 777 is STAR's default, written out for reproducibility
    STAR --genomeDir "$genome_dir" \
        --readFilesIn *1.fastq.gz *2.fastq.gz --readFilesCommand zcat \
        --outSAMattributes NH HI nM AS GX GN \
        --outSAMtype BAM SortedByCoordinate \
        --runThreadN "$THREADS" \
        --outMultimapperOrder Random --outSAMmultNmax 1 --runRNGseed 777 \
        --alignIntronMax 1 --alignMatesGapMax 1000 \
        --outFilterScoreMin 30 \
        --alignSJDBoverhangMin 999 \
        --outFileNamePrefix "$out_dir/$i/"
    cd "$out_dir/$i" || exit 1
    # Index BAM
    samtools index -@ "$THREADS" Aligned.sortedByCoord.out.bam
    # Create fragment file, barcode = read name up to the first colon
    sinto fragments -b Aligned.sortedByCoord.out.bam -p "$THREADS" -f fragments.bed --barcode_regex "[^:]*"
    # Sort, compress and index fragment file (sort uses up to 30 GB memory)
    sort -k 1,1 -k2,2n --parallel "$THREADS" -S 30G fragments.bed > frags.sort.bed
    bgzip -@ "$THREADS" frags.sort.bed
    tabix -p bed frags.sort.bed.gz
    rm -f fragments.bed
    cd "$raw_dir" || exit 1
done
