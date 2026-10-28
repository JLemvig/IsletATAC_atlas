#!/bin/bash
# ==============================================================================
# align_10x_atac.sh
# Purpose : Align 10x snATAC-seq (and the ATAC part of 10x multiome) with STARsolo
#           and create fragment files with sinto.
# Inputs  : <raw_dir>: one subfolder per sample with gzipped FASTQ files named
#           *R1*, *R2*, *R3* (R1 = read 1, R2 = cell barcode, R3 = read 2; lanes are
#           concatenated); data/annotation/GRCh38_index (STAR index);
#           data/whitelist/ (barcode whitelists, see below)
# Outputs : <out_dir>/<sample>/: Aligned.sortedByCoord.out.bam(.bai), Solo.out/,
#           Log.final.out, frags.sort.bed.gz(.tbi);
#           <out_dir>/finallogout.txt (STAR log summary) and <out_dir>/barcodes.txt
#           (barcode stats), used for Supplementary Table S2
# Usage:
#   bash pipeline/01_alignment/align_10x_atac.sh <raw_dir> <out_dir> [threads]
#   threads defaults to 24
#
# Per sample, the barcode whitelist is chosen as the one with most hits among the
# first 1M barcode reads (10x ATAC forward/reverse, multiome ATAC
# forward/reverse), and both mates are clipped by 48 bases if reads are longer
# than 97 bases.
#
# Whitelists in data/whitelist/:
#   whitelist_atac            10x scATAC v1 barcodes (737K-cratac-v1.txt)
#   whitelist_atac_rev        reverse complement of whitelist_atac
#   737K-arc-v1-ATAC.txt      10x multiome ATAC barcodes (cellranger-arc)
#   737K-arc-v1-revATAC.txt   reverse complement of 737K-arc-v1-ATAC.txt
#
# Tools: STAR 2.7.10a, samtools 1.17, sinto, htslib 1.21 (bgzip, tabix).
# ==============================================================================
set -euo pipefail

genome_dir="data/annotation/GRCh38_index"
whitelist_dir="data/whitelist"

[[ $# -ge 2 ]] || { echo "Usage: bash pipeline/01_alignment/align_10x_atac.sh <raw_dir> <out_dir> [threads]" >&2; exit 1; }

raw_dir=$(realpath "$1")
out_dir=$(realpath "$2")
THREADS="${3:-24}"
genome_dir=$(realpath "$genome_dir")
whitelist_dir=$(realpath "$whitelist_dir")

# Start the per-run log summaries empty, so a rerun does not append duplicate rows
: > "$out_dir/finallogout.txt"
: > "$out_dir/barcodes.txt"

cd "$raw_dir" || exit 1
VAR=$(ls)
for i in $VAR
do
    mkdir -p "$out_dir/$i"
    cd "$raw_dir/$i" || exit 1
    # Concatenate lanes (remove leftovers of an aborted run first, so they are not
    # matched by the globs below)
    rm -f R1.fq.gz R2.fq.gz R3.fq.gz
    cat *R1* > R1.fq.gz
    cat *R2* > R2.fq.gz
    cat *R3* > R3.fq.gz
    # The read-count pipelines below stop reading early (head) and grep exits 1 on
    # zero hits; "|| true" keeps set -o pipefail from aborting, the values are unchanged.
    # Determine read length from the first read, if above 97, set CUT to 48
    CUT=$(gunzip -c R1.fq.gz | head -n 4 | awk 'NR==2{if(length>97) print 48; else print 0}' || true)
    # Barcode starts at base 1 of the barcode read (base 9 for 737K-arc-v1-revATAC.txt, see below)
    CBSTART=1
    # Count barcode hits in forward and reverse whitelists, using the first 1M reads
    # (4,000,000 FASTQ lines)
    revcount=$(gunzip -c R2.fq.gz | head -n 4000000 | awk 'NR%4==2 { print $0 }' | grep -f "$whitelist_dir/whitelist_atac_rev" - | wc -l || true)
    fwdcount=$(gunzip -c R2.fq.gz | head -n 4000000 | awk 'NR%4==2 { print $0 }' | grep -f "$whitelist_dir/whitelist_atac" - | wc -l || true)
    arcataccount=$(gunzip -c R2.fq.gz | head -n 4000000 | awk 'NR%4==2 { print $0 }' | grep -f "$whitelist_dir/737K-arc-v1-ATAC.txt" - | wc -l || true)
    arcatacrevcount=$(gunzip -c R2.fq.gz | head -n 4000000 | awk 'NR%4==2 { print $0 }' | grep -f "$whitelist_dir/737K-arc-v1-revATAC.txt" - | wc -l || true)
    # Select whitelist based on highest counts
    if [[ $revcount -gt $fwdcount ]] && [[ $revcount -gt $arcataccount ]] && [[ $revcount -gt $arcatacrevcount ]];
    then
        WL="$whitelist_dir/whitelist_atac_rev"
    elif [[ $fwdcount -gt $revcount ]] && [[ $fwdcount -gt $arcataccount ]] && [[ $fwdcount -gt $arcatacrevcount ]];
    then
        WL="$whitelist_dir/whitelist_atac"
    elif [[ $arcataccount -gt $revcount ]] && [[ $arcataccount -gt $fwdcount ]] && [[ $arcataccount -gt $arcatacrevcount ]];
    then
        WL="$whitelist_dir/737K-arc-v1-ATAC.txt"
    else
        WL="$whitelist_dir/737K-arc-v1-revATAC.txt"
        CBSTART=9
    fi
    # Align using STARsolo
    #   --soloType CB_samTagOut: only add the corrected cell barcode (CB tag) to the BAM
    #   --soloCellFilter None: all barcodes are kept, cells are filtered later in ArchR
    #   --alignIntronMax 1, --alignSJDBoverhangMin 999: no spliced alignments (ATAC)
    #   --alignMatesGapMax 1000: maximum fragment gap between mates
    #   --outMultimapperOrder Random --outSAMmultNmax 1: one random alignment per
    #   multimapper; --runRNGseed 777 is STAR's default, written out for reproducibility
    STAR --genomeDir "$genome_dir" \
        --readFilesIn R1.fq.gz R3.fq.gz R2.fq.gz --readFilesCommand zcat \
        --soloType CB_samTagOut --soloCBmatchWLtype 1MM --soloCBwhitelist "$WL" \
        --outSAMattributes NH HI nM AS CR UR CB GX GN sS sQ sM CY \
        --outSAMtype BAM SortedByCoordinate \
        --soloBarcodeReadLength 0 \
        --runThreadN "$THREADS" \
        --outMultimapperOrder Random --outSAMmultNmax 1 --runRNGseed 777 \
        --alignIntronMax 1 --alignMatesGapMax 1000 \
        --soloCellFilter None \
        --outFilterScoreMin 30 \
        --alignSJDBoverhangMin 999 \
        --clip3pNbases "$CUT" "$CUT" \
        --soloCBstart "$CBSTART" \
        --outFileNamePrefix "$out_dir/$i/"
    rm R1.fq.gz R2.fq.gz R3.fq.gz
    cd "$out_dir/$i" || exit 1
    # Index BAM
    samtools index -@ "$THREADS" Aligned.sortedByCoord.out.bam
    # Create fragment file, cell barcode from the CB tag
    sinto fragments -b Aligned.sortedByCoord.out.bam -p "$THREADS" -f fragments.bed -t CB
    # Sort, compress and index fragment file (sort uses up to 30 GB memory)
    sort -k 1,1 -k2,2n --parallel "$THREADS" -S 30G fragments.bed > frags.sort.bed
    bgzip -@ "$THREADS" frags.sort.bed
    tabix -p bed frags.sort.bed.gz
    rm -f fragments.bed
    # Combine logs
    printf '\n%s' "$i" >> "$out_dir/finallogout.txt"
    awk -F '|' '{print $2}' Log.final.out | tr '\n' ' ' >> "$out_dir/finallogout.txt"
    cd Solo.out/ || exit 1
    printf '\n%s' "$i" >> "$out_dir/barcodes.txt"
    sed -n 's/[^0-9]*\([0-9]\+\).*/\1/p' Barcodes.stats | tr '\n' ' ' >> "$out_dir/barcodes.txt"
    cd "$raw_dir" || exit 1
done
