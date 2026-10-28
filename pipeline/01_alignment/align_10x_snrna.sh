#!/bin/bash
# ==============================================================================
# align_10x_snrna.sh
# Purpose : Download and align 10x snRNA-seq (RNA part of the Wang et al. multiome,
#           PRJNA822928) with STARsolo. Used for the snRNA comparison of the cycling
#           trajectories (Figure 2B) and gene expression.
# Inputs  : <download_manifest>: tab-separated file, one line per FASTQ file:
#           donor, file name, URL, md5 (one R1.fq.gz / R2.fq.gz pair per donor), e.g.
#             <donor>  <donor>_R1.fq.gz  <url>  <md5>
#           data/annotation/GRCh38_gencode35_index (STAR index, GENCODE v35);
#           data/whitelist/737K-arc-v1.txt (10x multiome gene expression barcodes)
# Outputs : <out_dir>/<donor>/: Solo.out/ (Gene, GeneFull, GeneFull_Ex50pAS;
#           multimappers by EM), STAR logs, failed.downloads (empty on success)
# Usage:
#   bash pipeline/01_alignment/align_10x_snrna.sh <download_manifest> <out_dir> [threads]
#   threads defaults to 60
#
# Tools: STAR 2.7.10a (GENCODE v35 index), aria2.
# ==============================================================================
set -euo pipefail

genome_dir="data/annotation/GRCh38_gencode35_index"
whitelist="data/whitelist/737K-arc-v1.txt"

[[ $# -ge 2 ]] || { echo "Usage: bash pipeline/01_alignment/align_10x_snrna.sh <download_manifest> <out_dir> [threads]" >&2; exit 1; }

manifest=$(realpath "$1")
out_dir=$(realpath "$2")
THREADS="${3:-60}"
genome_dir=$(realpath "$genome_dir")
whitelist=$(realpath "$whitelist")

mkdir -p "$out_dir"
work_dir=$(mktemp -d)
cd "$work_dir" || exit 1
# Unique donors, sorted
Donors=$(cut -f 1 "$manifest" | sort | uniq)

# Load STAR genome into shared memory, used by all donors (LoadAndKeep below)
STAR --genomeDir "$genome_dir" --genomeLoad LoadAndExit
rm -f Log* Aligned.out.sam SJ.out.tab
# Remove the genome from shared memory also if the script stops early
trap 'STAR --genomeDir "$genome_dir" --genomeLoad Remove >/dev/null 2>&1 || true' EXIT

# Loop over donors
for z in $Donors; do
    echo "Processing donor: $z"

    # Files to download for this donor
    awk -v VAR="$z" '$1 == VAR { print $0 }' "$manifest" > Donor
    awk '{ print $3"\n out="$2"\n checksum=md5="$4 }' Donor > Download

    # Download using aria2, rerun (at most 5 times) until all files pass the checksum.
    # --save-session lists the unfinished downloads; an empty file means success.
    # "|| true": aria2c exits non-zero on failed downloads, which are retried below.
    aria2c -i Download -j10 --check-integrity=true --save-session failed.downloads || true
    has_error=$(wc -l < failed.downloads)
    retries=0
    while [ "$has_error" -gt 0 ]; do
        if [ "$retries" -ge 5 ]; then
            echo "Downloads for donor $z still fail after 5 retries, see failed.downloads" >&2
            exit 1
        fi
        retries=$((retries + 1))
        echo "Still has $has_error errors, rerun aria2 to download ..."
        mv failed.downloads Download
        aria2c -i Download -j10 --check-integrity=true -c --save-session failed.downloads || true
        has_error=$(wc -l < failed.downloads)
        sleep 10
    done

    # One R1/R2 pair per donor: STAR takes one file per mate here
    n_r1=$(ls *R1.fq.gz 2>/dev/null | wc -l || true)
    n_r2=$(ls *R2.fq.gz 2>/dev/null | wc -l || true)
    if [[ $n_r1 -ne 1 ]] || [[ $n_r2 -ne 1 ]]; then
        echo "Donor $z: expected one R1.fq.gz and one R2.fq.gz file, found $n_r1 and $n_r2" >&2
        exit 1
    fi

    # Align using STARsolo (R2 = cDNA, R1 = cell barcode + UMI)
    #   --soloFeatures: Gene (exonic reads), GeneFull (exonic + intronic reads),
    #   GeneFull_Ex50pAS (as GeneFull, reads with >50% exonic overlap prioritised,
    #   fully exonic antisense reads not counted)
    #   --soloCellFilter None: all barcodes are kept, cells are filtered downstream
    #   --soloMultiMappers EM: multimapping reads are distributed by EM
    #   --soloUMIlen 12: 16 bp cell barcode + 12 bp UMI (10x multiome GEX)
    #   --outMultimapperOrder Random --outSAMmultNmax 1: one random alignment per
    #   multimapper; --runRNGseed 777 is STAR's default, written out for reproducibility
    STAR --genomeLoad LoadAndKeep --genomeDir "$genome_dir" \
        --readFilesIn *R2.fq.gz *R1.fq.gz --readFilesCommand zcat \
        --soloType CB_UMI_Simple \
        --soloFeatures Gene GeneFull GeneFull_Ex50pAS \
        --soloCellFilter None \
        --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
        --clipAdapterType CellRanger4 \
        --outFilterScoreMin 30 \
        --soloMultiMappers EM \
        --soloUMIfiltering MultiGeneUMI_CR --soloUMIdedup 1MM_CR \
        --runThreadN "$THREADS" \
        --outMultimapperOrder Random --outSAMmultNmax 1 --runRNGseed 777 \
        --soloCBwhitelist "$whitelist" --soloBarcodeReadLength 0 --soloUMIlen 12

    # Move results to donor folder
    mkdir -p "$out_dir/$z/"
    mv failed.downloads Solo.out Log* "$out_dir/$z/"
    rm -f Aligned.out.sam SJ.out.tab *.fq.gz Donor Download
done

# Unload the genome from shared memory
STAR --genomeDir "$genome_dir" --genomeLoad Remove
trap - EXIT
rm -f Log* Aligned.out.sam SJ.out.tab
rmdir "$work_dir"
