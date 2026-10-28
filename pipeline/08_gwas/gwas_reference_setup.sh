#!/bin/bash
# ==============================================================================
# gwas_reference_setup.sh
# Purpose : Install LDSC and TORUS and download the LDSC reference files and GWAS
#           summary statistics used by ldsc_pc_peaks.R and torus_pc_peaks.R.
# Inputs  : env/ldsc.yml (conda env 'ldsc', python2);
#           data/project/gwas_environment/1000G_Phase3_baselineLD_v2.3_ldscores.zip
# Outputs : data/project/gwas_environment/:
#             miniconda3/ (with env 'ldsc'), ldsc/ (ldsc.py),
#             torus/torus-master/src/torus,
#             1000G_EUR_Phase3_plink/, 1000G_Phase3_frq/,
#             1000G_Phase3_weights_hm3_no_MHC/, hm3_no_MHC.list.txt,
#             baselineLD_v2.3/, sumstats/ (UKB_460K.*.sumstats.gz)
# Usage:
#   bash pipeline/08_gwas/gwas_reference_setup.sh
#
# TORUS needs the GSL, Boost, zlib and liblzma development libraries, e.g. on Ubuntu:
#   apt-get install liblzma-dev libgsl-dev libcurl4-openssl-dev libboost-all-dev
# Suzuki.ALL-T2D.sumstats.gz is not downloaded here, see README.md (External and
# manual steps).
#
# Tools: Miniconda3, LDSC (commit aa33296), TORUS (master branch), wget, unzip.
# ==============================================================================
set -euo pipefail

tools_dir="data/project/gwas_environment"
ldsc_yml=$(realpath "env/ldsc.yml")

mkdir -p "$tools_dir"
cd "$tools_dir" || exit 1

## LDSC: conda env 'ldsc' in tools_dir/miniconda3 (path used by ldsc_pc_peaks.R)
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b -u -p "$PWD/miniconda3"
rm Miniconda3-latest-Linux-x86_64.sh
source miniconda3/etc/profile.d/conda.sh
conda env create -f "$ldsc_yml"
git clone https://github.com/bulik/ldsc.git
git -C ldsc checkout aa33296abac9569a6422ee6ba7eb4b902422cc74

## TORUS
mkdir -p torus
wget https://github.com/xqwen/torus/archive/refs/heads/master.zip -O torus/master.zip
unzip -q torus/master.zip -d torus
make -C torus/torus-master/src

## LDSC reference files and UK Biobank summary statistics (https://zenodo.org/records/10515792)
for f in hm3_no_MHC.list.txt 1000G_Phase3_plinkfiles.tgz 1000G_Phase3_frq.tgz \
         1000G_Phase3_weights_hm3_no_MHC.tgz sumstats.tgz; do
    wget "https://zenodo.org/api/records/10515792/files/$f/content" -O "$f"
done
for f in 1000G_Phase3_plinkfiles.tgz 1000G_Phase3_frq.tgz \
         1000G_Phase3_weights_hm3_no_MHC.tgz sumstats.tgz; do
    tar xzf "$f"
done

## baselineLD v2.3 LD scores (1000G Phase 3, hg19), unpacks to baselineLD_v2.3/
# 1000G_Phase3_baselineLD_v2.3_ldscores.zip is not downloaded here; place it in tools_dir first
unzip -q 1000G_Phase3_baselineLD_v2.3_ldscores.zip -x '__MACOSX/*'
