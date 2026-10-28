# 10_chrombpnet

ChromBPNet models of per-cell-type pseudo-bulk accessibility, fine-tuned on non-diabetic (ND)
and type 2 diabetic (T2D) donors, and motif ablation (Methods; Fig. 5, ESM Fig. 5).
Notebooks are run from the repository root. Paths are built from `data/public`, `data/project`
and `data/project/chrombpnet/` (`chrombpnet_dir`), which mirrors the working folder of the original
analysis: `chrombpnet_dir/data/` holds the reference files, bigWigs, `fold<k>/` model folders and
`peaks/`. `figures/fig5/Figure5_clean.ipynb` reads the same tree with `data/...` relative to
`chrombpnet_dir`.

## Run order

| Step | File | Content |
|---|---|---|
| 1 | `setup_chrombpnet_inputs.sh` | Zenodo folds, ENCODE blacklist, GRCh38 genome, motif collection, blacklist-filtered peak set with 250 bp summits (`peaks_no_blacklist_summits.bed`), copy of the PC-associated peak BEDs (alpha, beta) |
| 2 | `alpha_beta_train_heldout_metadata.R` | Train / held-out cell lists for alpha and beta cells (was `Subtypes.R`, lines 1-91 and 94-165, without the xgboost subtype model) |
| 3 | `train_celltype_models.ipynb` (was `Preprocess.ipynb`) | Pseudo-bulk bigWigs per cell type, bigWig QC, GC-matched negatives, bias model and per-cell-type models per fold (5 folds, seed 1234); `training_metrics.csv` (ESM Fig. 5b) |
| 4 | `prepare_finetune_peaks.ipynb` (cell moved from `Interpretation All Cell types.ipynb`) | `{alpha,beta}_correlated.bed` from the first-run PC peak BEDs (see "Peak sets") |
| 5 | `finetune_nd_t2d_models.ipynb` (was `Finetune.ipynb`) | Alpha and beta models fine-tuned on ND or T2D donors; evaluation on `{alpha,beta}_correlated.bed`; `finetune_metrics.csv`, `finetune_eval.pdf`, `{alpha,beta}_model_evaluation.pdf` (ESM Fig. 5b-c). Models are read by `figures/fig5/Figure5_clean.ipynb` |
| 6 | `motif_ablation_alpha_beta.ipynb` (was `Interpretation.ipynb`, latest version, only the cells that produce files) | FIMO scan of all peaks and PC-associated peaks, motif enrichment (`motif_enrichments.pkl`) and motif ablation on the ND/T2D models (`peaks/correlated/old/{midx}.pkl`, `motif_ablations.pkl`). The same code is repeated inline in `figures/fig5/Figure5_clean.ipynb`, which skips motifs whose `old/{midx}.pkl` exists |
| 7 | `motif_ablation_all_celltypes.ipynb` (was `Interpretation All Cell types.ipynb`) | Models for all cell types (`celltype_models/{celltype}.pt`), genome-wide motif scan (`all_peaks.motifs`), motif ablation (`motif_activity.pkl`), `MarginPeaks_top1000.csv`; the last cells make the PBX1 panels of Fig. 5 (`PBX1_protein.pdf`, `PBX1_motif_GWAS.pdf`) |
| 8 | `motif_activity_heatmap.ipynb` (was `Motif activity.ipynb`) | Motif by cell type heatmap from `motif_activity.pkl` |

## Peak sets (what each step used)

- Training (all cell types): `peaks_no_blacklist_summits.bed` = `macs2_peaks.narrowPeak` of the final
  ArchR project without ENCODE blacklist regions (ENCFF356LFX, extended by 1057 bp); per cell type
  `{celltype}_filtered.bed` is derived from it in step 3.
- Fine-tuning evaluation (step 5): `{alpha,beta}_correlated.bed`, built in step 4 from the **first**
  correlation run (not the DESeq2 LRT peak sets of the paper): peaks associated with PC4, PC9, PC10,
  PC11 and PC13 (`pc<k>.bed`, correlation of DESeq2-normalised pseudo-bulk counts with the PC, adjusted
  for PC1, PC2, PC3, PC5, PC6, PC7 and PC18) from
  `data/public/scATACseq/correlations/<celltype>/pc1pc2pc3pc5pc6pc7pc18/DESeq_normalized_nogc/bed/`.
  These BEDs have five columns (chromosome, start, end, peak id, GC). They were written by
  `120_IsletEpigenetics_notes.R` ("To Jesper" block, `create_bed`) from the `pc<k>.csv` files of
  `122_IsletEpigenetics_create_correlation.R` (internal scripts, not in this repository), whose inputs
  (`peakset.rds`, `peakmatrix_DESeq_normalized_nogc.csv`, `open_peak_fraction.csv`) come from
  `121_IsletEpigenetics_peaknormalization.R`. The code is published as it was run: this upstream
  chain is not in the repository and the inputs of that run are no longer available, so step 4 cannot
  be reproduced from the scripts here; the published `{alpha,beta}_correlated.bed` are the files used.
- Motif scan, enrichment and ablation (step 6 and Fig. 5): `{celltype}_filtered.bed` plus the
  PC-associated peaks `PC<k>_{positive,negative}.bed` of the DESeq2 LRT run
  (`pipeline/08_gwas/ldsc_pc_peaks.R`; `*_combined.bed` is not used). The files used for Fig. 5 were
  written on 2025-11-28; 38 of the 46 alpha and beta files are identical to the final ones from
  2025-12-03 and 8 differ by one or two peaks (alpha PC4, PC11, PC12, PC16 positive; beta PC4 and PC11
  positive, PC11 and PC14 negative). `setup_chrombpnet_inputs.sh` copies the final files.
- Ablation results read by `Figure5_clean.ipynb`: `data/peaks/correlated/old/{midx}.pkl`, 637 files
  written 2025-12-08 to 2025-12-10 by the ablation cell of step 6 (then writing to
  `peaks/correlated/{midx}.pkl`; the folder was renamed `old/` later).

## Versions and downloads

- Training and fine-tuning (steps 3 and 5): chrombpnet 1.0.1, Python 3.8, TensorFlow 2.8.0 (conda-forge,
  cudatoolkit 11.2, cudnn 8.1), samtools, bedtools, ucsc-bedgraphtobigwig, pybigwig, meme, pysam, seaborn.
  `bedtools` and `bedGraphToBigWig` must be on the PATH.
- Interpretation (steps 6-8): Python 3.12.3, PyTorch, bpnet-lite 1.0.0 (Figure5_clean) / 1.1.0
  (latest Interpretation), tangermeme 1.5.0 / 1.6.0, memelite 0.5.0, numba 0.60.0, numpy 2.0.2, pandas,
  pysam 0.24.1, seaborn 0.13.2, statsmodels 0.15.0, scipy, tqdm, h5py. Step 7 was last run with
  bpnet-lite 1.0.0, tangermeme 1.0.3, memelite 0.2.0 (pyranges 0.1.4, liftover 1.3.3, xgboost 3.2.0,
  statsmodels 0.14.6, pysam 0.23.3); further packages used by
  step 7: anndata, scanpy, adpbulk, edgepython, pyranges.
  The versions of the December 2025 ablation run were not recorded.
- Zenodo 7445373 (`fold_{0..4}.json`), ENCODE ENCFF356LFX (blacklist), GRCh38_EBV.chrom.sizes,
  GRCh38 no-alt analysis set (GCA_000001405.15), motif-clustering v2.1beta consensus PWMs and
  `metadata.tsv` (https://resources.altius.org/~jvierstra/projects/motif-clustering-v2.1beta/): all
  downloaded by `setup_chrombpnet_inputs.sh`.
- Wang et al. snRNA-seq (GSE200044), used in the RNA sections of `figures/fig4/Figure4_clean.ipynb`
  and `figures/fig5/Figure5_clean.ipynb`:
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE200nnn/GSE200044/suppl/GSE200044%5Fcell%5Fcluster.tar.gz,
  extracted to `data/project/chrombpnet/cell_cluster/`.

## Inputs without a producer in the repository (TODO)

- `data/public/scATACseq/ArchR/macs2_peaks.narrowPeak`: the peak set of the final ArchR project's
  PeakMatrix (456,191 peaks, `addReproduciblePeakSet` in `pipeline/03_doublets/041`) exported as
  narrowPeak; the export line itself is not preserved.
- First-run PC peak BEDs `pc<k>.bed` (step 4): `120_IsletEpigenetics_notes.R`,
  `122_IsletEpigenetics_create_correlation.R` and `121_IsletEpigenetics_peaknormalization.R` are not
  in the repository, and the files `peakset.rds` / `peakmatrix_DESeq_normalized_nogc.csv` of that run
  are no longer available (the current versions of 121 and 122 are from 2026 and read the final
  project).
- `IC_Final_Raw_AnnData.h5ad` (step 7, peak count AnnData of all cells): `ExportArchR.R` and
  `CreateH5.ipynb` (internal scripts, not in this repository).
- `donor_metadata.tsv` and `cell_cluster/multiome_RNA*.h5ad` (RNA sections of the Fig. 4 and Fig. 5
  notebooks), `EnhancerPredictionsAllPutative.tsv` (ABC calls combined over cell types,
  `pipeline/09_abc`; Fig. 4 and Fig. 5 notebooks), `GWAs/Suzuki.ALL-T2D.sumstats` and `PBX1.protein.csv` (PBX1 panels of Fig. 5, made in step 7:
  `PBX1_protein.pdf`, `PBX1_motif_GWAS.pdf`).
- ESM Fig. 5 plotting code; the Dec 2025 pickles under `old/` also have a JASPAR
  variant (`{midx}_JASPAR.pkl`) that no notebook in this repository reads.
