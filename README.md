# Charting single-nucleus chromatin accessibility in human pancreatic islets

Code for the paper:

Lemvig J, Ernst IVS, Madsen JGS. *Charting single-nucleus chromatin accessibility in human pancreatic islets*. Submitted to Diabetologia (the reference and DOI will be added on publication).

## Study overview

The cis-regulatory landscape of human islets remains incompletely characterised. To build an integrated chromatin accessibility atlas of human islets, we reprocessed and integrated public single-nucleus ATAC-seq (snATAC-seq) data from eight studies. The atlas comprises 526,082 nuclei from 98 donors and resolves 15 cell populations, including rare epsilon and Schwann cells. The atlas can be browsed at https://www.isletcartography.com/.

Cycling cells followed trajectories towards all four major endocrine cell types, and their proportions were associated with age, BMI, sex and type 2 diabetes. A semi-supervised variational autoencoder (PeakANVI) learned separate cell and donor embeddings, and the donor components were tested for association with chromatin accessibility and donor traits. Peaks associated with the type 2 diabetes components were linked to target genes with the activity-by-contact (ABC) model and tested for enrichment of GWAS variants (LDSC, TORUS). Candidate transcription factors were prioritised by motif enrichment, ChromBPNet models and gene expression, which nominated PBX1 as a candidate regulator of beta cell dysfunction in type 2 diabetes.

<!-- TODO: illustration, e.g. the graphical abstract as illustrations/graphical_abstract.png, with a short italic caption below it -->

# Brief description of folders

- `pipeline/`: analysis steps in run order, one folder per Methods section
- `figures/`: scripts and notebooks that make the paper figures (`figN` = main figure, `esm_figN` = ESM figure)
- `env/`: environment files
- `data/`, `results/`: not tracked, see "Data"

All scripts are run from the repository root.

# Description of scripts

### `pipeline/01_alignment/`: alignment

- `align_10x_atac.sh`: align 10x snATAC-seq (and the ATAC part of 10x multiome) with STARsolo, fragment files with sinto
- `align_combinatorial.sh`: align combinatorial-indexing snATAC-seq with STAR, fragment files with sinto
- `align_10x_snrna.sh`: download and align the Wang et al. multiome snRNA-seq with STARsolo (Fig. 2b)

### `pipeline/02_archr_qc/`: ArchR project and QC

- `011_IsletEpigenetics_CreateArchRProject.R`: ArrowFiles and the initial ArchR project
- `012_IsletEpigenetics_CreateQCplots.R`: per-sample QC plots used to choose the cutoffs in `cutoffvalues.csv`
- `013_IsletEpigenetics_FilterQCandExclusions.R`: remove low-quality cells with the per-sample cutoffs
- `021_IsletEpigenetics_CreatePseudobulkpeaks.R`: LSI, Harmony, clustering, UMAP and the first MACS2 peak set
- `022_IsletEpigenetics_CreateFRIPQCplots.R`: per-sample FRIP plots used to choose the FRIP cutoffs
- `023_IsletEpigenetics_FilterQC.R`: remove cells below the per-sample FRIP cutoff
- `cutoffvalues.csv`: per-sample QC cutoffs

### `pipeline/03_doublets/`: doublet removal

- `031_IsletEpigenetics_CreateClustersforDoubletfinding.R`: LSI, Harmony and clusters for doublet detection
- `032_IsletEpigenetics_Doubletfinding-scDblFinder.R`: scDblFinder, 10 runs per sample
- `033_IsletEpigenetics_Doubletfinding-Amulet.R`: AMULET
- `034_IsletEpigenetics_Doubletfinding-ArchRDoubletfinding.R`: ArchR doublet scores
- `035_IsletEpigenetics_Doubletfinding_RemoveDoubletClusters.R`: remove doublet-enriched clusters
- `041_IsletEpigenetics_Second_Pseudobulk_and_peakcalling.R`: second pseudobulk and MACS2 peak calling
- `042_IsletEpigenetics_Doubletfinding_scDblFinder_2ndtime.R`: second scDblFinder run
- `043_IsletEpigenetics_endocrine_cluster_removal.R`, `043_IsletEpigenetics_non_endocrine_cluster_removal.R`: per-cluster doublet and QC overviews
- `044_IsletEpigenetics_remove_bad_clusters.R`: remove doublet and low-quality clusters, final project

### `pipeline/04_annotation/`: cell type annotation

- `051_isletepigenetics_labeling.R`: clustering, UMAP and cell type labels from gene score markers
- `052_isletepigenetics_PairedInsertionCounting.R`: paired-insertion counts (PICsnATAC) in the project peak set
- `053_isletepigenetics_export_to_anndata.R`: export to one AnnData file per sample
- `054_IsletEpigenetics_metadata.R`: harmonised donor metadata
- `benaglio_labels.R`, `chiou_labels.R`, `hpap_labels.R`, `wang_labels.R`: comparison with published cell type labels
- `maping_clusters_to_labels.R`: majority label per cluster
- `sex_validation2.R`: donor sex check from chrX accessibility and XIST
- `export_celltype_bigwigs.R`: BigWig tracks per cell type

### `pipeline/05_cycling/`: cycling cells and trajectories

- `cycling_labels.R`: cycling cell clusters in the final project
- `create_cycling_projects.R`: trajectory project of cycling and endocrine cells
- `cycling_clusters.R`: cycling cluster groups from gene score markers
- `module_scores.R`: endocrine marker and cell cycle module scores
- `run_slingshot.R`: Slingshot trajectories from the cycling cells to each endocrine cell type
- `wang_rna_mapping.R`: Wang et al. snRNA-seq mapped to the ATAC cells of the project

### `pipeline/06_donor_model/`: donor model

- `build_anndata.ipynb`: AnnData object and peak filter for the model
- `PeakANVI_donor_model.ipynb`: PeakANVI model, cell and donor embeddings
- `donor_pca_covariate_tests.R`: PCA of the donor embeddings and tests against donor covariates
- `donor_model_evaluation.R`: LISI and random forest evaluation of the donor model

### `pipeline/07_peak_associations/`: peak-donor PC associations

- `peak_pc_association.R`: per cell type DESeq2 likelihood ratio test of each donor PC against peak accessibility
- `significant_peaks.R`: significant peak-PC associations and their direction

### `pipeline/08_gwas/`: GWAS enrichment

- `gwas_reference_setup.sh`: install LDSC and TORUS, download reference files and UK Biobank summary statistics
- `ldsc_pc_peaks.R`: LDSC partitioned heritability of the PC-associated peaks
- `torus_pc_peaks.R`: TORUS enrichment of the PC-associated peaks

### `pipeline/09_abc/`: activity-by-contact model

- `131_IsletEpigenetics_ABC.R`: inputs for the ABC model

### `pipeline/10_chrombpnet/`: ChromBPNet

- `setup_chrombpnet_inputs.sh`: reference files and peak sets
- `alpha_beta_train_heldout_metadata.R`: training and held-out cells for the ND and T2D fine-tuning
- `train_celltype_models.ipynb`: one ChromBPNet model per cell type
- `prepare_finetune_peaks.ipynb`, `finetune_nd_t2d_models.ipynb`: ND and T2D fine-tuning of the alpha and beta models
- `motif_ablation_alpha_beta.ipynb`, `motif_ablation_all_celltypes.ipynb`: motif scan, enrichment and ablation (Fig. 5)
- `motif_activity_heatmap.ipynb`: motif activity across cell types
- `README.md`: run order, peak sets, versions and downloads of this step

### `figures/`

| Figure | Location |
|---|---|
| Fig. 1 | `figures/fig1/` (011 donor heatmap 1A; 012 donor violins 1B-F; 013 UMAP 1G; 014 gene score heatmap 1H; 015 browser tracks 1I) |
| ESM Fig. 1 | `figures/esm_fig1/` (histogram S1A; 016s and 017s label comparison S1F; 018s bar plots S1E); also 013 (S1B, S1C) and 015 (S1D) in `fig1/` |
| Fig. 2 | `figures/fig2/` (UMAP_genescore and cycling_umap_pseudotime 2A; cycling_genescore_* 2B (ATAC lines; `cycling_genescore_pseudotime_wang_rna.R` snRNA lines); DA 2C-E and S2B; barplot_cycling S2D-E) |
| ESM Fig. 2 | `figures/esm_fig2/` (026s S2A; tricycle_pseudotime S2B; 023 S2C) |
| Fig. 3 | `figures/fig3/` (umap_model 3B; n_correlated_peaks_pr_PC and donor_pca_covariate_heatmap 3C; latentspace_PC_plots 3D; gwas_heatmap 3F; `Figure3_clean.ipynb` enrichment and chromatin annotation) |
| ESM Fig. 3 | `figures/esm_fig3/` (lisi S3D); also umap_latent_space (S3C) and gwas_heatmap (S3I-J) in `fig3/` |
| Fig. 4 | `figures/fig4/Figure4_clean.ipynb` |
| Fig. 5 | `figures/fig5/Figure5_clean.ipynb`; PBX1 panels in `pipeline/10_chrombpnet/motif_ablation_all_celltypes.ipynb` |

`figures/colours.R`: shared colour palettes.

# Data

All datasets in the atlas are public, openly or through managed access:

| Study | Source | Accession |
|---|---|---|
| Augsornworawat | GEO | GSE199636 |
| Benaglio | GEO | GSE205853 |
| Chiou | GEO | GSE160472, GSE163160 |
| HPAP | PANC-DB (Human Pancreas Analysis Program) | https://hpap.pmacs.upenn.edu/ |
| Rai (FUSION) | dbGaP, controlled access | phs001188.v2.p1 |
| Thibodeau | GEO | GSE165212 |
| Wang | GEO (snATAC-seq and multiome); ENA (multiome snRNA-seq FASTQ) | GSE169453, GSE200044; PRJNA822928 |
| Weng | GEO | GSE195523 |

Only samples from untreated, purified human islets of adult donors were used; see the Methods of the paper.

# Software

- Alignment: STAR 2.7.10a (STARsolo), samtools 1.17, sinto, htslib 1.21, sra-tools, aria2; genome GRCh38 with GENCODE v35 for snRNA-seq.
- R analysis: ArchR (hg38) with MACS2 2.2.9.1, Harmony, Seurat 4.0.0 (Matrix 1.6-1.1), SeuratDisk, scDblFinder, PICsnATAC, slingshot, tricycle, UCell, mgcv, DESeq2, mapgen, liftOver. TODO: versions of the other R packages (`sessionInfo()` of the analysis environment).
- PeakANVI donor model and AnnData notebooks: Python with PyTorch, scanpy and anndata.
- GWAS: LDSC (commit aa33296, `env/ldsc.yml`, Python 2.7) with baselineLD v2.3, TORUS; installed by `pipeline/08_gwas/gwas_reference_setup.sh`.
- ChromBPNet: chrombpnet 1.0.1 (Python 3.8, TensorFlow 2.8.0) for training; bpnet-lite, tangermeme and memelite for interpretation. Versions in `pipeline/10_chrombpnet/README.md`.
- ABC: ABC-Enhancer-Gene-Prediction Snakemake workflow (https://github.com/GreenleafLab/ABC-Enhancer-Gene-Prediction-CustomRegions).

# License

MIT, see `LICENSE`.

# Contact

Jesper Grud Skat Madsen (corresponding author), Department of Biochemistry and Molecular Biology, University of Southern Denmark: jgsm@bmb.sdu.dk
