# ==============================================================================
# donor_model_evaluation.R
# Purpose : Evaluation of the donor model (T1D donors excluded): sex separation on
#           donor PC12 and PC15 (ESM Fig. 3e), LISI of donor, study and disease on the
#           donor PCs (printed), and random forest prediction of study, sex and
#           disease from the donor PCs, the LSI, Harmony and the PeakANVI cell latent
#           space (ESM Fig. 3d).
# Inputs  : data/public/scATACseq/ArchR/IsletProject;
#           data/public/scATACseq/correlations/pca.csv (donor_pca_covariate_tests.R);
#           data/public/scATACseq/JM/cells.csv (cell latent space, exported from the
#           trained PeakANVI model, PeakANVI_donor_model.ipynb)
# Outputs : data/public/scATACseq/JM/S3E.pdf, S3D.pdf
# ==============================================================================

data_dir      <- normalizePath("data/public", mustWork = FALSE)     # public/processed data

library(ArchR)
library(lisi)
library(dplyr)
library(ggplot2)
library(ranger)
library(pROC)

proj <- loadArchRProject(file.path(data_dir, "scATACseq/ArchR/IsletProject/"))

## Sex separation on donor PC12 and PC15
pca <- read.csv(file.path(data_dir, "scATACseq/correlations/pca.csv"), row.names = 1)
pca <- pca[rownames(pca) %in% unique(proj$Sample), ]
pca <- pca[order(rownames(pca)), ]
ccd <- as.data.frame(proj@cellColData)
pca <- pca[ !(rownames(pca) %in% ccd[ ccd$disease_harmonized_hba1c == "T1D", "Sample"]),]
ccd <- ccd[ match(rownames(pca), ccd$Sample),]

pdf(file.path(data_dir, "scATACseq/JM/S3E.pdf"), width = 10, height = 10, useDingbats = FALSE)
plot(x = pca$PC12, y = pca$PC15, col = ifelse(ccd$gender == "M", "darkgreen", "darkred"), pch = 16)
dev.off()

## LISI on the donor PCs, scaled to 0 (no mixing) - 1 (full mixing)
lisi <- compute_lisi(X = pca, meta_data = ccd, label_colnames = c("donor_harmonized", "study", "disease_harmonized"), perplexity = 3)

median((lisi[,"study"] - 1) / (length(unique(ccd$study)) - 1))
median((lisi[,"disease_harmonized"] - 1) / (length(unique(ccd$disease_harmonized)) - 1))
median((lisi[,"donor_harmonized"] - 1) / (length(unique(ccd$donor_harmonized)) - 1))

## Random forest prediction of study, sex and disease from each embedding space
# Stratified k-fold split
make_stratified_folds <- function(target, k = 5, seed = 1) {
  set.seed(seed)
  target <- as.factor(target)
  n <- length(target)
  fold_id <- integer(n)
  for (cls in levels(target)) {
    idx <- which(target == cls)
    idx <- sample(idx, length(idx))
    fold_id[idx] <- rep(1:k, length.out = length(idx))
  }

  fold_id
}

pca <- read.csv(file.path(data_dir, "scATACseq/correlations/pca.csv"), row.names = 1)
pca <- pca[rownames(pca) %in% unique(proj$Sample), ]
pca <- pca[order(rownames(pca)), ]

cells <- read.csv(file.path(data_dir, "scATACseq/JM/cells.csv"), check.names = FALSE, stringsAsFactors = FALSE)
ccd <- as.data.frame(proj@cellColData)
rownames(cells) <- rownames(ccd)

# First 20 dimensions of each space (LSI without component 1)
LSI = proj@reducedDims[[1]]$matSVD[,2:21]
Harm = proj@reducedDims[[2]]$matDR[,1:20]
Zcells = cells[,1:20]

pca <- pca[ !(rownames(pca) %in% ccd[ ccd$disease_harmonized_hba1c == "T1D", "Sample"]),]
LSI <- LSI[ !(rownames(LSI) %in% rownames(ccd[ ccd$disease_harmonized_hba1c == "T1D",])),]
Harm <- Harm[ !(rownames(Harm) %in% rownames(ccd[ ccd$disease_harmonized_hba1c == "T1D",])),]
Zcells <- Zcells[ !(rownames(Zcells) %in% rownames(ccd[ ccd$disease_harmonized_hba1c == "T1D",])),]
ccd <- ccd[ ccd$disease_harmonized_hba1c != "T1D",]

inputs = list(
  "Zdonor" = list("data" = pca, "metadata" = ccd[ match(rownames(pca), ccd$Sample),]),
  "LSI" = list("data" =  LSI, "metadata" = ccd[ match(rownames(LSI), rownames(ccd)),]),
  "Harmony" = list("data" = Harm, "metadata" = ccd[ match(rownames(Harm), rownames(ccd)),]),
  "Zcells" = list("data" = Zcells, "metadata" = ccd)
)

targets <- c("study", "gender", "disease_harmonized_hba1c")
k <- 5
first = TRUE

# ranger draws its seed from R's random number generator
set.seed(1)
for (inputidx in 1:4) {
  for (target in targets) {
    dat <- as.data.frame(inputs[[inputidx]]$data)
    md <- inputs[[inputidx]]$metadata
    dat$target <- md[,target]
    keep <- names(which(table(dat$target) >= 2))
    dat <- dat[dat$target %in% keep, ]
    dat$target <- factor(dat$target)
    folds <- make_stratified_folds(dat$target, k = k, seed = 123)
    truth <- c()

    for(i in 1:k) {
      train_data <- dat[folds != i, ]
      test_data  <- dat[folds == i, ]

      # Scale with the training fold mean and SD
      means <- apply(train_data[,1:20], 2, FUN = "mean")
      sds <- apply(train_data[,1:20], 2, FUN = "sd")
      train_data[,1:20] = t(((t(train_data[,1:20]) - means) / sds))
      test_data[,1:20] = t(((t(test_data[,1:20]) - means) / sds))

      fit = ranger(y = train_data$target, x = train_data[,1:20], num.threads = 15, probability = TRUE)
      pred <- predict(fit, num.threads = 15, data = test_data[,1:20])$predictions
      if (i == 1) {
        pred_mat = pred
      } else {
        pred_mat <- rbind(pred_mat, pred)
      }
      truth <- c(truth, as.vector(test_data$target))
    }

    # Multiclass AUC over the held-out predictions
    df_tmp = data.frame(AUC = multiclass.roc(response = truth, predictor = pred_mat)$auc[1], Variable = target, Space = names(inputs)[inputidx])
    if (first) {
      res_df = df_tmp
      first = FALSE
    } else {
      res_df = rbind(res_df, df_tmp)
    }
  }
}
RANGER <- res_df
# Score: 1 = no better than chance (AUC 0.5), 0 = perfect prediction (AUC 1)
RANGER$score <- 1 - (RANGER$AUC - 0.5) / 0.5

pdf(file.path(data_dir, "scATACseq/JM/S3D.pdf"), width = 10, height = 10, useDingbats = FALSE)
print(ggplot(data = RANGER, aes(x = Space, y = score)) + geom_bar(stat = "identity") + facet_wrap(~Variable))
dev.off()
