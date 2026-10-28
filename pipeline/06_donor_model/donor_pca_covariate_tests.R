# ==============================================================================
# donor_pca_covariate_tests.R
# Purpose : PCA of the donor embeddings from the PeakANVI model and tests of each donor
#           PC against donor covariates (linear model, Wilcoxon, Kruskal-Wallis and
#           pairwise Wilcoxon tests).
# Inputs  : data/public/scATACseq/JM/samples.csv (donor embedding, exported from the
#           trained PeakANVI model, PeakANVI_donor_model.ipynb);
#           data/public/scATACseq/Project/meta.csv (donor metadata, same as ESM Table 1)
# Outputs : data/public/scATACseq/correlations/pca.csv (donor PC scores, used by
#           donor_model_evaluation.R, 07_peak_associations and the Fig. 3 scripts);
#           data/public/scATACseq/correlations/PC_covariate_correlation.rds
# ==============================================================================

data_dir      <- normalizePath("data/public", mustWork = FALSE)     # public/processed data

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(purrr)
  library(vroom)
})

# -----------------------------
# Load + build donor matrix
# -----------------------------
samples <- read.csv(file.path(data_dir, "scATACseq/JM/samples.csv"), stringsAsFactors = FALSE)

donor_data <- as.matrix(apply(samples[, 1:20, drop = FALSE], 2, as.numeric))
rownames(donor_data) <- as.character(samples[, 21])

donor_meta0 <- as.data.frame(samples[, 21:28, drop = FALSE], stringsAsFactors = FALSE)
colnames(donor_meta0) <- c("ID", "study", "disease", "gender", "age", "bmi", "hba1c", "ethnicity")
donor_meta0$ID    <- as.character(donor_meta0$ID)
donor_meta0$study <- as.character(donor_meta0$study)
donor_meta0$study[donor_meta0$study == "FUSION"] <- "Rai"

meta <- vroom::vroom(file.path(data_dir, "scATACseq/Project/meta.csv"), show_col_types = FALSE)

donor_meta <- donor_meta0 %>%
  dplyr::select(ID, study) %>%
  dplyr::left_join(meta, by = "ID") %>%
  dplyr::mutate(
    BMI   = as.numeric(dplyr::coalesce(BMI, donor_meta0$bmi[match(ID, donor_meta0$ID)])),
    Age   = as.numeric(dplyr::coalesce(Age, donor_meta0$age[match(ID, donor_meta0$ID)])),
    hba1c = as.numeric(dplyr::coalesce(hba1c, donor_meta0$hba1c[match(ID, donor_meta0$ID)])),
    Gender = as.factor(dplyr::coalesce(Gender, donor_meta0$gender[match(ID, donor_meta0$ID)])),
    ethnicity_harmonized = as.factor(dplyr::coalesce(ethnicity_harmonized, donor_meta0$ethnicity[match(ID, donor_meta0$ID)])),

    study = as.factor(study),
    harmonized_disease = as.factor(harmonized_disease),
    harmonized_disease_hba1c = as.factor(harmonized_disease_hba1c),
    metformin_harmonized = as.factor(metformin_harmonized),

    fraction_endocrine = as.numeric(fraction_endocrine),
    beta_fraction = as.numeric(beta_fraction),
    beta_fraction_endocrine = as.numeric(beta_fraction_endocrine)
  )

# Align
keep <- rownames(donor_data) %in% donor_meta$ID
donor_data <- donor_data[keep, , drop = FALSE]
donor_meta <- donor_meta[match(rownames(donor_data), donor_meta$ID), , drop = FALSE]
stopifnot(identical(rownames(donor_data), donor_meta$ID))

# PCA
donor_pca <- prcomp(donor_data, scale. = TRUE)
dir.create(file.path(data_dir, "scATACseq/correlations"), recursive = TRUE, showWarnings = FALSE)
write.csv(donor_pca$x, file.path(data_dir, "scATACseq/correlations/pca.csv"))
pcs <- 1:20
X <- donor_pca$x[, pcs, drop = FALSE]  # all donors x PCs

# Stats utilities
prep_factor <- function(g, level_order = NULL, drop_levels = NULL) {
  g <- as.character(g)
  g[g == ""] <- NA_character_
  if (!is.null(drop_levels)) g[g %in% drop_levels] <- NA_character_

  if (!is.null(level_order)) {
    g <- factor(g, levels = unique(level_order))
  } else {
    g <- factor(g)
  }
  droplevels(g)
}

lm_stats <- function(y, x) {
  ok <- is.finite(y) & is.finite(x)
  y <- y[ok]; x <- x[ok]
  n <- length(y)
  if (n < 8) return(NULL)

  fit <- lm(y ~ x)
  sm  <- summary(fit)

  tibble(
    test = "lm",
    comparison = NA_character_,
    n = n,
    estimate = unname(coef(fit)[["x"]]),
    se = unname(sm$coefficients["x", "Std. Error"]),
    statistic = unname(sm$coefficients["x", "t value"]),
    df = unname(fit$df.residual),
    r2 = unname(sm$r.squared),
    median_diff = NA_real_,
    cliffs_delta = NA_real_,
    epsilon2 = NA_real_,
    p_value = unname(sm$coefficients["x", "Pr(>|t|)"])
  )
}

wilcox_stats <- function(y, g, level_order = NULL, drop_levels = NULL) {
  ok <- is.finite(y) & !is.na(g)
  y <- y[ok]
  g <- prep_factor(g[ok], level_order = level_order, drop_levels = drop_levels)

  # IMPORTANT: drop NA introduced by level_order filtering
  keep <- !is.na(g)
  y <- y[keep]
  g <- droplevels(g[keep])

  if (nlevels(g) != 2) return(NULL)
  if (!is.null(level_order)) g <- factor(g, levels = unique(level_order))

  n1 <- sum(g == levels(g)[1], na.rm = TRUE)
  n2 <- sum(g == levels(g)[2], na.rm = TRUE)
  if (n1 < 2 || n2 < 2) return(NULL)

  wt <- suppressWarnings(wilcox.test(y ~ g))  # old default behavior
  W  <- unname(wt$statistic)

  U1 <- W - n1 * (n1 + 1) / 2
  delta <- 1 - 2 * (U1 / (n1 * n2))

  med1 <- stats::median(y[g == levels(g)[1]], na.rm = TRUE)
  med2 <- stats::median(y[g == levels(g)[2]], na.rm = TRUE)

  tibble::tibble(
    test = "wilcox",
    comparison = paste0(levels(g)[2], " vs ", levels(g)[1]),
    n = length(y),
    estimate = NA_real_,
    se = NA_real_,
    statistic = unname(wt$statistic),
    df = NA_real_,
    r2 = NA_real_,
    median_diff = med2 - med1,
    cliffs_delta = delta,
    epsilon2 = NA_real_,
    p_value = unname(wt$p.value)
  )
}

kruskal_stats <- function(y, g, level_order = NULL, drop_levels = NULL) {
  ok <- is.finite(y) & !is.na(g)
  y <- y[ok]
  g <- prep_factor(g[ok], level_order = level_order, drop_levels = drop_levels)

  # IMPORTANT: drop NA introduced by level_order filtering
  keep <- !is.na(g)
  y <- y[keep]
  g <- droplevels(g[keep])

  if (nlevels(g) < 2 || length(g) < 3) return(NULL)

  kt <- suppressWarnings(kruskal.test(y ~ g))
  H  <- unname(kt$statistic)
  k  <- nlevels(g)
  n  <- length(g)
  eps2 <- (H - k + 1) / (n - k)

  tibble::tibble(
    test = "kruskal",
    comparison = "overall",
    n = n,
    estimate = NA_real_,
    se = NA_real_,
    statistic = H,
    df = unname(kt$parameter),
    r2 = NA_real_,
    median_diff = NA_real_,
    cliffs_delta = NA_real_,
    epsilon2 = eps2,
    p_value = unname(kt$p.value)
  )
}

pairwise_wilcox_3lvl <- function(y, g, level_order, drop_levels = NULL) {
  ok <- is.finite(y) & !is.na(g)
  y <- y[ok]
  g <- prep_factor(g[ok], level_order = level_order, drop_levels = drop_levels)

  # drop NA introduced by level_order filtering
  keep <- !is.na(g)
  y <- y[keep]
  g <- droplevels(g[keep])

  lv <- levels(g)
  if (length(lv) < 2) return(NULL)

  pairs <- combn(lv, 2, simplify = FALSE)
  out <- purrr::map_dfr(pairs, function(p) {
    kk <- g %in% p
    y2 <- y[kk]
    g2 <- factor(g[kk], levels = p)

    s <- wilcox_stats(y2, g2, level_order = p)
    if (is.null(s)) return(NULL)
    s$comparison <- paste0(p[1], " vs ", p[2])
    s
  })

  if (nrow(out) == 0) return(NULL)
  out %>% dplyr::mutate(p_adj_bh = p.adjust(p_value, method = "BH"))
}

# Test specification
spec_overall <- tibble(
  condition = c(
    "Age","BMI","hba1c","fraction_endocrine","beta_fraction_endocrine",
    "Gender","disease_ND_vs_T2D","metformin","ethnicity","study",
    "disease_hba1c_overall"
  ),
  var = c(
    "Age","BMI","hba1c","fraction_endocrine","beta_fraction_endocrine",
    "Gender","harmonized_disease","metformin_harmonized","ethnicity_harmonized","study",
    "harmonized_disease_hba1c"
  ),
  method = c(
    "lm","lm","lm","lm","lm",
    "wilcox","wilcox","wilcox","kruskal","kruskal",
    "kruskal"
  ),
  level_order = list(
    NULL, NULL, NULL, NULL, NULL,
    c("F","M"),
    c("ND","T2D"),                    # drops T1D implicitly
    NULL,
    NULL,
    NULL,
    c("ND","PRE-T2D","T2D")
  ),
  drop_levels = list(
    NULL, NULL, NULL, NULL, NULL,
    NULL,
    NULL,
    c("Unknown","unknown"),
    c("Unknown","unknown"),
    NULL,
    NULL
  )
)

spec_pairwise <- tibble(
  condition   = "disease_hba1c_pairwise",
  var         = "harmonized_disease_hba1c",
  level_order = list(c("ND","PRE-T2D","T2D")),
  drop_levels = list(NULL)
)

dispatch <- list(
  lm      = function(y, meta, var, level_order, drop_levels) lm_stats(y, meta[[var]]),
  wilcox  = function(y, meta, var, level_order, drop_levels) wilcox_stats(y, meta[[var]], level_order, drop_levels),
  kruskal = function(y, meta, var, level_order, drop_levels) kruskal_stats(y, meta[[var]], level_order, drop_levels)
)

# Run tests per PC
overall_results <- purrr::map_dfr(pcs, function(i) {
  y <- X[, i]
  purrr::pmap_dfr(
    spec_overall,
    function(condition, var, method, level_order, drop_levels) {
      fn <- dispatch[[method]]
      s <- fn(y, donor_meta, var, level_order, drop_levels)
      if (is.null(s)) return(NULL)

      s %>%
        mutate(pc = paste0("PC", i), condition = condition) %>%
        dplyr::select(pc, condition, test, comparison, n, df, statistic, p_value,
                      estimate, se, r2, median_diff, cliffs_delta, epsilon2)
    }
  )
})

pairwise_results <- purrr::map_dfr(pcs, function(i) {
  y <- X[, i]
  purrr::pmap_dfr(
    spec_pairwise,
    function(condition, var, level_order, drop_levels) {
      pw <- pairwise_wilcox_3lvl(y, donor_meta[[var]], level_order, drop_levels)
      if (is.null(pw)) return(NULL)

      pw %>%
        mutate(pc = paste0("PC", i), condition = condition) %>%
        dplyr::select(pc, condition, test, comparison, n, statistic, p_value, p_adj_bh,
                      median_diff, cliffs_delta)
    }
  )
})

res <- list(overall = overall_results, pairwise = pairwise_results)

saveRDS(res, file.path(data_dir, "scATACseq/correlations/PC_covariate_correlation.rds"))
