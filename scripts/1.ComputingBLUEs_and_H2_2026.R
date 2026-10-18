# ==============================================================================
# AUGMENTED RCBD GENOMIC SELECTION PIPELINE (2026 PROCESSING TOMATO TRIAL)
# Directory: C:/Users/techz/OneDrive/Documents/PhD-GS25/GP5/
# Pheno File: data_final.csv (275 plots: 95 Check plots + 180 Test plots)
# Geno File:  PhDGS25Geno176ReadyForGS.csv (176 genotypes x 1905 SNPs)
# ==============================================================================

# 1. Setup Environment & Load Packages -----------------------------------------
work_dir <- "C:/Users/techz/OneDrive/Documents/PhD-GS25/GP5"
setwd(work_dir)

required_pkgs <- c("tidyverse", "lme4", "rrBLUP", "data.table", "ggplot2")
new_pkgs <- required_pkgs[!(required_pkgs %in% installed.packages()[, "Package"])]
if (length(new_pkgs)) install.packages(new_pkgs)

library(tidyverse)
library(lme4)
library(rrBLUP)
library(data.table)
library(ggplot2)

# 2. Load Data & Define Target Traits ------------------------------------------
pheno_raw <- read.csv("data_final.csv", stringsAsFactors = FALSE)
geno_raw  <- fread("PhDGS25Geno176ReadyForGS.csv", data.table = FALSE)

traits <- c("TotCH", "Fv_Fm", "Sm", "Pi_Abs", "plant_weight_g", "fruit_weight_g")

# Clean strings and prepare Augmented RCBD grouping factors
# Note: Year is ignored (all 2026); Column (A-K) and Block (1-11) are identical
pheno_df <- pheno_raw %>%
  select(ID, Genotype, Type, Column, Row, Block, all_of(traits)) %>%
  mutate(
    ID       = trimws(as.character(ID)),
    Genotype = trimws(as.character(Genotype)),
    Type     = factor(trimws(Type), levels = c("Check", "Test")),
    Block    = as.factor(Block),
    Row      = as.factor(Row),
    # Helper factors for Augmented RCBD mixed model outlier screening & H2:
    # Check_ID separates each of the 7 replicated Checks from the Test population
    Check_ID = ifelse(Type == "Check", Genotype, "Test_Population"),
    # Test_ID assigns random genetic effects to unreplicated Test lines
    Test_ID  = ifelse(Type == "Test", Genotype, paste0("Check_", Genotype))
  )

# ==============================================================================
# STAGE 1: BIOLOGICAL FILTERING & AUGMENTED RCBD OUTLIER CLEANING
# ==============================================================================

# 1A. Biological hard bounds (removes 0 g fruit weight in 26-4094 ELBA & misfires)
pheno_clean <- pheno_df %>%
  mutate(
    Fv_Fm          = ifelse(Fv_Fm < 0.40 | Fv_Fm > 0.88, NA, Fv_Fm),
    TotCH          = ifelse(TotCH <= 0, NA, TotCH),
    Sm             = ifelse(Sm <= 0, NA, Sm),
    Pi_Abs         = ifelse(Pi_Abs <= 0, NA, Pi_Abs),
    plant_weight_g = ifelse(plant_weight_g <= 0, NA, plant_weight_g),
    fruit_weight_g = ifelse(fruit_weight_g <= 0, NA, fruit_weight_g)
  )

# 1B. Model-based Studentized Residual Outlier Removal (|s_res| > 3.0 SD)
outlier_log <- data.frame()

for (t in traits) {
  form_out <- as.formula(paste(t, "~ Check_ID + (1|Block) + (1|Test_ID)"))
  fit_out  <- lmer(form_out, data = pheno_clean, na.action = na.exclude,
                   control = lmerControl(optimizer = "bobyqa", calc.derivs = FALSE))
  
  # Standardized Pearson residuals
  s_res <- resid(fit_out, type = "pearson") / sigma(fit_out)
  out_idx <- which(abs(s_res) > 3.0)
  
  if (length(out_idx) > 0) {
    for (idx in out_idx) {
      outlier_log <- rbind(outlier_log, data.frame(
        Trait    = t,
        ID       = pheno_clean$ID[idx],
        Genotype = pheno_clean$Genotype[idx],
        Type     = as.character(pheno_clean$Type[idx]),
        Block    = as.character(pheno_clean$Block[idx]),
        Raw_Val  = pheno_clean[idx, t],
        Std_Res  = round(s_res[idx], 2)
      ))
    }
    message(sprintf("[%s] Removed %d outlier(s): %s", t, length(out_idx),
                    paste(pheno_clean$ID[out_idx], collapse = ", ")))
    pheno_clean[out_idx, t] <- NA
  } else {
    message(sprintf("[%s] 0 outliers detected (all |residuals| <= 3.0 SD).", t))
  }
}

if (nrow(outlier_log) > 0) {
  write.csv(outlier_log, "Stage1_Removed_Outliers_Log.csv", row.names = FALSE)
}

# ==============================================================================
# STAGE 2: EXTRACT BLOCK-ADJUSTED BLUEs & HERITABILITY (H2)
# ==============================================================================

blues_list <- list()
stage2_summary <- data.frame(
  Trait      = traits,
  N_Genos    = NA,
  Mean_BLUE  = NA,
  SD_BLUE    = NA,
  CV_BLUE_pct= NA,
  Sigma2_G   = NA,
  Sigma2_E   = NA,
  H2_Aug     = NA,
  H2_Cullis  = NA
)

for (i in seq_along(traits)) {
  t <- traits[i]
  sub_df <- pheno_clean %>% filter(!is.na(.data[[t]]))
  
  # --- Model A: Genotype FIXED, Block RANDOM -> Extract Exact BLUEs ---
  form_blue <- as.formula(paste(t, "~ 0 + Genotype + (1|Block)"))
  fit_blue  <- lmer(form_blue, data = sub_df,
                    control = lmerControl(optimizer = "bobyqa", calc.derivs = FALSE))
  
  fe <- fixef(fit_blue)
  names(fe) <- sub("^Genotype", "", names(fe))
  blues_list[[t]] <- data.frame(Genotype = names(fe), val = as.numeric(fe))
  colnames(blues_list[[t]])[2] <- t
  
  # --- Model B: Augmented Mixed Model -> Genetic Variance & Heritability ---
  # Check_ID fixed + Genotype (Test genetic variance) random + Block random
  form_vc <- as.formula(paste(t, "~ Check_ID + (1|Block) + (1|Genotype)"))
  fit_vc  <- lmer(form_vc, data = sub_df,
                  control = lmerControl(optimizer = "bobyqa", calc.derivs = FALSE))
  
  vc    <- as.data.frame(VarCorr(fit_vc))
  var_g <- vc$vcov[vc$grp == "Genotype"]
  var_e <- vc$vcov[vc$grp == "Residual"]
  
  # Standard Augmented Broad-Sense Heritability
  h2_aug <- ifelse((var_g + var_e) > 0, var_g / (var_g + var_e), 0)
  
  # Cullis Heritability: 1 - (mean PEV / Sigma2_G)
  re_geno <- ranef(fit_vc, condVar = TRUE)$Genotype
  pev     <- attr(re_geno, "postVar")[1, 1, ]
  h2_cullis <- ifelse(var_g > 0, max(0, 1 - (mean(pev) / var_g)), 0)
  
  stage2_summary$N_Genos[i]     <- length(fe)
  stage2_summary$Mean_BLUE[i]   <- round(mean(fe, na.rm = TRUE), 2)
  stage2_summary$SD_BLUE[i]     <- round(sd(fe, na.rm = TRUE), 2)
  stage2_summary$CV_BLUE_pct[i] <- round((sd(fe, na.rm = TRUE) / abs(mean(fe, na.rm = TRUE))) * 100, 2)
  stage2_summary$Sigma2_G[i]    <- round(var_g, 4)
  stage2_summary$Sigma2_E[i]    <- round(var_e, 4)
  stage2_summary$H2_Aug[i]      <- round(h2_aug, 3)
  stage2_summary$H2_Cullis[i]   <- round(h2_cullis, 3)
}

# Merge BLUEs across all 6 traits
BLUE_df <- reduce(blues_list, full_join, by = "Genotype")
geno_meta <- pheno_df %>% distinct(Genotype, Type)
BLUE_df <- right_join(geno_meta, BLUE_df, by = "Genotype")

write.csv(BLUE_df, "Stage2_Tomato_AugmentedRCBD_BLUEs.csv", row.names = FALSE)

# ==============================================================================
# STAGE 3: MAP 25-XXXX GENO IDs TO 26-XXXX PHENO IDs & BUILD VANRADEN G-MATRIX
# ==============================================================================

# Map '25-4000' -> '26-4000' -> Genotype Name ('VOLARE', etc.)
geno_id_26 <- sub("^25-", "26-", trimws(as.character(geno_raw$Genotype)))
id_lookup  <- setNames(pheno_df$Genotype, pheno_df$ID)
mapped_geno_names <- unname(id_lookup[geno_id_26])

# Verify mapping completeness
if (any(is.na(mapped_geno_names))) {
  warning("Some 25-XXXX IDs did not map to 26-XXXX IDs in data_final.csv!")
}

# Build numeric marker matrix M_mat coded as {-1, 0, 1} for rrBLUP
# PhDGS25Geno176ReadyForGS.csv is coded as {0, 1, 2}, so subtract 1
M_mat <- as.matrix(geno_raw[, -1]) - 1
rownames(M_mat) <- mapped_geno_names

# Align BLUEs and SNP Matrix (retains the 172 genotyped lines that have BLUEs)
common_genos <- intersect(BLUE_df$Genotype, rownames(M_mat))
cat(sprintf("\nSuccessfully aligned %d genotyped lines with field BLUEs (7 Checks + %d Tests).\n",
            length(common_genos), sum(BLUE_df$Type[match(common_genos, BLUE_df$Genotype)] == "Test")))

BLUE_aligned <- BLUE_df[match(common_genos, BLUE_df$Genotype), ]
rownames(BLUE_aligned) <- common_genos
M_aligned    <- M_mat[common_genos, ]

# Z-score standardize BLUEs (Mean = 0, SD = 1) for multi-scale numerical stability
Y_scaled <- scale(as.matrix(BLUE_aligned[, traits]))
rownames(Y_scaled) <- common_genos

# Compute VanRaden Additive Relationship Matrix (G) with MAF >= 0.05
G_mat <- A.mat(M_aligned, min.MAF = 0.05, max.missing = 0.20, impute.method = "mean")
G_mat <- G_mat + diag(1e-4, nrow(G_mat)) # Small ridge for numerical stability

# ==============================================================================
# STAGE 4: GBLUP CROSS-VALIDATION (5-FOLD 20 REPS + LEAVE-ONE-OUT CV)
# ==============================================================================

set.seed(2026)
n_folds <- 5
n_reps  <- 20

cv_5fold_all  <- matrix(NA, nrow = n_reps, ncol = length(traits), dimnames = list(NULL, traits))
cv_5fold_test <- matrix(NA, nrow = n_reps, ncol = length(traits), dimnames = list(NULL, traits))
loocv_all     <- setNames(rep(NA, length(traits)), traits)
loocv_test    <- setNames(rep(NA, length(traits)), traits)

is_test <- BLUE_aligned$Type == "Test"

for (t in traits) {
  y <- Y_scaled[, t]
  valid_idx <- which(!is.na(y)) # 172 for physiology traits, 171 for weights
  
  # --- 4A. Repeated 5-Fold Cross-Validation (20 Reps) ---
  for (rep in 1:n_reps) {
    folds <- rep(NA, length(y))
    folds[valid_idx] <- sample(rep(1:n_folds, length.out = length(valid_idx)))
    gebv_5f <- rep(NA, length(y))
    
    for (f in 1:n_folds) {
      val_idx <- which(folds == f)
      y_train <- y
      y_train[val_idx] <- NA
      
      fit_gs <- mixed.solve(y = y_train, K = G_mat)
      gebv_5f[val_idx] <- as.numeric(fit_gs$u[val_idx])
    }
    cv_5fold_all[rep, t]  <- cor(gebv_5f, y, use = "complete.obs")
    cv_5fold_test[rep, t] <- cor(gebv_5f[is_test], y[is_test], use = "complete.obs")
  }
  
  # --- 4B. Leave-One-Out Cross-Validation (LOOCV) ---
  gebv_loo <- rep(NA, length(y))
  for (idx in valid_idx) {
    y_train <- y
    y_train[idx] <- NA
    fit_loo <- mixed.solve(y = y_train, K = G_mat)
    gebv_loo[idx] <- as.numeric(fit_loo$u[idx])
  }
  loocv_all[t]  <- cor(gebv_loo, y, use = "complete.obs")
  loocv_test[t] <- cor(gebv_loo[is_test], y[is_test], use = "complete.obs")
}

# ==============================================================================
# STAGE 5: EXPORT CONFERENCE SUMMARY TABLE, GEBVs, & VISUALIZATIONS
# ==============================================================================

final_gs_table <- stage2_summary %>%
  mutate(
    r_5Fold_All      = round(colMeans(cv_5fold_all, na.rm = TRUE), 3),
    SD_5Fold_All     = round(apply(cv_5fold_all, 2, sd, na.rm = TRUE), 3),
    r_5Fold_TestOnly = round(colMeans(cv_5fold_test, na.rm = TRUE), 3),
    SD_5Fold_TestOnly= round(apply(cv_5fold_test, 2, sd, na.rm = TRUE), 3),
    r_LOOCV_All      = round(loocv_all[Trait], 3),
    r_LOOCV_TestOnly = round(loocv_test[Trait], 3)
  )

print(final_gs_table)
write.csv(final_gs_table, "Conference_Paper_GS_Summary_2026.csv", row.names = FALSE)

# Fit full-dataset GBLUP models to export final GEBVs for all 176 genotyped lines
G_full <- A.mat(M_mat, min.MAF = 0.05) + diag(1e-4, nrow(M_mat))
gebv_full_df <- data.frame(
  Genotype_25ID = geno_raw$Genotype,
  Plot_26ID     = geno_id_26,
  Genotype      = rownames(M_mat),
  Type          = pheno_df$Type[match(rownames(M_mat), pheno_df$Genotype)]
)

for (t in traits) {
  y_full <- BLUE_df[[t]][match(rownames(M_mat), BLUE_df$Genotype)]
  fit_full <- mixed.solve(y = y_full, K = G_full)
  gebv_full_df[[paste0("GEBV_", t)]] <- round(as.numeric(fit_full$u) + as.numeric(fit_full$beta), 3)
}
write.csv(gebv_full_df, "Stage5_All176_Genomic_Breeding_Values_GEBVs.csv", row.names = FALSE)

# Plot 5-Fold Cross-Validation Predictive Ability Boxplot
cv_long <- as.data.frame(cv_5fold_all) %>%
  pivot_longer(cols = everything(), names_to = "Trait", values_to = "r_MG") %>%
  mutate(Trait = factor(Trait, levels = traits))

p_box <- ggplot(cv_long, aes(x = Trait, y = r_MG, fill = Trait)) +
  geom_boxplot(width = 0.55, alpha = 0.85, outlier.shape = 16) +
  geom_jitter(width = 0.12, alpha = 0.4, size = 1.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  theme_minimal(base_size = 13) +
  labs(
    title = "Genomic Prediction Ability (GBLUP) Across 6 Tomato Traits",
    subtitle = "172 Genotyped Lines x 1,901 SNPs (MAF >= 0.05) | 5-Fold CV (20 Reps)",
    x = "Trait",
    y = "Predictive Ability (r_MG)"
  ) +
  theme(legend.position = "none", panel.grid.minor = element_blank())

ggsave("Tomato_GS_Predictive_Ability_Boxplot.png", plot = p_box, width = 8.5, height = 5.2, dpi = 300)