library(here)
library(tidyverse)
library(lme4)
library(rrBLUP)
library(data.table)
library(ggplot2)
library(augmentedRCBD)
library(sommer)
library(heritable)
########################################### script arguments ########
pheno_f<-here("data/pheno/data_final.csv") ## data 2026

######################################## main script ################
pheno_raw <- read.csv(pheno_f, stringsAsFactors = FALSE)
traits <- c("TotCH", "Fv_Fm", "Sm", "Pi_Abs", "plant_weight_g", "fruit_weight_g")

# Clean strings and prepare Augmented RCBD grouping factors
# Note: Year is ignored (all 2026); Column (A-K) and Block (1-11) are identical
pheno_df <- pheno_raw %>%
  dplyr::select(ID, Genotype, Type, Column, Row, Block, all_of(traits)) %>%
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
#pheno_clean[,7:12][pheno_clean$Genotype=="VOLARE",]<-NA
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
  H2_Cullis  = NA,
  H2_AugmentedRCBD     = NA,
  H2_sommer =NA
)
heri<-data.frame(
  Trait=traits,
  Cullis=NA,
  Oakey=NA,
  Delta=NA,
  Piepho=NA,
  Standard=NA
)


for (i in seq_along(traits)) {
  t <- traits[i]
  print (paste0("Analyzing trait ",t," ......."))
  sub_df <- pheno_clean %>% filter(!is.na(.data[[t]]))
  form_blue <- as.formula(paste(t, "~ 1 + Genotype + (1|Block)"))
  fit_blue  <- lmer(form_blue, data = sub_df,
                    control = lmerControl(optimizer = "bobyqa", calc.derivs = FALSE))
  
  form_blup <- as.formula(paste(t, "~ 1 + (1|Genotype) + (1|Block)"))
  fit_blup  <- lmer(form_blup, data = sub_df,
                    control = lmerControl(optimizer = "bobyqa", calc.derivs = FALSE))
  
  fe <- fixef(fit_blue)
  names(fe) <- sub("^Genotype", "", names(fe))
  blues_list[[t]] <- data.frame(Genotype = names(fe), val = as.numeric(fe))
  colnames(blues_list[[t]])[2] <- t
  vc    <- as.data.frame(VarCorr(fit_blup))
  var_g <- vc$vcov[vc$grp == "Genotype"]
  var_e <- vc$vcov[vc$grp == "Residual"]
  
  # Standard Augmented Broad-Sense Heritability
  h2_aug <- ifelse((var_g + var_e) > 0, var_g / (var_g + var_e), 0)
  
  # Cullis Heritability: 1 - (mean PEV / Sigma2_G)
  re_geno <- ranef(fit_blup, condVar = TRUE)$Genotype
  pev     <- attr(re_geno, "postVar")[1, 1, ]
  h2_cullis <- ifelse(var_g > 0, max(0, 1 - (mean(pev) / var_g)), 0)
  sub_df$Block<-as.factor(sub_df$Block)
  sub_df$Genotype<-as.factor(sub_df$Genotype)
  
  # Compute broad sense heritability computed using the package augmentedRCBD
  ### this loop removes the missing checks to avoid errors in augmentedRCBD
  checks=c("DF0175","UG11227","H1301","H2206","H3402" ,"N6438") #no volare
  for (check in checks){
    #print (check)
    if ( (length((sub_df$Block[sub_df$Check_ID==check])))< 11){
      print (paste0("The following check is eliminated from the analysis of trait ",t
                    ,": ",check))
      checks=checks[-which(checks==check)]
    }
    
  }
  #print (checks)
  out1 <- augmentedRCBD(sub_df$Block, sub_df$Genotype, sub_df[[t]], method.comp = "lsd",
                        alpha = 0.05, group = TRUE, console = TRUE,
                        checks = checks)
  
  #### using sommer for computing heritability
  mod <- mmes(as.formula(paste(t," ~ 1")), 
              random = ~ Genotype+Block , 
              data = sub_df)
  h2_sommer <- vpredict(mod, h2 ~ V1 / (V1 + V3))
  stage2_summary$N_Genos[i]     <- length(fe)
  stage2_summary$Mean_BLUE[i]   <- round(mean(fe, na.rm = TRUE), 2)
  stage2_summary$SD_BLUE[i]     <- round(sd(fe, na.rm = TRUE), 2)
  stage2_summary$CV_BLUE_pct[i] <- round((sd(fe, na.rm = TRUE) / abs(mean(fe, na.rm = TRUE))) * 100, 2)
  stage2_summary$Sigma2_G[i]    <- round(var_g, 4)
  stage2_summary$Sigma2_E[i]    <- round(var_e, 4)
  stage2_summary$H2_Aug[i]      <- round(h2_aug, 3)
  stage2_summary$H2_Cullis[i]   <- round(h2_cullis, 3)
  stage2_summary$H2_AugmentedRCBD[i]  <- round(gva.augmentedRCBD(out1)$hBS/100,3)
  stage2_summary$H2_sommer[i]  <- round(h2_sommer$Estimate,3)
  print ("Computing heritability with heritable package")
  #H2(fit_blup , target = "Genotype",method = c("Cullis", "Oakey", "Delta", "Piepho", "Standard"))
  s=H2(fit_blup , target = "Genotype",method = c("Cullis", "Oakey", "Delta", "Piepho", "Standard"))
  heri$Cullis[i]=round(s[1],3)
  heri$Oakey[i]=round(s[2],3)
  heri$Delta[i]=round(s[3],3)
  heri$Piepho[i]=round(s[4],3)
  heri$Standard[1]=round(s[5],3)
}


write.csv(stage2_summary, "Results/1.Tomato_AugmentedRCBD_BLUEs_bis.csv", row.names = FALSE)
write.csv(heri, "Results/1.Tomato_AugmentedRCBD_BLUEs_H2_heritable.csv", row.names = FALSE)

