# this script is based on 1.Computing BLUEs_and_H2_2026.R and computes
# BLUEs and heritability using the new dataset (all years and trials) provided
# by Anwar on 27th October 2026

library(here)
library(tidyverse)
library(lme4)
library(rrBLUP)
library(data.table)
library(ggplot2)
library(augmentedRCBD)
library(sommer)
library(heritable)
library(openxlsx)
library(emmeans)
########################################### script arguments ########
pheno_f<-here("data/pheno/All_years_and_trials/PhDGS_Fidenza_Foggia_data_v1.xlsx") ## data 2026
out_table=here(paste0("Results/All_years_and_trials/"))

#################################################### function ##########
BLUE_H2<-function(trait_f, REMOVE_VOLARE,pheno){
  if (isTRUE(REMOVE_VOLARE)){
    pheno[,7:length(pheno[1,])][pheno$Genotype=="VOLARE",]<-NA
  }
  blues_list_lmer <- list()
  blues_list_aug <-list()
  stage2_summary <- data.frame(
    Trait      = trait_f,
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
    Trait=trait_f,
    Cullis=NA,
    Oakey=NA,
    Delta=NA,
    Piepho=NA,
    Standard=NA
  )
  
  
  for (i in seq_along(trait_f)) {
    t <- trait_f[i]
    print (paste0("Analyzing trait ",t," ......."))
    sub_df <- pheno%>% filter(!is.na(.data[[t]]))
    form_blue <- as.formula(paste(t, "~ 1 + Genotype + (1|Block)"))
    fit_blue  <- lmer(form_blue, data = sub_df,
                      control = lmerControl(optimizer = "bobyqa", calc.derivs = FALSE))
    
    form_blup <- as.formula(paste(t, "~ 1 + (1|Genotype) + (1|Block)"))
    fit_blup  <- lmer(form_blup, data = sub_df,
                      control = lmerControl(optimizer = "bobyqa", calc.derivs = FALSE))
    fe <- fixef(fit_blue)
    names(fe) <- sub("^Genotype", "", names(fe))
    #blues_list[[t]] <- data.frame(Genotype = names(fe), val = as.numeric(fe))
    blues_list_lmer[[t]] <-data.frame(emmeans(object = fit_blue,~Genotype ))
    colnames(blues_list_lmer[[t]])[2] <- t
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
    checks=c("DF0175","UG11227","H1301","H2206","H3402" ,"N6438","VOLARE") #no volare
    if (isTRUE(REMOVE_VOLARE)){
      checks<-checks[-which(checks=="VOLARE")]
    }
    for (check in checks){
      #print (check)
      if ( (length((sub_df$Block[sub_df$Check_ID==check]))) <= length(unique(as.character(sub_df$Block)))){
        print (paste0("The following check is eliminated from the analysis of trait ",t
                      ,": ",check))
        checks=checks[-which(checks==check)]
      }
      
    }
    print (checks)
    out1 <- augmentedRCBD(sub_df$Block, sub_df$Genotype, sub_df[[t]], method.comp = "lsd",
                          alpha = 0.05, group = FALSE, console = TRUE,
                          checks = checks)
    
    blues_list_aug[[t]]<-out1$Means
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
    if (t == "Peri_dry_g") {
      s=H2(fit_blup , target = "Genotype",method = c("Standard"))
    } else {
      s=H2(fit_blup , target = "Genotype",method = c("Cullis", "Oakey", "Delta", "Standard"))
    }
    
    heri$Cullis[i]=round(s[1],3)
    heri$Oakey[i]=round(s[2],3)
    heri$Delta[i]=round(s[3],3)
    #heri$Piepho[i]=round(s[4],3)
    heri$Standard[1]=round(s[4],3)
  }
  return(list("BLUE_aug"=blues_list_aug,"BLUE_lmer"=blues_list_lmer,"HeritabilityStand"=stage2_summary,
              "HeritabiliyHeritable"=heri))
}




######################################## main script ##################
### these are the priority traits coded in the new file:
traits<-c("LBF","Plant_g","Red_g","SS_ratio","Maturity_pct","DAS","Firm_avg",
          "Brix","Peri_dry_g","SPAD")

pheno_raw=list()
pheno_df=list()
pheno_clean=list()
for (tr in c("FD25","FD26","FG26")){
  pheno_raw[[tr]]<-read.xlsx(xlsxFile = pheno_f,sheet = tr)
  length(traits)
  traits<-traits[traits %in% colnames(pheno_raw[[tr]])]
  # Clean strings and prepare Augmented RCBD grouping factors
  # Note: Year is ignored (all 2026); Column (A-K) and Block (1-11) are identical
  pheno_df[[tr]] <- pheno_raw[[tr]] %>%
    dplyr::select(ID, Genotype, Type, Col, Row, Block, all_of(traits)) %>%
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
  pheno_clean[[tr]] <- pheno_df[[tr]] %>%
    mutate(
      Plant_g = ifelse(Plant_g <= 0, NA, Plant_g),
      Red_g = ifelse(Red_g <= 0, NA, Red_g)
    )
  traits<-c("LBF","Plant_g","Red_g","SS_ratio","Maturity_pct","DAS","Firm_avg",
            "Brix","Peri_dry_g","SPAD")

}

FD25_no_volare=BLUE_H2(trait_f = traits[6],REMOVE_VOLARE = TRUE,pheno = pheno_clean$FD25) # fail
FD25=BLUE_H2(trait_f = traits,REMOVE_VOLARE = FALSE,pheno = pheno_clean$FD25) # fail  
FD26_no_volare=BLUE_H2(trait_f = traits[6],REMOVE_VOLARE = TRUE,pheno = pheno_clean$FD26)
FD26=BLUE_H2(trait_f = traits,REMOVE_VOLARE = FALSE,pheno = pheno_clean$FD26)

FG26_no_volare=BLUE_H2(trait_f = traits,REMOVE_VOLARE = TRUE,pheno = pheno_clean$FG26)
FG26=BLUE_H2(trait_f = traits,REMOVE_VOLARE = FALSE,pheno = pheno_clean$FG26)


write.xlsx(x = FD25_no_volare[c("HeritabilityStand","HeritabiliyHeritable")], file = paste0(out_table,"/","H2_FD25_no_VOLARE_all",".xlsx"))
write.xlsx(x = FD26_no_volare[c("HeritabilityStand","HeritabiliyHeritable")], file = paste0(out_table,"/","H2_FD26_no_VOLARE_all",".xlsx"))

write.xlsx(x = FD25[c("HeritabilityStand","HeritabiliyHeritable")], file = paste0(out_table,"/","H2_FD25_with_VOLARE_all",".xlsx"))
write.xlsx(x = FD26[c("HeritabilityStand","HeritabiliyHeritable")], file = paste0(out_table,"/","H2_FD26_with_VOLARE_all",".xlsx"))


for (t in traits){
  write.xlsx(x = FD25_no_volare$BLUE_aug[[t]], file = paste0(out_table,"/",t,"_BLUE_Aug_FD25_no_VOLARE_all",".xlsx"))
  write.xlsx(x = FD26_no_volare$BLUE_lmer[[t]], file = paste0(out_table,"/",t,"_BLUE_lmer_FD26_no_VOLARE_all",".xlsx"))
  write.xlsx(x = FD25$BLUE_aug[[t]], file = paste0(out_table,"/",t,"_BLUE_Aug_FD25_with_VOLARE_all",".xlsx"))
  write.xlsx(x = FD26$BLUE_lmer[[t]], file = paste0(out_table,"/",t,"_BLUE_lmer_FD25_with_VOLARE_all",".xlsx"))
  
  
}

# sub_df$TestGenotype <- ifelse(
#   sub_df$Check_ID == "Test_Population",
#   as.character(sub_df$Genotype),
#   NA
# )
# sub_df$CheckGenotype <- ifelse(
#   sub_df$Check_ID != "Test_Population",
#   as.character(sub_df$Genotype),
#   NA
# )
# sub_df$TestGenotype <- factor(sub_df$TestGenotype)
# sub_df$CheckGenotype <- factor(sub_df$CheckGenotype)
# sub_df <- sub_df %>%
#   mutate(
#     Genotype = factor(Genotype),
#     Block = factor(Block),
#     Type = ifelse(
#       Genotype %in% out1$Details$`Check treatments`,
#       "Check",
#       "Test"
#     )
#   )
# fit_test <- lmer(
#   Peri_dry_g ~ Type + 
#     (1 | Block) +
#     (1 | TestGenotype),
#   data = sub_df,
#   REML = TRUE
# )
# fit_test <- lmer(
#   Peri_dry_g ~ 0 + Genotype +
#     (1 | Block) +
#     (1 | TestGenotype),
#   data = sub_df,
#   REML = TRUE
# )
# fit_test <- lmer(
#   Peri_dry_g ~ CheckGenotype +
#     (1 | Block) +
#     (1 | TestGenotype),
#   data = sub_df,
#   REML = TRUE,
#   subset = Type == "Test" | Type == "Check"
# )
