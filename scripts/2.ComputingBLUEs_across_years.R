# this script is based on 1.Computing BLUEs_and_H2_2026.R and computes
# BLUEs and heritability using the new dataset (all years and trials) provided
# by Anwar on 27th October 2026
library(agriutilities)
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
library(ggpubr)
########################################### script arguments ########
pheno_f<-here("data/pheno/All_years_and_trials/PhDGS_Fidenza_Foggia_data_v1.xlsx") ## data 2026
out_table=here(paste0("Results/All_years_and_trials/"))

#################################################### functions ##########
BLUE_H2<-function(trait_f, REMOVE_VOLARE,pheno, normalization=FALSE){
  if (isTRUE(REMOVE_VOLARE)){
    pheno[,7:length(pheno[1,])][pheno$Genotype=="VOLARE",]<-NA
  }
  blues_list_lmer <- list()
  blues_list_aug <-list()
  model_BLUE=list()
  model_BLUP=list()
  stage2_summary <- data.frame(
    Trait      = trait_f,
    N_Genos    = NA,
    Mean_BLUE  = NA,
    SD_BLUE    = NA,
    CV_BLUE_pct= NA,
    Sigma2_G   = NA,
    Sigma2_E   = NA,
    H2_Cullis  = NA
  )
  for (i in seq_along(trait_f)) {
    t <- trait_f[i]
    print (paste0("Analyzing trait ",t," ......."))
    sub_df <- pheno%>% filter(!is.na(.data[[t]]))
    if (isTRUE(normalization)){
      sub_df[[t]]<-scale(x =sub_df[[t]], center = T,scale = T)
    }
    form_blue <- as.formula(paste(t, "~ 1 + Genotype + year + (1|year:Block)+ (1|year:Genotype)"))
    fit_blue  <- lmer(form_blue, data = sub_df,
                      control = lmerControl(optimizer = "bobyqa", calc.derivs = FALSE))
    
    form_blup <- as.formula(paste(t, "~ 1 + (1|Genotype) + year + (1|year:Block)+ (1|year:Genotype)"))
    fit_blup  <- lmer(form_blup, data = sub_df,
                      control = lmerControl(optimizer = "bobyqa", calc.derivs = FALSE))
    ## reminder (1|year:Genotype) and (1|Genotype:year) are equivalent in lmer
    model_BLUP[[t]]<- fit_blup
    model_BLUE[[t]]<- fit_blue
    blues_list_lmer[[t]] <-data.frame(emmeans(object = fit_blue,~Genotype ))
    colnames(blues_list_lmer[[t]])[2] <- t
    vc    <- as.data.frame(VarCorr(fit_blup))
    var_g <- vc$vcov[vc$grp == "Genotype"]
    var_e <- vc$vcov[vc$grp == "Residual"]
    fe <- fixef(fit_blue)

    # Cullis Heritability: 1 - (mean PEV / Sigma2_G)
    re_geno <- ranef(fit_blup, condVar = TRUE)$Genotype
    pev     <- attr(re_geno, "postVar")[1, 1, ]
    h2_cullis <- ifelse(var_g > 0, max(0, 1 - (mean(pev) / var_g)), 0)
    sub_df$Block<-as.factor(sub_df$Block)
    sub_df$Genotype<-as.factor(sub_df$Genotype)
    stage2_summary$N_Genos[i]     <- length(fe)
    stage2_summary$Mean_BLUE[i]   <- round(mean(fe, na.rm = TRUE), 2)
    stage2_summary$SD_BLUE[i]     <- round(sd(fe, na.rm = TRUE), 2)
    stage2_summary$CV_BLUE_pct[i] <- round((sd(fe, na.rm = TRUE) / abs(mean(fe, na.rm = TRUE))) * 100, 2)
    stage2_summary$Sigma2_G[i]    <- round(var_g, 4)
    stage2_summary$Sigma2_E[i]    <- round(var_e, 4)
    stage2_summary$H2_Cullis[i]   <- round(h2_cullis, 3)
    h2_cullis <- h_cullis(
      model = fit_blup,
      genotype = "Genotype"
    )
    h2_cullis
  }
  return(list("BLUE_aug"=blues_list_aug,"BLUE_lmer"=blues_list_lmer,"model_BLUE"=model_BLUE,"model_BLUP"=model_BLUP,"HeritabilityStand"=stage2_summary))
}
model_diagnostics<-function(res,REMOVE_VOLARE,df_all){
  plt<-list()
  qq<-list()
  var_com<-list()
  random_effects<-list()
  my_data_clean<-list()
  observation_outliers<-list()
  genotype_outliers<-list()
  for (trait in names(res$model_BLUP)) {
    df_diag <- data.frame(
      fitted = fitted(res$model_BLUP[[trait]]),
      resid = resid(res$model_BLUP[[trait]])
    )
    plot(res$model_BLUP[[trait]], which = 1)
    plt[[trait]]<-ggplot(df_diag, aes(x = fitted, y = resid)) +
      geom_point(alpha = 0.6) +
      geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
      geom_smooth(method = "loess", se = FALSE, color = "blue") +
      labs(
        x = "Fitted Values)",
        y = "Residuals",
        title = trait
      ) +
      theme_minimal()
    qq[[trait]]<-ggplot(df_diag, aes(sample = resid)) +
      stat_qq(alpha = 0.6) +
      stat_qq_line(color = "red", linetype = "dashed") +
      labs(
        title = "Q-Q Plot of residuals",
        x = "Theoretical distribution",
        y = "Observed distribution"
      ) +
      theme_minimal()
    
    var_com[["var"]][[trait]]<-as.data.frame(VarCorr(res$model_BLUP[[trait]]))
    
    var_com[["singular"]][[trait]]<-isSingular(res$model_BLUP[[trait]], tol = 1e-4)
    random_effects[[trait]]<- ranef(res$model_BLUP[[trait]])$Genotype
    df_all2=df_all
    if (isTRUE(REMOVE_VOLARE)){
      df_all2[,7:length(df_all2[1,])][df_all2$Genotype=="VOLARE",]<-NA
    }
    my_data_clean[[trait]] <- df_all2[,c("year","Genotype","Block",trait)] %>%
      filter(!is.na(!!sym(trait)))  %>%
      mutate(
        # standardize residual values (ration betwenn residual and its standard deviation)
        std_resid = residuals(res$model_BLUP[[trait]], type = "pearson"),
        #threshold for outliers = 3)
        is_outlier_obs = abs(std_resid) > 3
      )
    genotype_outliers[[trait]] <- random_effects[[trait]] %>%
      mutate(
        Genotype_ID = rownames(.),
        z_score = (`(Intercept)` - mean(`(Intercept)`)) / sd(`(Intercept)`)
      ) %>%
      filter(abs(z_score) > 3) # threshold used is 3
    
    observation_outliers[[trait]] <- my_data_clean[[trait]] %>% 
      filter(is_outlier_obs == TRUE)
    
  }
  return(list("residuals"=plt,"qqplot"=qq,"variance_components"=var_com,
              "outliers"=observation_outliers,"outliers_geno"=genotype_outliers))
}


######################################## main script ##################
### these are the priority traits coded in the new file:
traits<-c("LBF","Plant_g","Red_g","SS_ratio","Maturity_pct","DAS","Firm_avg",
          "Brix","Peri_dry_g","SPAD")

pheno_raw=list()
pheno_df=list()
pheno_clean=list()

for (tr in c("FD25","FD26")){
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

pheno_clean$FD25$year<-"2025"
pheno_clean$FD26$year<-"2026"
df_all=rbind(pheno_clean$FD25,pheno_clean$FD26)
res=BLUE_H2(trait_f = traits,REMOVE_VOLARE = TRUE,pheno = df_all,normalization = F)
#res2=BLUE_H2(trait_f = traits,REMOVE_VOLARE = FALSE,pheno = df_all)

dia<-model_diagnostics(res = res, REMOVE_VOLARE = T,df_all = df_all)
for (trait in names(res$model_BLUP)){
  write.table(x = res$BLUE_lmer[[trait]],file = paste0("Results/All_years_and_trials/model_diagnostics/",trait,"_BLUEs.txt"),sep = "\t",col.names = T,dec = ".", row.names = F)
}



## Outlier a livello di Singola Osservazione (Residui) :dia$outliers 
## Outlier a livello di Genotipo :dia$outliers_geno
plt<-list()
for (trait in names(res$model_BLUP)){
  plt[[trait]]<-ggarrange(dia$residuals[[trait]],dia$qqplot[[trait]],nrow = 2)
  ggsave(filename = paste0(trait,"_normalized_QQ_plots.png"),device = "png",plot = plt[[trait]],width = 90,height = 90,units = "mm",dpi = 1200,limitsize = F,bg = "white",path = "Results/All_years_and_trials/model_diagnostics/")
  write.table(x = dia$variance_components$var[[trait]],file = paste0("Results/All_years_and_trials/model_diagnostics/",trait,"_normalized_var_com.txt"),sep = "\t",col.names = T,dec = ".", row.names = F)
  write.table(x = dia$outliers[[trait]],file = paste0("Results/All_years_and_trials/model_diagnostics/",trait,"_normalized_outliers.txt"),sep = "\t",col.names = T,dec = ".", row.names = F)
  write.table(x = dia$outliers_geno[[trait]],file = paste0("Results/All_years_and_trials/model_diagnostics/",trait,"_normalized_outliers_geno.txt"),sep = "\t",col.names = T,dec = ".", row.names = F)
  }


rs_var<-do.call(rbind.data.frame, dia$variance_components$var)



