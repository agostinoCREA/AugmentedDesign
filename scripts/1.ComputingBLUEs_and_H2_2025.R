## this script examines the tomato field trial obtained in 2025 in Fidenza
## by PanoraSeeds

##################################### load libraries ###########
library("openxlsx")
library("ggplot2")
library("dplyr")
library("lme4")
library("emmeans")
library("augmentedRCBD")
library("here")
library("heritable")
library("sommer")
################################## script arguments ###################
file=here("data/2025/FieldMap.xlsx")  ## 20
pheno=here("data/2025//PhD_GS_25_Dataset1_Final.xlsx")
maps2025<-here("Results/2025/Field_map_2025.png")
out_table1=here("Results/2025/1.AugmentedRCBD_BLUEs_H2")
out_table2=here("Results/2025/1.AugmentedRCBD_BLUEs_H2_heritable")
REMOVE_VOLARE=FALSE

################################## main script #######################
df=read.xlsx(xlsxFile = file, sheet = "Sheet1")
colnames(df)[5]<-"Genotype"
colnames(df)[3]<-"Block"
df$Genotype<-ifelse(test = is.na(df$Genotype),yes = df$ISI.CODE,no = df$Genotype)
df$Block <- as.numeric(as.factor(df$Block))

plt=ggplot(df, aes(x = Block, y = ROW, fill = PLOT, label = Genotype)) +
  geom_tile(color = "black") +     # field layout
  geom_text(size = 2) +            # plot labels
  scale_y_reverse() +              # so row 1 is at the top
  scale_x_continuous(breaks = 1:length(unique(df$COL)), 
                     labels = unique(sort(unique(df$COL)))) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = "Blocks", y = "Row", title = "Field Map 2025")+theme_bw(base_size = 15)+theme(legend.position = "none")

ggsave(filename = maps2025,plot = plt,device = "png",width = 180,height = 100,units = "mm",dpi = 1200,limitsize = F,bg = "white")

## Merging field trial data with field map
df_p<-read.xlsx(xlsxFile = pheno, sheet = 1, na.strings = ".")
colnames(df_p)[1]<-"PLOT"
df_all=merge(df,df_p,by = "PLOT")
df_all$Type<-NA
df_all$Type[df_all$Genotype %in% c("VOLARE","H1301","N6438","H2206","DF0175","H3402","UG11227","UG11227","DF0175")]<-"Check"
df_all$Type[is.na(df_all$Type)]<-"Variety"




traits<-c("plant.wt_Kg","green.wt.Kg","red.wt.Kg","red.fruit.wt.g","pericarp.wt.g","water.wt.g","avg_firm",
          "med_firm","Brix")
colnames(df_all)[12:20]<-traits


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

for (i in seq_along(traits)){
  t <- traits[i]
  print (paste0("Analyzing trait ",t," ......."))
  sub_df <- df_all %>% filter(!is.na(.data[[t]]))
  if (isTRUE(REMOVE_VOLARE)){
    sub_df$Genotype[sub_df$Genotype=="VOLARE"]<-NA
  }
  form_blue <- as.formula(paste(t, "~ 1 + Genotype + (1|Block)"))
  fit_blue  <- lmer(formula = form_blue, data = sub_df,control = lmerControl(optimizer = "bobyqa", calc.derivs = FALSE))
  
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
  h2_aug
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
  # #print (checks)
  # out1 <- augmentedRCBD(sub_df$Block, sub_df$Genotype, sub_df[[t]], method.comp = "lsd",
  #                       alpha = 0.05, group = TRUE, console = TRUE,
  #                       checks = checks)
  # 
  # #### using sommer for computing heritability
  # mod <- mmes(as.formula(paste(t," ~ 1")), 
  #             random = ~ Genotype+Block , 
  #             data = sub_df)
  # h2_sommer <- vpredict(mod, h2 ~ V1 / (V1 + V3))
  h2_sommer<-NA
  stage2_summary$N_Genos[i]     <- length(fe)
  stage2_summary$Mean_BLUE[i]   <- round(mean(fe, na.rm = TRUE), 2)
  stage2_summary$SD_BLUE[i]     <- round(sd(fe, na.rm = TRUE), 2)
  stage2_summary$CV_BLUE_pct[i] <- round((sd(fe, na.rm = TRUE) / abs(mean(fe, na.rm = TRUE))) * 100, 2)
  stage2_summary$Sigma2_G[i]    <- round(var_g, 4)
  stage2_summary$Sigma2_E[i]    <- round(var_e, 4)
  stage2_summary$H2_Aug[i]      <- round(h2_aug, 3)
  stage2_summary$H2_Cullis[i]   <- round(h2_cullis, 3)
  #stage2_summary$H2_AugmentedRCBD[i]  <- round(gva.augmentedRCBD(out1)$hBS/100,3)
  stage2_summary$H2_AugmentedRCBD[i]<-NA
  #stage2_summary$H2_sommer[i]  <- round(h2_sommer$Estimate,3)
  stage2_summary$H2_sommer[i]<-NA
  print ("Computing heritability with heritable package")
  #H2(fit_blup , target = "Genotype",method = c("Cullis", "Oakey", "Delta", "Piepho", "Standard"))
  s=H2(fit_blup , target = "Genotype",method = c("Cullis", "Standard","Oakey","Delta")) #c("Cullis", "Oakey", "Delta", "Piepho", "Standard")
  heri$Cullis[i]=round(s[1],3)
  heri$Oakey[i]=round(s[2],3)
  heri$Delta[i]=round(s[3],3)
  #heri$Piepho[i]=round(s[4],3)
  heri$Standard[i]=round(s[5],3)
}
write.csv(stage2_summary, paste0(out_table1,"_REMOVED_VOLARE_",as.character(REMOVE_VOLARE),".csv"), row.names = FALSE)
write.csv(heri, paste0(out_table2,"_REMOVED_VOLARE_",as.character(REMOVE_VOLARE),".csv"), row.names = FALSE)


# summary(model1)
# ranef(model1)    
# vc    <- as.data.frame(VarCorr(model1))
# var_g <- vc$vcov[vc$grp == "Genotype"]
# var_e <- vc$vcov[vc$grp == "Residual"]
# 
# # Standard Augmented Broad-Sense Heritability
# h2_aug <- ifelse((var_g + var_e) > 0, var_g / (var_g + var_e), 0)
# 
# h2_aug
# s=H2(model1, target = "Genotype",method = c("Cullis", "Oakey", "Delta", "Piepho", "Standard"))
# 
# library(emmeans)
# geno <- emmeans(model2, ~CODE)
# geno=data.frame(geno)
# re=data.frame(ranef(model1)$CODE)
# re$CODE=rownames(re)
# geno=merge(re,geno,by="CODE")
# geno$adjmeans=geno$emmean+geno$X.Intercept.
# ch=c("VOLARE","H1301","N6438","H2206","DF0175","H3402","UG11227","UG11227","DF0175")
# #df_all$PlantWt[is.na(df_all$PlantWt)]<-mean(df_all$PlantWt,na.rm=T)
# #df_all=df_all[!is.na(df_all$`plant.wt.(kg)`),]
# 
# model2=augmentedRCBD(block = as.factor(df_all$COL.x),treatment = as.factor(df_all$CODE),y = as.vector(df_all$`plant.wt.(kg)`),checks =ch )
# augmentedRCBD::gva.augmentedRCBD(model2)
# 
# summary(model2)
# adj2=model2$Means ## here the adjusted means computed with the package augumenedRCBD
# write.table(x = model2$Means,file = "C:/Users/Agostino.Fricano/OneDrive - CREA/CREA/Projects&Proposals/2024_ISI_sementi_PhD/Anwar/GenomicPrediction/Practice_II/Pheno_Data/adjusted_means_red.fruit.wt.txt",quote = F,sep = "\t",row.names = F,col.names = T)
# gva.augmentedRCBD(model2)
# 
