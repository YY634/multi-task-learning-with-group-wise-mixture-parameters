df<-read.csv("/mnt/project/Data/Final_phenotype_covariate.csv")
pheno<-df[,c(1:69)]
colnames(pheno)<-c("ID",
                   "L_bankssts", "R_bankssts","L_caudalanteriorcingulate", "R_caudalanteriorcingulate",
                   "L_caudalmiddlefrontal", "R_caudalmiddlefrontal","L_cuneus", "R_cuneus",
                   "L_entorhinal", "R_entorhinal","L_frontalpole", "R_frontalpole",
                   "L_fusiform", "R_fusiform","L_inferiorparietal", "R_inferiorparietal",
                   "L_inferiortemporal", "R_inferiortemporal","L_insula", "R_insula",
                   "L_isthmuscingulate", "R_isthmuscingulate","L_lateraloccipital", "R_lateraloccipital",
                   "L_lateralorbitofrontal", "R_lateralorbitofrontal","L_lingual", "R_lingual",
                   "L_medialorbitofrontal", "R_medialorbitofrontal","L_middletemporal", "R_middletemporal",
                   "L_paracentral", "R_paracentral","L_parahippocampal", "R_parahippocampal",
                   "L_parsopercularis", "R_parsopercularis","L_parsorbitalis", "R_parsorbitalis",
                   "L_parstriangularis", "R_parstriangularis","L_pericalcarine", "R_pericalcarine",
                   "L_postcentral", "R_postcentral","L_posteriorcingulate", "R_posteriorcingulate",
                   "L_precentral", "R_precentral","L_precuneus", "R_precuneus",
                   "L_rostralanteriorcingulate", "R_rostralanteriorcingulate","L_rostralmiddlefrontal", "R_rostralmiddlefrontal",
                   "L_superiorfrontal", "R_superiorfrontal","L_superiorparietal", "R_superiorparietal",
                   "L_superiortemporal", "R_superiortemporal",  "L_supramarginal", "R_supramarginal",
                   "L_transversetemporal", "R_transversetemporal","L_temporalpole","R_temporalpole")


hcp_column_order <- c("ID","L_bankssts", "L_caudalanteriorcingulate", "L_caudalmiddlefrontal",
                      "L_cuneus", "L_entorhinal", "L_fusiform", "L_inferiorparietal",
                      "L_inferiortemporal", "L_isthmuscingulate", "L_lateraloccipital", "L_lateralorbitofrontal",
                      "L_lingual", "L_medialorbitofrontal", "L_middletemporal", "L_parahippocampal",
                      "L_paracentral", "L_parsopercularis", "L_parsorbitalis", "L_parstriangularis",
                      "L_pericalcarine", "L_postcentral", "L_posteriorcingulate", "L_precentral",
                      "L_precuneus", "L_rostralanteriorcingulate", "L_rostralmiddlefrontal", "L_superiorfrontal",
                      "L_superiorparietal", "L_superiortemporal", "L_supramarginal", "L_frontalpole",
                      "L_temporalpole", "L_transversetemporal", "L_insula",
                      "R_bankssts", "R_caudalanteriorcingulate", "R_caudalmiddlefrontal", "R_cuneus",
                      "R_entorhinal", "R_fusiform", "R_inferiorparietal", "R_inferiortemporal",
                      "R_isthmuscingulate", "R_lateraloccipital", "R_lateralorbitofrontal", "R_lingual",
                      "R_medialorbitofrontal", "R_middletemporal", "R_parahippocampal", "R_paracentral",
                      "R_parsopercularis", "R_parsorbitalis", "R_parstriangularis", "R_pericalcarine",
                      "R_postcentral", "R_posteriorcingulate", "R_precentral", "R_precuneus",
                      "R_rostralanteriorcingulate", "R_rostralmiddlefrontal", "R_superiorfrontal", "R_superiorparietal",
                      "R_superiortemporal", "R_supramarginal", "R_frontalpole", "R_temporalpole",
                      "R_transversetemporal", "R_insula")

# Reorder columns of data frame
pheno <- pheno[ , hcp_column_order]
sum(colnames(pheno)==hcp_column_order)

pheno_scale<-pheno
pheno_scale[,2:69]<-scale(pheno[,2:69],center=T,scale=F)

write.csv(pheno_scale,file="Final_phenotype_scaled.csv",row.names = F,quote = F)
system("dx upload Final_phenotype_scaled.csv")

#check the distribution between HCP and UKB
load("original_brain_regVs.RData")
c1 <- rgb(0.1, 0.7, 0.2, alpha = 0.4);  c2 <- rgb(1, 0.6, 0, alpha = 0.5)
par(mfrow=c(3,2))
for(i in 2:4){
  hist(pheno[,i], breaks=20)                    
  hist(d[,i], breaks=20)   
}


system("dx upload process_Y.R")



