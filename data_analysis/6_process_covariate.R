df<-read.csv("/mnt/project/Data/Final_phenotype_covariate.csv")
cov<-df[,c(1,70:83)]
colnames(cov)<-c("ID","age","sex","handedness","batch_group",paste0("PC",1:10))

table(cov$sex)
cov$sex<-ifelse(cov$sex == "Male", 1, 0)
table(cov$sex)

table(cov$handedness)
cov$handedness <- ifelse(cov$handedness == "Left-handed", -1,
                  ifelse(cov$handedness == "Right-handed", 1, 0))
table(cov$handedness)

sum(grepl("UKBiLEVEAX_", cov$batch_group))
sum(grepl("Batch_b", cov$batch_group))
cov$batch_group <- ifelse(grepl("UKBiLEVEAX_", cov$batch_group), 0, 1)
table(cov$batch_group)

cov_scale<-cov
cov_scale[,c(2,6:15)]<-scale(cov[,c(2,6:15)],center=T,scale=T) 

write.csv(cov_scale,file="Final_covariate_scaled.csv",row.names = F,quote = F)
system("dx upload Final_covariate_scaled.csv")
system("dx upload process_covariate.R")
