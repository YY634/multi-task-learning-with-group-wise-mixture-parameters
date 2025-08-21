#preprocess the phenotype and covariate
df<-read.csv("/mnt/project/imaging_covariate_participant.csv")
df_sub<-df[,c(1,3:72,75,76,78:87)]
sum(is.na(df_sub)) #p26904_i2 has one NA
sum(df_sub=="",na.rm=T) #p1707_i0 has 9 ""
summary(df_sub)

df_sub[df_sub==""]<-NA
df_sub<-df_sub[complete.cases(df_sub),]
df_sub<-df_sub[df_sub$p1707_i0!="Prefer not to answer",] #remove 9 participants
summary(df_sub)
sum(is.na(df_sub))
sum(df_sub=="")
table(df_sub$p31)
table(df_sub$p1707_i0)
sum(table(df_sub$p22000))

write.csv(df_sub,file="Final_phenotype_covariate.csv",row.names = F,quote = F)
system("dx upload Final_phenotype_covariate.csv")

keep_sample<-df_sub[,c(1,1)]
colnames(keep_sample)<-c("FID","IID")
write.table(keep_sample,file="keep_sample.txt",row.names = F,quote=F)
system("dx upload keep_sample.txt")

#preprocess the genotype
system("dx download -r /plink/")
system("chmod 775 plink/plink")

system("for chr in {1..22}; do
plink/plink \
--bfile /mnt/project/Bulk/'Genotype Results'/'Genotype calls'/ukb22418_c${chr}_b0_v2 \
--keep keep_sample.txt \
--geno 0.1 \
--hwe 0.000001 \
--maf 0.01 \
--make-bed \
--mind 0.1 \
--out qc_geno/ukb_qc_chr${chr}
done")

system("dx upload -r qc_geno")


