system("dx download -r /annovar/")

#1.Select one individual ID and transform the .bed .bim .fam format into .vcf format via PLINK.

system("dx download -r /plink/")
system("chmod 775 plink/plink")
keep_sample<-read.table("/mnt/project/Data/keep_sample.txt",header=T)
select_sample<-keep_sample[1,]
write.table(select_sample,file="select_sample_for_annovar.txt",row.names = F,quote=F)
system("dx upload select_sample_for_annovar.txt")

system("for chr in {1..22}; do
plink/plink \
--bfile /mnt/project/qc_geno/ukb_qc_chr${chr} \
--keep select_sample_for_annovar.txt \
--export vcf \
--out annovar_geno/ukb_annovar_chr${chr}
done")


#2. Transform the .vcf format into .avinput, use:
system("chmod 775 annovar/convert2annovar.pl")

system("for chr in {1..22}; do
annovar/convert2annovar.pl \
-format vcf4 -allsample -withfreq \
annovar_geno/ukb_annovar_chr${chr}.vcf > annovar_result/ukb_annovar_chr${chr}.avinput 
done")


#3. Run the ANNOVAR command: 
system("chmod 775 annovar/annotate_variation.pl")

system("for chr in {1..22}; do
annovar/annotate_variation.pl \
-buildver hg19 -geneanno \
-dbtype refGene \
annovar_result/ukb_annovar_chr${chr}.avinput annovar/humandb/ \
-out annovar_result/ukb_annovar_chr${chr}
done")


#4. Save the results to worksapce
system("dx upload -r annovar_geno")
system("dx upload -r annovar_result")
system("dx upload run_annovar.R")




