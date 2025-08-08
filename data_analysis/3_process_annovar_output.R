# annotate_variation.pl -geneanno -buildver hg19 -dbtype refGene XXXX.avinput humandb/ -out XXXX 
# produces three output files: .variant_function/  .exonic_variant_function/  .log

for(chr in 1:22){
  # read in the .variant_function file that is output by ANNOVAR
  vf <- read.table(paste0("/mnt/project/Results/annovar_result/ukb_annovar_chr",chr,".variant_function"), sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  colnames(vf)[1:7] <- c("Func", "Gene", "Chr", "Start", "End", "Ref", "Alt")
  vf$Chr <- gsub("^chr", "", vf$Chr); vf <- subset(vf, nchar(Ref) == 1 & nchar(Alt) == 1)
  
  #remove the intergenic SNPs
  use <- vf[vf$Func != "intergenic", ];  rm(vf)
  
  # read in the .bim file. The .bim file should be the .bim you used to generate the .vcf file (input into ./plink convert to vcf format). 
  BIM <- read.table(paste0("/mnt/project/Data/qc_geno/ukb_qc_chr",chr,".bim"), sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  colnames(BIM) <- c("Chr", "SNP", "cM", "Pos", "Allele1", "Allele2")
  BIM$Chr <- as.character(BIM$Chr)
  
  # get overlap part between the BIM and vf
  use$key <- paste(use$Chr, use$Start, sep = ":"); U_uniq <- use[!duplicated(use$key) & !duplicated(use$key, fromLast = TRUE), ]
  BIM$key <- paste(BIM$Chr, BIM$Pos, sep = ":");  B_uniq <- BIM[!duplicated(BIM$key) & !duplicated(BIM$key, fromLast = TRUE), ];  rm(BIM)
  merged <- merge(U_uniq, B_uniq, by = "key") 
  M <- merged[, c( "Chr.x", "Pos", "SNP", "Func", "Gene", "Ref", "Alt")]
  colnames(M) <- c( "Chr", "Pos", "SNP", "Func", "Gene", "Ref", "Alt")
  M <- M[order(M$Chr, M$Pos), ]
  rm(B_uniq, U_uniq, merged, use)
  
  # deal with SNPs that belong to two genes. 
  gene_clean <- gsub("\\s*\\(.*?\\)", "", M$Gene); M$Gene <- gene_clean; rm(gene_clean)
  cind <- grep(",", M$Gene);  sind <- grep(";", M$Gene);  oind <- intersect(cind, sind)
  cind <- cind[!cind %in% oind]; sind <- sind[!sind %in% oind]
  # remove SNPs that belong to too many genes, those annotations are probably not correct.
  if(length(oind)==0){
    M1 <- M
  }else{
    M1 <- M[-oind,]
  }
  
  if(length(unique(M1$SNP)) == dim(M1)[1]){
    check1 <- M$Gene[cind];  check2 <- M$Gene[sind]
    print(paste0("chr:",chr,", #snp:",nrow(M1)))
    file_name <- paste0("SNP_group/chr", chr,"_SNP_group.RData")
    save(M1, file = file_name)
  }else{
    stop()
  } 
  
}

#Upload the results to worksapce
system("dx upload -r SNP_group")
system("dx upload process_annovar_output_revise.R")
system("dx upload process_annovar_output.R")

"chr:1, #snp:26595"
"chr:2, #snp:23477"
"chr:3, #snp:22240"
"chr:4, #snp:17198"
"chr:5, #snp:17071"
"chr:6, #snp:20911"
"chr:7, #snp:17695"
"chr:8, #snp:15538"
"chr:9, #snp:13188"
"chr:10, #snp:16230"
"chr:11, #snp:16362"
"chr:12, #snp:15524"
"chr:13, #snp:9217"
"chr:14, #snp:9446"
"chr:15, #snp:10540"
"chr:16, #snp:11396"
"chr:17, #snp:12523"
"chr:18, #snp:8353"
"chr:19, #snp:11113"
"chr:20, #snp:8247"
"chr:21, #snp:4386"
"chr:22, #snp:5945"


#use$key <- paste(use$Chr, use$Start, use$Ref, use$Alt, sep = ":")
#BIM$key1 <- paste(BIM$Chr, BIM$Pos, BIM$Allele1, BIM$Allele2, sep = ":")
#BIM$key2 <- paste(BIM$Chr, BIM$Pos, BIM$Allele2, BIM$Allele1, sep = ":")  # account for strand flip
#U1 <- use[use$key %in% BIM$key1, ];  U1_uniq <- U1[!duplicated(U1$key) & !duplicated(U1$key, fromLast = TRUE), ]
#B1 <- BIM[BIM$key1 %in% use$key,];  B1_uniq <- B1[!duplicated(B1$key1) & !duplicated(B1$key1, fromLast = TRUE), ]
#M1 <- merge(U1_uniq, B1_uniq, by="key")
#U2 <- use[use$key %in% BIM$key2, ];  U2_uniq <- U2[!duplicated(U2$key) & !duplicated(U2$key, fromLast = TRUE), ] 
#B2 <- BIM[BIM$key2 %in% use$key,];  B2_uniq <- B2[!duplicated(B2$key2) & !duplicated(B2$key2, fromLast = TRUE), ]
#M2 <- merge(U2_uniq, B2_uniq, by="key")

