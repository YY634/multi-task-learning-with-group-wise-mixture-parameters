library(readxl)

###### code below is for the chromosome (chr18 in UKB) which does not go through group-combinations

load("/mnt/project/Results/SNP_group/chr18_SNP_group.RData")
SNP_group<-M1; rm(M1) #make the naming consistent
all(diff(SNP_group$Pos) > 0) # answer should be TRUE
load("/mnt/project/Data/X_Rdata/chr18_X.RData"); rm(s_PC_X, n, pc_d, PC_group_sizes, PC_grouping, pc_low, pc_upp)
SNP_group <- SNP_group[which(SNP_group$SNP %in% colnames(X)),]
q == length(unique(SNP_group$Gene)) # answer should be TRUE

groups <- unique(SNP_group$Gene); M <- NULL
for (i in 1:q){
  Mt <- SNP_group[which(SNP_group$Gene == groups[i]), ];  Mt <- Mt[order(Mt$Pos), ]
  M <- rbind(M, Mt)
}
SNP_group <- M;  rm(Mt, M, groups, X)
SNP_group <- subset(SNP_group, select = -c(Pos, Func, Ref, Alt));  
save(SNP_group, file="clusters/chr18_genes_SNPs.RData")

Gene <- cbind(unique(SNP_group$Gene), c(1:q)); colnames(Gene) <- c("gene", "group_label")

load("/mnt/project/Results/final/chr18_result.RData")
L <- length(Betas); clasif <- classifications[[L]]; cnums <- cluster_nums[[L]]
rm(Betas, Blabels, classifications, cluster_nums, covariances, means, L)

clusters <- NULL
for (l in 1:q){
  cn <- cnums[l];  M <- NULL
  if (cn>0 & cn<6){
    for (k in 1:cn){
      members <- as.numeric( clasif[l,] == k )
      if (length(which(members==1)) >3){ ctp <- c(Gene[l,], k, members); M <- rbind(M, ctp) }
    }
    clusters <- rbind(clusters, M)
  }
}       
Bregion_label <- read_excel("/mnt/project/Data/brainregion_number.xlsx")
colnames(clusters) <- c("gene", "group_label", "cluster_ID", Bregion_label$roi)
save(clusters, file = "clusters/chr18_clusters.RData")

###### code above is for the chromosome (chr18 in UKB) which does not go through group-combinations
rm(list=ls())
for(chr in c(1:17,19:22)){
  print(chr)
  load(paste0("/mnt/project/Results/SNP_group/chr",chr,"_SNP_group.RData"))
  SNP_group<-M1; rm(M1) #make the naming consistent
  if(!all(diff(SNP_group$Pos) > 0)){
    stop("all(diff(SNP_group$Pos) > 0) FALSE")
  } # answer should be TRUE
  load(paste0("/mnt/project/Data/X_Rdata/chr",chr,"_X.RData")) 
  rm(n_s_PC_X, n, n_X_group_sizes, n_pc_d, n_PC_group_sizes, n_PC_grouping, n_pc_low, n_pc_upp, n_q)
  SNP_group <- SNP_group[which(SNP_group$SNP %in% colnames(n_X)),]
  if(!(q==length(unique(SNP_group$Gene)))){
    stop("q==length(unique(SNP_group$Gene)) FALSE")
  } # answer should be TRUE
  
  groups <- unique(SNP_group$Gene); M <- NULL
  for (i in 1:q){
    Mt <- SNP_group[which(SNP_group$Gene == groups[i]), ];  Mt <- Mt[order(Mt$Pos), ]
    M <- rbind(M, Mt)
  }
  SNP_group <- M;  rm(Mt, M, groups, n_X)
  SNP_group <- subset(SNP_group, select = -c(Pos, Func, Ref, Alt));  
  save(SNP_group, file=paste0("clusters/chr",chr,"_genes_SNPs.RData"))
  
  if(!all(diff(old_new_groups$old_group_label) > 0)){
    stop("all(diff(old_new_groups$old_group_label) > 0) FALSE")
  }
  if(!all(unique(SNP_group$Gene)==old_new_groups$gene)){
    stop("all(unique(SNP_group$Gene)==old_new_groups$gene) FALSE")
  }
  
  TMP <- old_new_groups[order(old_new_groups$new_group_label),] 
  rm(SNP_group, i, q, old_new_groups)
  
  load(paste0("/mnt/project/Results/final/chr",chr,"_result.RData"))  
  # the "q" in result1.RData corresponds to the "n_q" in chr09_X.RData;  the "q0" in result1.RData corresponds to the "q" in chr09_X.RData
  L <- length(Betas); clasif <- classifications[[L]]; cnums <- cluster_nums[[L]]
  rm(Betas, Blabels, classifications, cluster_nums, covariances, means, L)
  #ind <- which(cnums != 0); check <- rep(0, q); for (i in 1:q){ check[i] <- !all(clasif[i,] ==0) }; all(ind == which(check !=0)) 
  # the answer should be TRUE
  # rm(ind, check, p, i, d)
  Clstrs <- NULL
  for (l in 1:q){
    cn <- cnums[l];  M <- NULL
    if (cn>0 & cn<6){
      for (k in 1:cn){
        members <- as.numeric( clasif[l,] == k )
        if (length(which(members==1)) >3){ ctp <- c(l, k, members); M <- rbind(M, ctp) }
      }
      Clstrs <- rbind(Clstrs, M)
    }
  }
  Bregion_label <- read_excel("/mnt/project/Data/brainregion_number.xlsx")
  colnames(Clstrs) <- c("new_group_label", "cluster_ID", Bregion_label$roi); rm(Bregion_label, clasif)
  Clstrs <- as.data.frame(Clstrs)
  REF <- TMP[which(TMP$new_group_label %in% Clstrs$new_group_label),];  REF <- REF[order(REF$new_group_label),]
  rm(TMP) #clusters <- merge(REF, Clstrs, by = "new_group_label")
  
  grp_labels <- unique(REF$new_group_label);  clusters <- NULL
  for (i in 1:length(grp_labels)){
    lab <- grp_labels[i]
    ind1 <- which(REF$new_group_label == lab); ind2 <- which(Clstrs$new_group_label == lab)
    if (length(ind1) > length(ind2)){#tmp <- left_join(REF[ind1,], Clstrs[ind2,], by="new_group_label")
      U <- length(ind1); L <- length(ind2)
      tmpa <- matrix(-1, nrow =U-L, ncol = 69);  colnames(tmpa) <- colnames(Clstrs)[-1]
      tmpb <- rbind(Clstrs[ind2, -1], tmpa)  
      tmp <- as.data.frame(cbind(REF[ind1,], tmpb));  colnames(tmp) <- c(colnames(REF), colnames(tmpb))
    }else if (length(ind1) < length(ind2)) {
      U <- length(ind2); L <- length(ind1)
      tmpa <- matrix(NA, nrow = U-L, ncol = 2);  colnames(tmpa) <- c("gene", "old_group_label")
      tmpb <- rbind(REF[ind1, -3], tmpa)  
      tmp <- as.data.frame(cbind(tmpb, Clstrs[ind2,]));  colnames(tmp) <- c(colnames(tmpb), colnames(Clstrs))
    }else {
      tmp <- as.data.frame(cbind(REF[ind1, -3], Clstrs[ind2,]));  colnames(tmp) <- c("gene", "old_group_label", colnames(Clstrs))
    }
    clusters <- rbind(clusters, tmp)
  }
  
  save(clusters, file = paste0("clusters/chr",chr,"_clusters.Rdata"))
  rm(list=ls())
}

system("dx upload extract_cluster.R")

system("dx upload -r clusters")

