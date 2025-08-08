library(readxl)
library(dplyr)

path = "/mnt/project/Data/"
dk_to_yeo7_key = read.csv(paste0(path, "Yeo_DK_key.csv"), row.names = 1); dk_to_yeo7_key <- dk_to_yeo7_key[, -1]
Bregion_label <- read_excel(paste0(path, "brainregion_number.xlsx"))
#load("~/project/JASA_case/Y.RData");  brain_names = colnames(Y)
#brain_names = paste0(sapply(sapply(brain_names, function(x) strsplit(x, "[_]")[[1]][1]), function(xx) ifelse(xx == "L", "lh", "rh")), "_", sapply(brain_names, function(x) strsplit(x, "[_]")[[1]][2]))
#all(brain_names == Bregion_label$roi) should be TRUE
tmp1 <- merge(Bregion_label, dk_to_yeo7_key, by = "roi")

yeo7 = read.csv(paste0(path, "yeo7_NetworkNames.csv")); colnames(yeo7) = c("yeo_krienen", "NetworkNames")
yeo7_network_names <- yeo7$NetworkNames;  yeo7_network_names[4] <- "Salience"
dk_to_yeo7_key = merge(tmp1, yeo7, by = "yeo_krienen")
rm(yeo7, tmp1, Bregion_label)
# "label" in dk_to_yeo7_key is the index for brain regions in our results

yeo7composition <- list()
for (k in 1:7){
  indd <- which(dk_to_yeo7_key$yeo_krienen == k);  yeo7composition[[k]] <- dk_to_yeo7_key$label[indd]
}

###### code below is for the chromosome (chr18 in UKB) which does not go through group-combinations
load("/mnt/project/Results/clusters/chr18_clusters.RData")
L <- dim(clusters)[1];  overlap <- matrix(0, L, 7)
for (l in 1:L){
  mem_lab <- which(clusters[l, 4:71] == 1)
  for (k in 1:7){
    ovlp <- intersect(mem_lab, yeo7composition[[k]])
    overlap[l,k] <- length(ovlp)/length(mem_lab)
  }
}
overlap_to_yeo7 <- cbind(clusters[, 1:3], overlap);  rm(l,k)
colnames(overlap_to_yeo7) <- c(colnames(clusters)[1:3], yeo7_network_names)
save(overlap_to_yeo7, file = "overlap/chr18_overlap_to_yeo7.RData")

selected_clusters <- NULL
for (l in 1:L){
  poslength <- length(which(overlap[l,] >0)); top <- sort(overlap[l,], decreasing = TRUE)[1:4]
  if ( poslength<=4 | top[1] >= 0.5){
    selected_clusters <- rbind(selected_clusters, overlap_to_yeo7[l,])
  }
}
save(selected_clusters, file = "overlap/chr18_selected_clusters.RData")

###### code above is for the chromosome (chr18 in UKB) which does not go through group-combinations
rm(list = setdiff(ls(), c("dk_to_yeo7_key", "yeo7composition", "mart", "path", "yeo7_network_names")))

for(chr in c(1:17,19:22)){
  print(chr)
  load(paste0("/mnt/project/Results/clusters/chr",chr,"_clusters.Rdata"))
  L <- dim(clusters)[1];  overlap <- matrix(0, L, 7)
  for (l in 1:L){
    if (all(clusters[l, 5:72]== -1)){ overlap[l,] <- rep(-1, 7) }else{
      mem_lab <- which(clusters[l, 5:72] == 1)
      for (k in 1:7){
        ovlp <- intersect(mem_lab, yeo7composition[[k]])
        overlap[l,k] <- length(ovlp)/length(mem_lab)
      }
    }
  }
  overlap_to_yeo7 <- cbind(clusters[, 1:4], overlap)
  colnames(overlap_to_yeo7) <- c(colnames(clusters)[1:4], yeo7_network_names)
  save(overlap_to_yeo7, file = paste0("overlap/chr",chr,"_overlap_to_yeo7.RData"))
  
  
  selected_clusters <- NULL;  check <- NULL
  for (l in 1:L){
    if (all(overlap[l,]== -1)){ check <- c(check, l) }else{
      poslength <- length(which(overlap[l,] >0)); top <- sort(overlap[l,], decreasing = TRUE)[1:4]
      if ( poslength<=4 | top[1] >= 0.5){
        selected_clusters <- rbind(selected_clusters, overlap_to_yeo7[l,])
      }
    }
  }
  selected_clusters <- as.data.frame(selected_clusters)
  check1 <- overlap_to_yeo7$new_group_label[check]; comn <- intersect(check1, selected_clusters$new_group_label)
  if (length(comn) >0){
    check2 <- NULL
    for (l in 1:L){
      if ( (l %in% check) & (overlap_to_yeo7$new_group_label[l] %in% comn) ){ check2 <- c(check2, l) }
    }
    if (length(check2) >0){ 
      Mtp <- as.data.frame(rbind(selected_clusters, overlap_to_yeo7[check2,])) 
      selected_clusters <- Mtp[order(Mtp$new_group_label),]
    }
  }
  
  save(selected_clusters, file = paste0("overlap/chr",chr,"_selected_clusters.RData"))
  
  rm(list = setdiff(ls(), c("dk_to_yeo7_key", "yeo7composition", "mart", "path", "yeo7_network_names")))
}


system("dx upload cluster_yeo7_overlap.R")

system("dx upload -r overlap")

