#### This R script will compare our identified brain networks to the standard Yeo7 Networks. 
#### More specifically, the output is a table summarizing the overlap percentages between any identified brain network and the Yeo7 Networks. 
#### It will show that our identified brain networks are enriched in two or three Yeo7 Networks insteading of just randomly picking up brain regions.

library(readxl)
library(dplyr)

dk_to_yeo7_key = read.csv("./demo_datasets/Yeo_DK_key.csv", row.names = 1); dk_to_yeo7_key <- dk_to_yeo7_key[, -1]
Bregion_label <- read_excel("./demo_datasets/brainregion_number.xlsx")
tmp1 <- merge(Bregion_label, dk_to_yeo7_key, by = "roi")

yeo7 = read.csv("./demo_datasets/yeo7_NetworkNames.csv"); colnames(yeo7) = c("yeo_krienen", "NetworkNames")
yeo7_network_names <- yeo7$NetworkNames;  yeo7_network_names[4] <- "Salience"
dk_to_yeo7_key = merge(tmp1, yeo7, by = "yeo_krienen")
rm(yeo7, tmp1, Bregion_label)

yeo7composition <- list()
for (k in 1:7){ indd <- which(dk_to_yeo7_key$yeo_krienen == k);  yeo7composition[[k]] <- dk_to_yeo7_key$label[indd] }

load("./chr18_clusters.RData")
L <- dim(clusters)[1];  overlap <- matrix(0, L, 7)
for (l in 1:L){
  mem_lab <- which(clusters[l, 4:71] == 1)
  for (k in 1:7){
    ovlp <- intersect(mem_lab, yeo7composition[[k]])
    overlap[l,k] <- length(ovlp)/length(mem_lab)
  }
}
overlap_to_yeo7 <- cbind(clusters[, 1:4], overlap);  rm(l,k)
colnames(overlap_to_yeo7) <- c(colnames(clusters)[1:4], yeo7_network_names)
save(overlap_to_yeo7, file = "./chr18_overlap_to_yeo7.RData")

selected_clusters <- NULL
for (l in 1:L){
  poslength <- length(which(overlap[l,] >0)); top <- sort(overlap[l,], decreasing = TRUE)[1:4]
  if ( poslength<=4 | top[1] >= 0.5){
    selected_clusters <- rbind(selected_clusters, overlap_to_yeo7[l,])
  }
}
save(selected_clusters, file = "./chr18_selected_clusters.RData")

toshow <- selected_clusters[19,];  toshow[1] <- "NPC1"
gene_label <- rep("NPC1", 7); chr <- rep(18, 7); 
NetworkNames <- c("Visual", "Somatomotor", "Dorsal Attention", "Salience", "Limbic", "Control", "Default")
across_cluster <- as.numeric(toshow[5:11])
rois_cluster <- rep("lh_cuneus,lh_lateralorbitofrontal,lh_medialorbitofrontal,lh_parsopercularis,lh_rostralanteriorcingulate,lh_rostralmiddlefrontal,lh_superiorfrontal,rh_cuneus,rh_lateralorbitofrontal,rh_medialorbitofrontal,rh_middletemporal,rh_rostralanteriorcingulate,rh_rostralmiddlefrontal,rh_superiorfrontal", 7)
cluster_length <- rep(14, 7)

realdata_plot <- cbind(gene_label, chr, NetworkNames, across_cluster, rois_cluster, cluster_length)
#save(realdata_plot, file = "./demo_datasets/realdata_plot.RData")
#rm(list = setdiff(ls(), c("selected_clusters")))


