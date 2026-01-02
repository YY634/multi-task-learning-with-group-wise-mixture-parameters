#### This R script will extract identified brain subnetworks from the results of fitting Model 2.1.
#### It also formulates the identified brain networks in a readable format for subsequent steps. 
#### The output includes the gene name (network-inducing gene), which chromosome the gene lies within, the number of brain regions belonging to the network, the brain regions within this network


load("./chr18_result.RData")
L <- length(Betas); clasif <- classifications[[L]]; cnums <- cluster_nums[[L]]
rm(Betas, Blabels, classifications, cluster_nums, covariances, means, L)
load("./demo_datasets/chr18_genes.RData")

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
Bregions <- read.csv("./demo_datasets/Yeo_DK_key.csv") 
Bregion_label <- Bregions$roi
colnames(clusters) <- c("gene", "group_label", "cluster_ID", Bregion_label)

L <- dim(clusters)[1]; cluster_num <- rep(0, L)
for (i in 1:L){
  cluster_num[i] <- length(which(clusters[i,4:71] != 0)) 
}
#inddd <- which(cluster_num == max(cluster_num));  inddd = 62
clusters <- cbind(clusters[, 1:3], cluster_num, clusters[,4:71])

save(clusters, file = "./chr18_clusters.RData")



