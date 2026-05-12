
library(readxl)
library(mclust)

df <- read_excel("~/Documents/papers/groupwise mixture/JASA_casestudy/heritability_based_clustering/68regions_heritability_to_import.xlsx")
df1 <- df[, c(3, 5, 7:28)]

for (i in 1:22) {
  coln <- sprintf("chr_%02d", i)
  df1[[coln]] <- as.numeric(sub("\\(.*", "", as.character(df1[[coln]])))
}

df1$chr_03[11] <- 0.2794
df1$chr_10[11] <- 0.2414
df1$chr_19[12] <- -0.1255
df1$chr_06[13] <- 0.1307

colnames(df1)[1] <- "label"
#save(df1, file = "~/Documents/papers/groupwise mixture/JASA_casestudy/heritability_based_clustering/global_local_heritability.RData")

region_num <- read_excel("~/Documents/papers/groupwise mixture/JASA_casestudy/real_data_analysis/brainregion_number.xlsx")
colnames(region_num) <- c("num_label", "label")

temp_dat <- merge(df1, region_num, by = "label")
dat <- temp_dat[order(temp_dat$num_label), ]
rm(temp_dat, df1)

#col_names <- grep("^chr_", names(df1), value = TRUE)
#df1[col_names] <- lapply( col_names, function(nm) { as.numeric(sub("\\(.*", "", as.character(df1[[nm]]))) } ) 
global <- Mclust(dat$Heritability)
# global$G records the number of clusters;   global$classification records the cluster labels
local <- Mclust(dat[, 3:24])
cluster <- local$classification
dat <- cbind(cluster, dat);  rm(cluster, local, global)
dat <- dat[, c("cluster", "num_label", "label")]
save(dat, file = "~/Documents/papers/groupwise mixture/JASA_casestudy/heritability_based_clustering/global_local_heritability.RData")




#### create the collection of clusters for plotting 
plot_dat <- dat
dk_to_yeo7_key = read.csv("~/Documents/papers/groupwise mixture/JASA_casestudy/real_data_analysis/Yeo_DK_key.csv", row.names = 1)
yeo7 = read.csv("~/Documents/papers/groupwise mixture/JASA_casestudy/real_data_analysis/yeo7_NetworkNames.csv")
colnames(yeo7) = c("yeo_krienen", "NetworkNames")
dk_to_yeo7_key = merge(dk_to_yeo7_key, yeo7, by = "yeo_krienen")
rm(yeo7)  
colnames(dk_to_yeo7_key)[3] = "label"
dk_to_yeo7_key$id <- NULL
TMP <- merge(plot_dat, dk_to_yeo7_key, by="label")

library(dplyr)
library(RColorBrewer)
library(ggseg)
library(ggpubr)

yeo7_plot = dk_to_yeo7_key %>% 
  mutate(label_id = sapply(label , function(x) strsplit(as.character(x), "[_]")[[1]][2])) %>%
  # group_by(label_id) %>%
  ggseg(mapping=aes(fill=NetworkNames),atlas = dk, 
        position="stacked", colour="black") +
  labs(title = "Yeo7(2011) Networks" ) +
  guides(fill=guide_legend(title="Yeo7 Networks color scheme"))+
  theme(legend.justification = c(1, 0),
        legend.position = "bottom",
        legend.text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = brewer.pal(7, "Set3"), na.value = "white")


yeo7_color = data.frame(color = brewer.pal(7, "Set3"), NetworkNames = names(table(dk_to_yeo7_key$NetworkNames)), stringsAsFactors = F)

PLOT <- TMP %>% left_join(yeo7_color, by = "NetworkNames")

for (i in 1:4){
  C <- PLOT[PLOT$cluster==i, ] 
  colnames(C) <- c("label", "cluster", "num_label", "yeo_krienen", "NetworkNames", "color")      
  cluster_plot = C %>%
    ggseg(mapping=aes(fill= NetworkNames), atlas = dk, position="stacked", colour="black") +
    labs(title = paste0("Brain subnetwork ", i)) +
    theme(legend.justification = c(1, 0), legend.position = "bottom", legend.text = element_text(size = 5), plot.title = element_text(hjust = 0.5)) + 
    scale_fill_manual(values = sapply(sort(unique(C$NetworkNames)), function(x) yeo7_color$color[yeo7_color$NetworkNames %in% x]), na.value = "white")
  
  multi.page2 <- ggarrange(yeo7_plot, cluster_plot + rremove("legend"), ncol=2, nrow = 1, common.legend = TRUE, legend = "bottom") 
  
  ggsave(filename=paste0("~/Documents/papers/groupwise mixture/JASA_casestudy/heritability_based_clustering/cluster_", i, ".pdf"), multi.page2,width = 10, height=5)
}



#### compute the overlap percentages compared to the Yeo7 Networks

dk_to_yeo7_key = read.csv("~/Documents/papers/groupwise mixture/JASA_casestudy/real_data_analysis/Yeo_DK_key.csv", row.names = 1)
dk_to_yeo7_key <- dk_to_yeo7_key[, -1]
region_num <- read_excel("~/Documents/papers/groupwise mixture/JASA_casestudy/real_data_analysis/brainregion_number.xlsx")
colnames(region_num) <- c("num_label", "roi")
tmp1 <- merge(region_num, dk_to_yeo7_key, by = "roi")

yeo7 = read.csv("~/Documents/papers/groupwise mixture/JASA_casestudy/real_data_analysis/yeo7_NetworkNames.csv"); 
colnames(yeo7) = c("yeo_krienen", "NetworkNames")
yeo7_network_names <- yeo7$NetworkNames;  yeo7_network_names[4] <- "Salience"
dk_to_yeo7_key = merge(tmp1, yeo7, by = "yeo_krienen")
rm(yeo7, tmp1, region_num)

yeo7composition <- list()
for (k in 1:7){
  indd <- which(dk_to_yeo7_key$yeo_krienen == k);  yeo7composition[[k]] <- dk_to_yeo7_key$num_label[indd]
}


load("~/Documents/papers/groupwise mixture/JASA_casestudy/heritability_based_clustering/global_local_heritability.RData")

OLP <- matrix(0, 4, 7)
for (i in 1:4){
  P <- dat[dat$cluster ==i, ] 
  nt <- dim(P)[1]
  for (k in 1:7){ 
    setp <- intersect(P$num_label, yeo7composition[[k]]) 
    OLP[i,k] <- length(setp)/nt
  }
}

rm(P, nt, setp)

colnames(OLP) <- yeo7_network_names


















  
  
