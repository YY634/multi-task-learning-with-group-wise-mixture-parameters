rm(list = ls())
library(readxl)
library(dplyr)
library(RColorBrewer)
library(ggseg)
library(ggpubr)

dk_to_yeo7_key = read.csv(paste0(path, "Yeo_DK_key.csv"), row.names = 1)
yeo7 = read.csv(paste0(path, "yeo7_NetworkNames.csv"))
colnames(yeo7) = c("yeo_krienen", "NetworkNames")
dk_to_yeo7_key = merge(dk_to_yeo7_key, yeo7, by = "yeo_krienen")
rm(yeo7)

colnames(dk_to_yeo7_key)[3] = "label"
######## Visualze yeo7 Networks
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


load(paste0(path,"PLOT.RData"))
PLOT<-as.data.frame(PLOT)


load(paste0(path,"/dkt_gene_expr_modified.RData"))
colnames(dkt_gene)[5:72] = brain_names
rm(dk)
#yeo7_color = data.frame(color = brewer.pal(7, "Set3"), 
#                        NetworkNames = levels(dk_to_yeo7_key$NetworkNames),
#                        stringsAsFactors = F)

yeo7_color = data.frame(color = brewer.pal(7, "Set3"), 
                        NetworkNames = names(table(dk_to_yeo7_key$NetworkNames)),
                        stringsAsFactors = F)

select_cluster<- PLOT[,]
select_cluster <- select_cluster[c(1,8,15,22,29,36),]

# within_cluster_gene_expr = list()
# for(i in 1:nrow(select_cluster)){
#   g = sapply(strsplit(select_cluster[i, "gene"], "[,]")[[1]], function(x) gene_symbols[gene_symbols[, 1] %in% x, 2])
#   rois = strsplit(as.character(select_cluster[i, "rois_cluster"]), "[,]")[[1]]
#   clusterd_ind = rep(0, 68)
#   clusterd_ind[sapply(rois, function(x) which(brain_names %in% x))] = 1
#   tmp = dkt_gene[dkt_gene$gene_symbol %in% g, -c(1:4)]
#   tmp = as.data.frame(t(tmp))
#   colnames(tmp) = dkt_gene$gene_symbol[dkt_gene$gene_symbol %in% g]
#   tmp$clustered = clusterd_ind
#   within_cluster_gene_expr[[i]] = tmp
# }

within_cluster_gene_expr = list()
for(i in 1:nrow(select_cluster)){
  g = sapply(strsplit(select_cluster[i, "gene_label"], "[,]")[[1]], function(x) dkt_gene[dkt_gene[, 2] %in% x, 2])
  rois = strsplit(as.character(select_cluster[i, "rois_cluster"]), "[,]")[[1]]
  clusterd_ind = rep(0, 68)
  clusterd_ind[sapply(rois, function(x) which(brain_names %in% x))] = 1
  tmp = dkt_gene[dkt_gene$gene_symbol %in% g, -c(1:4)]
  tmp = as.data.frame(t(tmp))
  colnames(tmp) = dkt_gene$gene_symbol[dkt_gene$gene_symbol %in% g]
  tmp$clustered = clusterd_ind
  within_cluster_gene_expr[[i]] = tmp
}


for(i in 1:length(within_cluster_gene_expr)){
  results = within_cluster_gene_expr[[i]] %>%
    mutate(label = rownames(within_cluster_gene_expr[[i]])) %>%
    mutate(label_id = sapply(rownames(within_cluster_gene_expr[[i]]) , function(x) strsplit(x, "[_]")[[1]][2]))
  d = results %>% group_by(clustered) %>% 
    summarise(across(!c(label, label_id),mean))
  d = round(d, 4)
  p = unlist(results %>% 
               summarise(across(!c(clustered, label, label_id), ~ t.test(., clustered)$p.value)))
  p = formatC(p, format = "e", digits = 4)
  d = rbind(d, c("p", p))
  
  p_clusters = list()
  for(j in 1:(ncol(results) - 3)){
    p_clusters[[j]] <- results %>% 
      mutate(clustered = factor(clustered)) %>%
      ggplot(aes(x = clustered, y = results[, j], group = clustered, fill = clustered)) +
      geom_boxplot(alpha = 0.4) +
      labs(x = "Cluster or Not", title = paste0("Gene:", colnames(results)[j]), y = "Normalized Gene Expression")  +
      theme(legend.justification = c(1, 0),
            legend.position = "none",
            legend.text = element_text(size = 5),
            plot.title = element_text(hjust = 0.5),
            axis.title.x=element_blank(),
            axis.text.x=element_text(size = 10))+
      scale_x_discrete(labels=c("0" = "Outside subnetwork", "1" = "Within subnetwork"))+
      scale_fill_manual(values = c(  "#D9D9D9","#E41A1C"))
  }
  
  d_old <- d
  d <- d[,-1 ]
  #scale_fill_manual(values = c("#D9D9D9", "#E41A1C"))
  stable.p <- ggtexttable(d, rows = c("Outside Cluster", "Within Cluster", "p-value"))
  p_clusters[[length(p_clusters) + 1]] = stable.p
  multi.page <- ggarrange(plotlist = p_clusters, nrow = 1) # for one plot per page
  ggsave(filename=paste0("/Volumes/Students/yh567/Manuscripts/Project_for_YY/plot/select_row_", i, "_Gene_boxplot_v2.pdf"), multi.page,width=8,height=3)
  
  # multi.page <- ggarrange(stable.p)
  # ggexport(multi.page,
  #          width = 16, height = 10,
  #          filename=paste0("Y:/wd278/Projects with Others/YY/Results/Visualization/select_cluster/select_row_", i, "_Gene_boxplot_v2.pdf"))

   
  tmp = results %>% 
    left_join(dk_to_yeo7_key, by = "label") %>%
    filter(clustered == 1) %>%
    left_join(yeo7_color, by = "NetworkNames") %>%
    mutate(color = as.character(color))
  cluster_plot = tmp%>%
    ggseg(mapping=aes(fill= NetworkNames),atlas = dk, 
          position="stacked", colour="black") +
    labs(title = "Identified brain subnetwork") +
    theme(legend.justification = c(1, 0),
          legend.position = "bottom",
          legend.text = element_text(size = 5),
          plot.title = element_text(hjust = 0.5)) + 
  scale_fill_manual(values = sapply(sort(unique(tmp$NetworkNames)), function(x) yeo7_color$color[yeo7_color$NetworkNames %in% x]),
                    na.value = "white")
  
  multi.page2 <- ggarrange(yeo7_plot, 
                           cluster_plot + rremove("legend"), ncol=2, nrow = 1,
                           common.legend = TRUE, legend = "bottom") # for one plot per page
  #ggexport(multi.page2, 
  #         filename=paste0("/Volumes/Students/yh567/Manuscripts/Project_for_YY/plot/select_row_", i, "_cluster_boxplot.pdf"))
  ggsave(filename=paste0("/Volumes/Students/yh567/Manuscripts/Project_for_YY/plot/select_row_", i, "_cluster_boxplot.pdf"), 
         multi.page2,width = 10, height=5)
  
  blank_plot <- ggplot() + theme_void()
  
  multi.page3 <- ggarrange(blank_plot,multi.page2,
                           blank_plot,blank_plot,
                           blank_plot,multi.page,
                           ncol = 2,
                           nrow = 3,
                           labels=c("(a)","","(b)",""),
                           heights = c(5,0.35,3.5),
                           widths = c(0.5,10))
  
  ggsave(filename=paste0("/Volumes/Students/yh567/Manuscripts/Project_for_YY/plot/select_row_", i, "_cluster_gene_boxplot.pdf"), 
         multi.page3,width = 10, height=8)
}

