library(plyr)
library(Matrix)

chr<-18

print(paste0("start: chr",chr))
start.time<-Sys.time()

load(paste0("/mnt/project/Data/X_Rdata/chr",chr,"_SNPs.RData")) 
# "chr21_SNPs.RData" contains genotype, an n*d matrix records the d SNP values of n individuals

# match the individuals in Y (in "brain_volumes.rda") and X
#load("/gpfs/gibbs/project/zhang_heping/wd278/Project_with_others/Yisha/Data/Phenotype/brain_volumes.rda")
d<-read.csv("/mnt/project/Data/Final_phenotype_covariate.csv")
X <- as.matrix(genotype); X[is.na(X)] <- 0; X <- X[rownames(X) %in% d[,1],]; rm(d)

load(paste0("/mnt/project/Results/SNP_group/chr",chr,"_SNP_group.RData"))
SNP_group<-M1; rm(M1) #make the naming consistent
comSNP <- intersect( colnames(genotype), SNP_group$SNP )
X <- X[, which(colnames(X) %in% comSNP)]; SNP_group <- SNP_group[which(SNP_group$SNP %in% comSNP),]
rm(genotype, comSNP)

#check
if(dim(X)[2] != dim(SNP_group)[1]){
  stop("dim(X)[2] != dim(SNP_group)[1]")
}
if(length(colnames(X)) != length(unique(SNP_group$SNP))){
  stop("length(colnames(X)) != length(unique(SNP_group$SNP))")
}#the answers to these two commands should be "TRUE"
if(!all(diff(SNP_group$Pos) > 0)){
  stop("!all(diff(SNP_group$Pos) > 0)")
}#if the positions of the SNPs are in ascending order, the answer should be "TRUE"

# There may exist overlapping SNPs, that is, SNPs that belong to two adjacent genes. The overlapping part is a different group 
# key <- paste(SNP_group$Chr, SNP_group$Pos, sep = ":");  length(unique(key)) == dim(SNP_group)[1]

# to reorder the rows of SNP_group so that all the SNPs with the same gene stay together
groups <- unique(SNP_group$Gene);  q <- length(groups);  M <- NULL
for (i in 1:q){
  Mt <- SNP_group[which(SNP_group$Gene == groups[i]), ];  Mt <- Mt[order(Mt$Pos), ]
  M <- rbind(M, Mt)
}
#all(M$SNP == SNP_group$SNP)
#clean up a little bit, a very small set of SNPs removed
freq_M <- as.data.frame(table(M$Gene))
colnames(freq_M) <- c("Gene", "Freq")
M_freq <- join(M, freq_M)
checkind <- unique(c(grep(",", M$Gene), grep(";", M$Gene))); checkone <- which(M_freq$Freq == 1); removeind <- intersect(checkind, checkone)
M <- M[-removeind,];    rm(checkind, checkone, removeind);  SNP_group <- M;   rm(M, Mt, M_freq, freq_M)
#checkind <- c(grep("upstream", M$Func), grep("downstream", M$Func)); checkind <- unique(checkind);

# Re-order the columns of X.
X_olp <- X[, SNP_group$SNP]
X <- scale(X_olp, center = TRUE, scale = FALSE)
n <- nrow(X); d <- ncol(X); rm(X_olp)
#duplicated_cols <- duplicated(colnames(X)); check <- which(duplicated_cols)

# get the group sizes, and the starting and ending indices of the groups
freq <- as.data.frame(table(SNP_group$Gene))
colnames(freq) <- c("Gene", "Freq")
groups <- as.data.frame(unique(SNP_group$Gene))#the order of the groups are the same as in the original SNP_group; correct order
colnames(groups) <- c("Gene")  
group_freq <- join(groups, freq); q <- dim(group_freq)[1] #the order of group_freq will be the same as the order in groups 
rm(freq, groups)

X_group_sizes <- group_freq$Freq
low <- rep(0, q)
upp <- rep(0, q)
count <- 0
for (i in 1:q){
  low[i] <- count + 1
  upp[i] <- count+ X_group_sizes[i]
  count <- count + X_group_sizes[i]
}
X_grouping <- rep(1:q, X_group_sizes)
rm(i, count)


#to obtain the PCs of the groups
PC_X <- NULL
PC_group_sizes <- rep(0, q)
for (i in 1:q){
  if(i%%100==0){print(i)}
  pca <- prcomp(X[, low[i]:upp[i]], scale. = TRUE)
  props <- (pca$sdev)^2/sum((pca$sdev)^2)
  cumprops <- cumsum(props)
  cutoff <- which(cumprops > 0.8)[1] #set a low 0.8 to make sure that top PCs of different groups 
  PC <- pca$x[, 1:cutoff]
  PC_group_sizes[i] <- cutoff
  PC_X <- cbind(PC_X, PC)
}
rm(i, pca, props, cumprops, cutoff, PC)

PC_grouping <- rep(1:q, PC_group_sizes)
pc_d <- length(PC_grouping)
pc_low <- rep(0, q)
pc_upp <- rep(0, q)
count <- 0
for (i in 1:q){
  pc_low[i] <- count + 1
  pc_upp[i] <- count+ PC_group_sizes[i]
  count <- count + PC_group_sizes[i]
}
rm(i, count)
s_PC_X <- scale(PC_X, center = TRUE, scale = TRUE) #should used normalized design since the coefficients are largely affected by the scale of the corresponding covariates 
PC_cor <- t(s_PC_X) %*% s_PC_X/(n-1) - diag(pc_d)

##### check whether there are strong correlations among the groups of PCs. 
#max(abs(PC_cor)); length(which(abs(PC_cor) > 0.9))/2
print(paste0("max(abs(PC_cor)):",max(abs(PC_cor))))
print(paste0("length(which(abs(PC_cor) > 0.9))/2:",length(which(abs(PC_cor) > 0.9))/2))

# if max(abs(PC_cor)) <= 0.9 or length(which(abs(PC_cor) > 0.9, arr.ind = TRUE))/4 <= min(0.000001*pc_d*(pc_d-1)/2, 4)
# the group-combining steps below can be omitted. s_PC_X is the needed design matrix
if((max(abs(PC_cor)) <= 0.9)|length(which(abs(PC_cor) > 0.9, arr.ind = TRUE))/4 <= min(0.000001*pc_d*(pc_d-1)/2, 4)){
  save(X, PC_group_sizes, PC_grouping, pc_d, pc_low, pc_upp, s_PC_X, n, q, file = paste0("chr",chr,"_X.RData"))
  system(paste0("dx upload chr",chr,"_X.RData"))
}else{
  # If there exist strong correlations among the groups, execute the code below. It combines strongly correlated groups of SNPs.
  # to find the pairs of groups with strong correlations
  maxcors <- matrix(0, q, q)
  for (i in 1:(q-1)){
    for (l in (i+1):q){
      cor <- PC_cor[pc_low[i]:pc_upp[i], pc_low[l]:pc_upp[l]]
      maxcors[i,l] <- max(max(cor), abs(min(cor)))
    }
  } #maxcors[i,l] is the max correlation between the i-th group of PCs and the l-th group of PCs. It is upper triangle
  
  maxcor_temp <- maxcors
  maxc <- max(maxcor_temp)
  edges <- NULL
  while (maxc >= 0.3){
    maxind <- which(maxcor_temp==maxc, arr.ind=TRUE)
    edges <- rbind(edges, c(maxind[1,],maxc))
    maxcor_temp[maxind[1,1], maxind[1,2]] <- 0
    maxc <- max(maxcor_temp)
  }
  # the rows of edges record the (i,l) and the max correlation between the i-th group of PCs and the l-th group of PCs
  # and the pairwise group correlations are in descending order
  colnames(edges) <- c("row", "col", "max_corr")
  rm(i,l, maxcor_temp, maxc, maxind, cor)
  
  cliques <- list()
  cliques[[1]] <- edges[1, 1:2]
  m <- 1 #the number of cliques
  U <- cliques[[1]]
  
  cut7 <- max(which(edges[,3] > 0.7))
  for (i in 2:cut7){
    mem <- 0; mem1 <- 0; mem2 <- 0
    if (length(intersect(edges[i,1:2], U)) ==2){
      log1 <- sapply(seq_along(cliques), function(k) length(intersect(cliques[[k]], edges[i,1])) !=0)
      mem1 <- which(log1 == TRUE)
      log2 <- sapply(seq_along(cliques), function(k) length(intersect(cliques[[k]], edges[i,2])) !=0)
      mem2 <- which(log2 == TRUE)
      if (mem1 != mem2){
        cliques[[min(mem1, mem2)]]<- union(cliques[[mem1]], cliques[[mem2]])
        cliques[[max(mem1, mem2)]]<- NULL
        m <- m-1
      }
    } else if (length(intersect(edges[i,1:2], U)) ==1){
      logg <- sapply(seq_along(cliques), function(k) length(intersect(cliques[[k]], edges[i,1:2])) !=0)
      mem <- which(logg == TRUE)
      cliques[[mem]] <- union(cliques[[mem]], edges[i,1:2])
      U <- union(U, edges[i,1:2])
    } else {
      m <- m+1
      cliques[[m]] <- edges[i,1:2]
      U <- union(U, edges[i,1:2])
    }
  }
  rm(logg, mem, log1, mem1, log2, mem2)
  
  #to record the old-new grouping correspondence
  o_cliques <- sapply(seq_along(cliques), function(k) unname(sort(cliques[[k]])))
  first_ele <- rep(0, length(o_cliques))
  first_ele <- sapply(seq_along(o_cliques), function(k) first_ele[k] <- o_cliques[[k]][1])
  o_first_ele <- sort(first_ele, index.return=T)
  oo_cliques <- unname(o_cliques[o_first_ele$ix]) #order the cliques by making the first element of each clique in ascending
  rm(o_cliques, first_ele, o_first_ele)
  reorganized <- unlist(oo_cliques) 
  remain <- c(1:q)[is.na(pmatch(c(1:q),reorganized))]
  n_q <- length(remain) + length(oo_cliques) 
  #n_q is the new total number of groups
  
  n_X_group_sizes <- rep(0, n_q)
  old_new_groups <- matrix(0,length(remain),2) # first column has the old group indices, second column has the new group indices
  reordering_ind <- rep(0, d)
  n_low <- rep(0, n_q)
  n_upp <- rep(0, n_q)
  count <- 0
  for (l in 1:length(remain)){
    n_X_group_sizes[l] <- X_group_sizes[remain[l]]
    old_new_groups[l,] <- c(remain[l], l)
    n_low[l] <- count + 1
    n_upp[l] <- count + X_group_sizes[remain[l]]
    count <- n_upp[l]
    reordering_ind[n_low[l]:n_upp[l]] <- c(low[remain[l]]:upp[remain[l]]) 
  }
  reordering_ind <- reordering_ind[-which(reordering_ind==0)]
  for (i in 1:length(oo_cliques)){
    n_X_group_sizes[i+length(remain)] <- sum(X_group_sizes[oo_cliques[[i]]])
    block <- cbind(oo_cliques[[i]], rep(i+length(remain), length(oo_cliques[[i]])))
    old_new_groups <- rbind(old_new_groups, block)
    n_low[i+length(remain)] <- count + 1
    n_upp[i+length(remain)] <- count + n_X_group_sizes[i+length(remain)]
    count <- count + n_X_group_sizes[i+length(remain)]
    v <- oo_cliques[[i]]
    for (k in 1:length(v)){
      reordering_ind <- c(reordering_ind, low[v[k]]:upp[v[k]])
    }
  }
  rm(i, l, k, block, v, remain, reorganized, count, U)
  colnames(old_new_groups) <- c("old_group_label", "new_group_label") # "old_group_label" contains 1:q, while "new_group_label" have repeats
  old_new_groups <- as.data.frame(old_new_groups)
  df <- old_new_groups[order(old_new_groups$old_group_label), ]
  old_new_groups <- cbind(group_freq$Gene, df); colnames(old_new_groups)[1] <- "gene"
  
  
  # n_X is the new X with some groups combined and groups' orders switched. dim(n_X)=dim(X)
  # re-do PCA among the groups of n_X.
  n_X <- X[, reordering_ind]
  n_X_grouping <- rep(1:n_q, n_X_group_sizes)
  #to obtain the new PCs
  n_PC_X <- NULL
  n_PC_group_sizes <- rep(0, n_q)
  for (i in 1:n_q){
    if(i%%100==0){print(i)}
    pca <- prcomp(n_X[, n_low[i]:n_upp[i]], scale. = FALSE)
    props <- (pca$sdev)^2/sum((pca$sdev)^2)
    cumprops <- cumsum(props)
    cutoff <- which(cumprops > 0.8)[1] #set a low 0.8 to make sure that top PCs of different groups 
    PC <- pca$x[, 1:cutoff]
    n_PC_group_sizes[i] <- cutoff
    n_PC_X <- cbind(n_PC_X, PC)
  }
  rm(i, pca, props, cumprops, cutoff, PC)
  
  n_PC_grouping <- rep(1:n_q, n_PC_group_sizes)
  n_pc_d <- length(n_PC_grouping)
  n_pc_low <- rep(0, n_q)
  n_pc_upp <- rep(0, n_q)
  count <- 0
  for (i in 1:n_q){
    n_pc_low[i] <- count + 1
    n_pc_upp[i] <- count+ n_PC_group_sizes[i]
    count <- count + n_PC_group_sizes[i]
  }
  rm(i, count, reordering_ind)
  n_s_PC_X <- scale(n_PC_X, center = TRUE, scale = TRUE) 
  n_PC_cor <- t(n_s_PC_X) %*% n_s_PC_X/(n-1) - diag(n_pc_d)
  
  ##### check the correlations among the groups
  #length(which(n_PC_cor > 0.9, arr.ind = TRUE))/4; length(which(n_PC_cor > 0.8, arr.ind = TRUE))/4
  #which(n_PC_cor > 0.9, arr.ind = TRUE)
  # n_s_PC_X is what we need, the centered and normalized design matrix.
  
  save(n_X, n_X_group_sizes, old_new_groups, n_PC_group_sizes, 
       n_PC_grouping, n_pc_d, n_pc_low, n_pc_upp, n_s_PC_X, n, q, n_q,
       file = paste0("chr",chr,"_X.RData"))
  system(paste0("dx upload chr",chr,"_X.RData"))
}

print(paste0("finished: chr",chr))
print(Sys.time()-start.time)

system("dx upload process_X.R")
