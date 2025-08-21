
library(plyr)
library(MASS)
library(Matrix)
library(plus)
library(mclust)
library(grplasso) #library(doParallel)

start.time<-Sys.time()

options(warn = 1)
load("/mnt/project/Data/X_Rdata/chr20_X.RData")
Y<-read.csv("/mnt/project/Data/Final_phenotype_scaled.csv")
U<-read.csv("/mnt/project/Data/Final_covariate_scaled.csv")
rm(n_X, n_X_group_sizes) 

# In UKB, different chrs may have slightly different row. Make sure X, Y, U have common rows.
U<-U[match(rownames(n_s_PC_X),U$ID),]
Y<-Y[match(rownames(n_s_PC_X),Y$ID),]
sum(rownames(n_s_PC_X)!=U$ID)
sum(rownames(n_s_PC_X)!=Y$ID)
U<-U[,-1]
U<-as.matrix(U)
Y<-Y[,-1]
Y<-as.matrix(Y)
Y<-Y/100

X <- n_s_PC_X; design <- cbind(X, U); group_sizes <- n_PC_group_sizes; grouping <- n_PC_grouping;  low <- n_pc_low;  upp <- n_pc_upp
p <- dim(Y)[2]; n <- dim(X)[1]; d <- dim(X)[2]; q0 <- q; q <- n_q; m <- dim(U)[2] #m is the # of demographic factors 
#q0 is the original number of groups before combining strongly correlated groups. q0=number of genes on chr1.
rm(n_s_PC_X, n_pc_d, n_PC_group_sizes, n_PC_grouping, n_pc_low, n_pc_upp, n_q)
grouping_all <- c(grouping, c((q+1):(q+m)))

lam_0 <- opt_lam <- 6800;  lam_thred <- 5800

#the records of iterations
Betas <- list()  #Betas[[t]] <- matrix(0, n_pc_d, p)
Blabels <- list()  #Blabels[[t]] <- matrix(rep(1:p, n_q), ncol=p, byrow = T)
means <- list()  #means[[t]] <- matrix(0, n_pc_d, p)
classifications <- list()  #classifications[[t]] <- matrix(0, n_q, p)
cluster_nums <- list()  #cluster_nums[[t]] <- rep(0, n_q)
covariances <- list()  #covariances[[t]] is also a list.


#to obtain the starting point
Beta_0 <- matrix(0, d, p)
Beta_0_labels <- matrix(rep(1:p, q), ncol=p, byrow = T)
A_hat <- matrix(0, m, p) # A_hat is the starting point [alpha_1, ..., alpha_p]
for (j in 1:p){
  print(j)
  fit <- grplasso(design, y=Y[,j], index=grouping_all, lambda=lam_0, model=LinReg(), penscale=sqrt, control=grpl.control(update.hess="lambda", trace=0), center=FALSE)
  est <- fit$coefficients[,1]; bj_est <- est[1:d]; alpha_est <- est[(d+1):(d+m)]
  
  if (all(est==0)){}else {
    indicator <- c(1:q)
    for (k in 1:q){
      if (all(bj_est[low[k]:upp[k]]==0)){indicator[k] <- 0}
    }  
    indicator <- indicator[which(indicator!=0)]
    bj_list <- list()
    for (l in indicator){
      bj_list[[l]] <- c(low[l]:upp[l])
    } 
    bj_ind <- unlist(bj_list)
    alpha_ind <- NULL; for (h in 1:m){ if (alpha_est[h]!=0){ alpha_ind <- c(alpha_ind, h) } }
    
    if (length(bj_ind) ==0){ 
      A_hat[,j] <- alpha_est
    }else if (length(bj_ind) == 1){
      Xtp <- cbind(X[,bj_ind], U[,alpha_ind]) 
      tp <- solve( t(Xtp) %*% Xtp ) %*% t(Xtp) %*% Y[,j]
      Beta_0[bj_ind,j] <- tp[1:length(bj_ind)]
      if (length(alpha_ind) >0){ A_hat[alpha_ind,j] <- tp[(length(bj_ind)+1):(length(bj_ind)+length(alpha_ind))] }
    }else { 
      mcpX <- cbind(X[,bj_ind], U[,alpha_ind]) 
      mcp <- plus(mcpX, Y[,j], method = "mc", normalize = TRUE, intercept = FALSE) 
      tp <- mcp$beta[dim(mcp$beta)[1],]
      Beta_0[bj_ind,j] <- tp[1:length(bj_ind)]
      if (length(alpha_ind) >0){ A_hat[alpha_ind,j] <- tp[(length(bj_ind)+1):(length(bj_ind)+length(alpha_ind))] }
    }
  }
  for (i in 1:q){
    if (all(Beta_0[low[i]:upp[i],j]==0)){Beta_0_labels[i,j] <- 0}
  }
}
Betas[[1]] <- Beta_0
Blabels[[1]] <- Beta_0_labels #Blabels only gives information which groups are not zero not classification labels
SIG_hat <-  t(Y- X %*% Beta_0 - U %*% A_hat) %*% (Y- X %*% Beta_0 - U %*% A_hat)/n
siginv <- solve(SIG_hat)
alpha_hat <- solve(t(U) %*% U) %*% t(U) %*% (Y- X%*% Beta_0) %*% rowSums(siginv)/sum(siginv)
rm(i,j,k,l, fit, mcp, est, bj_list, bj_ind, bj_est, alpha_est, alpha_ind, indicator, mcpX, tp, Beta_0, Beta_0_labels, A_hat)


Beta_temp <- Betas[[1]]
Blabels_temp <- Blabels[[1]]  #Blabels_temp only gives information which groups are not zero but not classification labels
means_temp <- matrix(0, d, p)
covs_temp <- list() 
classifications_temp <- matrix(0, q, p) #this recrods the classification (including zeros)
cluster_nums_temp <- rep(0, q)
for (i in 1:q){
  if (length(which(Blabels_temp[i,]!=0)) <= 7){ 
    classifications_temp[i,] <- rep(0, p)  
    means_temp[low[i]:upp[i],] <- 0
    covs_temp[[i]] <- 0
    cluster_nums_temp[i] <- 0
  } else if ( length(which(Beta_temp[(low[i]:upp[i]),Blabels_temp[i,]] !=0))/length(Beta_temp[(low[i]:upp[i]),Blabels_temp[i,]]) < 1/2){ 
    classifications_temp[i,] <- rep(0, p)  
    means_temp[low[i]:upp[i],] <- 0
    covs_temp[[i]] <- 0
    cluster_nums_temp[i] <- 0
  } else{
    temp <- Beta_temp[(low[i]:upp[i]),Blabels_temp[i,]]
    mod <- Mclust(t(temp), G=1:6, control=emControl(eps=0.001)) # maybe add the shape of the clusters?
    index1 <- which(Blabels_temp[i,] ==0)
    if (length(index1) >0){
      part1 <- as.data.frame(cbind(index1, 0))
      colnames(part1) <- c("index", "label")
      index2 <- which(Blabels_temp[i,] !=0)
      part2 <- as.data.frame(cbind(index2, mod$classification))
      colnames(part2) <- c("index", "label")
      est <- rbind(part1, part2)
      est <- est[order(est$index),]
      classifications_temp[i,] <- est$label
    }else {
      classifications_temp[i,] <- mod$classification
    }
    clasind <- which(classifications_temp[i,] != 0)
    for (j in clasind){
      if (group_sizes[i]==1){
        means_temp[low[i]:upp[i],j] <- mod$parameters$mean[classifications_temp[i,j]]
      }else {
        means_temp[low[i]:upp[i],j] <- mod$parameters$mean[,classifications_temp[i,j]]
      }
    }
    covs_temp[[i]] <- mod$parameters$variance$sigma
    cluster_nums_temp[i] <- mod$G
  }
}
means[[1]] <- means_temp
covariances[[1]] <- covs_temp
classifications[[1]] <- classifications_temp
cluster_nums[[1]] <- cluster_nums_temp #best_models[[1]] <- best_models_temp
rm(i, temp, mod, index1, part1, index2, part2, est, Beta_temp, Blabels_temp) 

t <- 1
repeat{
  #update Beta.  
  opt_lam <- max(opt_lam*0.85, lam_thred)
  Beta_temp <- matrix(0, d, p)
  Blabels_temp <- matrix(rep(1:p, q), ncol=p, byrow = T)
  for (j in 1:p){
    #update the on-support part
    n_g_ind <- which(classifications_temp[,j]!=0)
    z_g_ind <- which(classifications_temp[,j]==0) 
    tildeyj <- Y[,j] - U %*% alpha_hat
    if (length(n_g_ind) >0){
      nonzero_list <- list()
      cov_list <- list()
      for (i in n_g_ind){ 
        nonzero_list[[i]] <- c(low[i]:upp[i])
        covs <- covs_temp[[i]]
        class_ind <- classifications_temp[i,j]  
        if (group_sizes[i] == 1){
          if (length(covs)==1){cov_list[[i]] <- 1/covs} else{cov_list[[i]] <- 1/covs[class_ind]}
        } else{cov_list[[i]] <- ginv(covs[,,class_ind])}
      }
      nonzero_ind <- unlist(nonzero_list)
      Xn <- X[, nonzero_ind]
      muj <- means_temp[nonzero_ind,j]
      cov_list[sapply(cov_list, is.null)] <- NULL;  covj <- as.matrix(bdiag(cov_list))
      
      leftmatrix <- (t(Xn) %*% Xn)/n + covj/3 
      rightmatrix <- t(Xn) %*% tildeyj/n + covj%*%muj/3
      Beta_temp[nonzero_ind, j] <- solve(leftmatrix) %*% rightmatrix
      
      zero_list <- list()
      for (k in z_g_ind){  zero_list[[k]] <- c(low[k]:upp[k])  }
      zero_ind <- unlist(zero_list)
      Xz <- X[, zero_ind];   rm(i,k)
      
    }else{ Xz <- X; zero_ind <- c(1:d) }
    
    #update the off-support part
    indicator <- c()
    res <- as.vector(tildeyj - X %*% Beta_temp[, j])
    fit1 <- grplasso(Xz, y=res, index=grouping[zero_ind], lambda=opt_lam, model=LinReg(), penscale=sqrt, control=grpl.control(update.hess="lambda", trace=0), center=FALSE)
    grp_est1 <- fit1$coefficients[,1]
    if (all(grp_est1==0)){}else {
      indicator <- z_g_ind
      count <- 0
      for (k in z_g_ind){
        l <- count + 1
        u <- count + group_sizes[k]
        count <- count + group_sizes[k]
        if (all(grp_est1[l:u]==0)){indicator[which(indicator==k)] <- 0}
      }  
      indicator <- indicator[which(indicator!=0)] #newly entered groups, that is previously zero groups but now nonzero
    }
    
    #update the support
    new_n_g_ind <- sort(c(n_g_ind, indicator))
    if (length(new_n_g_ind) ==0){}else if(length(new_n_g_ind) ==1){
      new_nonzero_ind <- c(low[new_n_g_ind]:upp[new_n_g_ind]);  xtp <- X[, new_nonzero_ind] 
      Beta_temp[new_nonzero_ind,j] <- solve( t(xtp) %*% xtp ) %*% t(xtp) %*% tildeyj
    }else{
      new_nonzero_list <- list()
      for (k in new_n_g_ind){
        new_nonzero_list[[k]] <- c(low[k]:upp[k]) 
      }
      new_nonzero_ind <- unlist(new_nonzero_list)
      Xnn <- X[, new_nonzero_ind]   
      fit2 <- grplasso(Xnn, y=tildeyj, index=grouping[new_nonzero_ind], lambda=opt_lam, model=LinReg(), penscale=sqrt, control=grpl.control(update.hess="lambda", trace=0), center=FALSE)
      grp_est2 <- fit2$coefficients[,1]
      indicatornew <- new_n_g_ind
      count <- 0
      for (k in new_n_g_ind){
        l <- count + 1
        u <- count + group_sizes[k]
        count <- count + group_sizes[k]
        if (all(grp_est2[l:u]==0)){indicatornew[which(indicatornew==k)] <- 0}
      }  
      indicatornew <- indicatornew[which(indicatornew !=0)] 
      Beta_temp[, j] <- rep(0, d)
      grp_list <- list()
      for (i in indicatornew){
        grp_list[[i]] <- c(low[i]:upp[i])
      }
      grp_ind <- unlist(grp_list)
      if (length(grp_ind)==0){}else if (length(grp_ind) ==1){
        xtp <- X[,grp_ind];  Beta_temp[grp_ind, j] <- solve( t(xtp) %*% xtp ) %*% t(xtp) %*% tildeyj
      }else {
        mcp <- plus(X[,grp_ind], as.vector(tildeyj), method = "mc", normalize = TRUE, intercept = FALSE)
        Beta_temp[grp_ind, j] <- mcp$beta[dim(mcp$beta)[1],]
      }
    }
    
    for (i in 1:q){ if (all(Beta_temp[low[i]:upp[i],j]==0)){Blabels_temp[i,j] <- 0} }
  }
  
  t <- t+1
  Betas[[t]] <- Beta_temp;  Blabels[[t]] <- Blabels_temp  #best_lambda[[t]] <- mean(best_lams)
  rm(n_g_ind, nonzero_list, cov_list, covs, 
     class_ind, nonzero_ind, Xn, muj, covj, z_g_ind, zero_list, zero_ind, Xz, 
     leftmatrix, rightmatrix, res, fit1, fit2, mcp, grp_list, grp_ind, indicator, indicatornew, i, j,
     count, grp_est1, grp_est2)
  
  #update SIG_hat and alpha_hat
  SIG_hat <-  t(Y- X %*% Beta_temp - U %*% (alpha_hat %*% t(rep(1, p)))) %*% (Y- X %*% Beta_temp - U %*% (alpha_hat %*% t(rep(1, p))))/n
  siginv <- solve(SIG_hat)
  alpha_hat <- solve(t(U) %*% U) %*% t(U) %*% (Y- X%*% Beta_temp) %*% rowSums(siginv)/sum(siginv)
  
  #clustering on Beta_temp
  means_temp <- matrix(0, d, p)
  covs_temp <- list() #best_models_temp <- c()
  classifications_temp <- matrix(0, q, p)
  cluster_nums_temp <- rep(0, q)
  for (i in 1:q){
    if (length(which(Blabels_temp[i,]!=0)) <= 7){ 
      classifications_temp[i,] <- rep(0, p)  
      means_temp[low[i]:upp[i],] <- 0
      covs_temp[[i]] <- 0
      cluster_nums_temp[i] <- 0
    } else if ( length(which(Beta_temp[(low[i]:upp[i]),Blabels_temp[i,]] !=0))/length(Beta_temp[(low[i]:upp[i]),Blabels_temp[i,]]) < 1/2){ 
      classifications_temp[i,] <- rep(0, p)  
      means_temp[low[i]:upp[i],] <- 0
      covs_temp[[i]] <- 0
      cluster_nums_temp[i] <- 0
    } else{
      temp <- Beta_temp[(low[i]:upp[i]),Blabels_temp[i,]]
      mod <- Mclust(t(temp), G=1:6, control=emControl(eps=0.001)) # maybe add the shape of the clusters?
      index1 <- which(Blabels_temp[i,] ==0)
      if (length(index1) >0){
        part1 <- as.data.frame(cbind(index1, 0))
        colnames(part1) <- c("index", "label")
        index2 <- which(Blabels_temp[i,] !=0)
        part2 <- as.data.frame(cbind(index2, mod$classification))
        colnames(part2) <- c("index", "label")
        est <- rbind(part1, part2)
        est <- est[order(est$index),]
        classifications_temp[i,] <- est$label
      }else {
        classifications_temp[i,] <- mod$classification
      }
      clasind <- which(classifications_temp[i,] != 0)
      for (j in clasind){
        if (group_sizes[i]==1){
          means_temp[low[i]:upp[i],j] <- mod$parameters$mean[classifications_temp[i,j]]
        }else {
          means_temp[low[i]:upp[i],j] <- mod$parameters$mean[,classifications_temp[i,j]]
        }
      }
      covs_temp[[i]] <- mod$parameters$variance$sigma
      cluster_nums_temp[i] <- mod$G
    }
  }
  means[[t]] <- means_temp
  covariances[[t]] <- covs_temp
  classifications[[t]] <- classifications_temp
  cluster_nums[[t]] <- cluster_nums_temp  #best_models[[t]] <- best_models_temp
  rm(i, temp, mod, index1, part1, index2, part2, est)
  
  if(sum((Beta_temp-Betas[[t-1]])^2)/(d*p) < 0.001 | t>4){break}
} 

print(Sys.time()-start.time)
run_time<-Sys.time()-start.time

save(n, d, p, q, q0, Betas, Blabels, classifications, means, covariances, cluster_nums, run_time,
     file = "chr20_result.RData")

system(paste0("dx upload chr20_result.RData"))

system("dx upload ukb_chr20.R")

system("dx terminate job-J1v05zQJbfG8PbYyyPBZPQpG")

#class <- classifications[[10]]; label <- Blabels[[10]]
#check <- rep(0, q)
#for (i in 1:q){ check[i] <- length(which(class[i,] !=0)) }
#length(which(check !=0)) 
