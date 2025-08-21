
library(plyr)
library(MASS)
library(Matrix)
library(plus)
library(mclust)
library(grplasso)

load("path_to/chr18_X.RData"); rm(X, pc_d) 
X <- s_PC_X; rm(s_PC_X)
load("path_to/UKB_U.RData") # these are environmental factors
n <- dim(X)[1];  d <- dim(X)[2];  p <- 68;  m <- dim(U)[2]
opt_lam <- 3200; lam_0 <- 4000;  lam_thred <- 1200 

group_sizes <- PC_group_sizes;  ind_low <- pc_low;  ind_upp <- pc_upp
grouping <- rep(1:q, group_sizes);  grouping_all <- c(grouping, q+c(1:m))
rm(PC_group_sizes, PC_grouping, pc_low, pc_upp)

load("path_to/result_18.RData")
L <- length(Betas)
Beta <- Betas[[L]]; Beta_cluster_num <- cluster_nums[[L]]; Beta_labels <- classifications[[L]]; true_means <- means[[L]]
alpha <- c(0.1, 0.2, 0, rep(0, 10)) 
rm(Betas, Blabels, classifications, cluster_nums, covariances, means, L)

E <- mvrnorm(n, mu=rep(0, p), Sigma=diag(p))
Y <- X %*% Beta + U %*% (alpha %*% t(rep(1, p))) + E
design <- cbind(X, U)

Betas <- list()  #Betas[[t]] <- matrix(0, n_pc_d, p)
Blabels <- list()  #Blabels[[t]] <- matrix(rep(1:p, n_q), ncol=p, byrow = T)
means <- list()  #means[[t]] <- matrix(0, n_pc_d, p)
classifications <- list()  #classifications[[t]] <- matrix(0, n_q, p)
cluster_nums <- list()  #cluster_nums[[t]] <- rep(0, n_q)
#best_models <- list()  #best_models[[t]] <- c()
covariances <- list()  #covariances[[t]] is also a list.

Beta_0 <- matrix(0, d, p)
Beta_0_labels <- matrix(rep(1:p, q), ncol=p, byrow = T)
A_hat <- matrix(0, m, p) # A_hat is the starting point [alpha_1, ..., alpha_p]
for (j in 1:p){
  fit <- grplasso(design, y=Y[,j], index=grouping_all, lambda=lam_0, model=LinReg(), penscale=sqrt, control=grpl.control(update.hess="lambda", trace=0), center=FALSE)
  est <- fit$coefficients[,1]; bj_est <- est[1:d]; alpha_est <- est[(d+1):(d+m)]
  
  if (all(est==0)){}else {
    indicator <- c(1:q)
    for (k in 1:q){
      if (all(bj_est[ind_low[k]:ind_upp[k]]==0)){indicator[k] <- 0}
    }  
    indicator <- indicator[which(indicator!=0)]
    bj_list <- list()
    for (l in indicator){
      bj_list[[l]] <- c(ind_low[l]:ind_upp[l])
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
    if (all(Beta_0[ind_low[i]:ind_upp[i],j]==0)){Beta_0_labels[i,j] <- 0}
  }
}
Betas[[1]] <- Beta_0
Blabels[[1]] <- Beta_0_labels #Blabels only gives information which groups are not zero not classification labels
SIG_hat <-  t(Y- X %*% Beta_0 - U %*% A_hat) %*% (Y- X %*% Beta_0 - U %*% A_hat)/n
siginv <- solve(SIG_hat)
alpha_hat <- solve(t(U) %*% U) %*% t(U) %*% (Y- X%*% Beta_0) %*% rowSums(siginv)/sum(siginv)
rm(i,j,k,l, fit, mcp, est, bj_list, bj_ind, bj_est, alpha_est, alpha_ind, indicator, Beta_0, Beta_0_labels, A_hat, mcpX, tp)


Beta_temp <- Betas[[1]]
Blabels_temp <- Blabels[[1]]  #Blabels_temp only gives information which groups are not zero
means_temp <- matrix(0, d, p)
covs_temp <- list() #best_models_temp <- c()
classifications_temp <- matrix(0, q, p) #this recrods the classification (including zeros)
cluster_nums_temp <- rep(0, q)
for (i in 1:q){
  if (length(which(Blabels_temp[i,]!=0)) <= 5){ 
    classifications_temp[i,] <- rep(0, p)  
    means_temp[ind_low[i]:ind_upp[i],] <- 0
    covs_temp[[i]] <- 0
    cluster_nums_temp[i] <- 0
  } else if ( length(which(Beta_temp[(ind_low[i]:ind_upp[i]),Blabels_temp[i,]] !=0))/length(Beta_temp[(ind_low[i]:ind_upp[i]),Blabels_temp[i,]]) < 1/3){ 
    classifications_temp[i,] <- rep(0, p)  
    means_temp[ind_low[i]:ind_upp[i],] <- 0
    covs_temp[[i]] <- 0
    cluster_nums_temp[i] <- 0
  } else{
    temp <- Beta_temp[(ind_low[i]:ind_upp[i]),Blabels_temp[i,]]
    mod <- Mclust(t(temp), G=1:6, control=emControl(eps=0.001)) # maybe add the shape of the clusters?
    index1 <- which(Blabels_temp[i,] ==0)
    part1 <- as.data.frame(cbind(index1, 0))
    colnames(part1) <- c("index", "label")
    index2 <- which(Blabels_temp[i,] !=0)
    part2 <- as.data.frame(cbind(index2, mod$classification))
    colnames(part2) <- c("index", "label")
    est <- rbind(part1, part2)
    est <- est[order(est$index),]
    classifications_temp[i,] <- est$label
    clasind <- which(classifications_temp[i,] != 0)
    for (j in clasind){
      if (group_sizes[i]==1){
        means_temp[ind_low[i]:ind_upp[i],j] <- mod$parameters$mean[classifications_temp[i,j]]
      }else {
        means_temp[ind_low[i]:ind_upp[i],j] <- mod$parameters$mean[,classifications_temp[i,j]]
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
  opt_lam <- max(opt_lam*0.9, lam_thred) 
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
        nonzero_list[[i]] <- c(ind_low[i]:ind_upp[i])
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
      for (k in z_g_ind){  zero_list[[k]] <- c(ind_low[k]:ind_upp[k])  }
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
      new_nonzero_ind <- c(ind_low[new_n_g_ind]:ind_upp[new_n_g_ind]);  xtp <- X[, new_nonzero_ind] 
      Beta_temp[new_nonzero_ind,j] <- solve( t(xtp) %*% xtp ) %*% t(xtp) %*% tildeyj
    }else{
      new_nonzero_list <- list()
      for (k in new_n_g_ind){
        new_nonzero_list[[k]] <- c(ind_low[k]:ind_upp[k]) 
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
        grp_list[[i]] <- c(ind_low[i]:ind_upp[i])
      }
      grp_ind <- unlist(grp_list)
      mcp <- plus(X[,grp_ind], as.vector(tildeyj), method = "mc", normalize = TRUE, intercept = FALSE)
      Beta_temp[grp_ind, j] <- mcp$beta[dim(mcp$beta)[1],]
    }
    
    for (i in 1:q){ if (all(Beta_temp[ind_low[i]:ind_upp[i],j]==0)){Blabels_temp[i,j] <- 0} }
  }
  
  t <- t+1
  Betas[[t]] <- Beta_temp
  Blabels[[t]] <- Blabels_temp
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
    if (length(which(Blabels_temp[i,]!=0)) <= 5){ 
      classifications_temp[i,] <- rep(0, p)  
      means_temp[ind_low[i]:ind_upp[i],] <- 0
      covs_temp[[i]] <- 0
      cluster_nums_temp[i] <- 0
    } else if ( length(which(Beta_temp[(ind_low[i]:ind_upp[i]),Blabels_temp[i,]] !=0))/length(Beta_temp[(ind_low[i]:ind_upp[i]),Blabels_temp[i,]]) < 1/3){ 
      classifications_temp[i,] <- rep(0, p)  
      means_temp[ind_low[i]:ind_upp[i],] <- 0
      covs_temp[[i]] <- 0
      cluster_nums_temp[i] <- 0
    } else{
      temp <- Beta_temp[(ind_low[i]:ind_upp[i]),Blabels_temp[i,]]
      mod <- Mclust(t(temp), G=1:6, control=emControl(eps=0.001)) # maybe add the shape of the clusters?
      index1 <- which(Blabels_temp[i,] ==0)
      part1 <- as.data.frame(cbind(index1, 0))
      colnames(part1) <- c("index", "label")
      index2 <- which(Blabels_temp[i,] !=0)
      part2 <- as.data.frame(cbind(index2, mod$classification))
      colnames(part2) <- c("index", "label")
      est <- rbind(part1, part2)
      est <- est[order(est$index),]
      classifications_temp[i,] <- est$label
      clasind <- which(classifications_temp[i,] != 0)
      for (j in clasind){
        if (group_sizes[i]==1){
          means_temp[ind_low[i]:ind_upp[i],j] <- mod$parameters$mean[classifications_temp[i,j]]
        }else {
          means_temp[ind_low[i]:ind_upp[i],j] <- mod$parameters$mean[,classifications_temp[i,j]]
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
  
  if(sum((Beta_temp-Betas[[t-1]])^2)/(d*p) < 0.001 & t>5){break}
}

# Beta_labels records the true classifications
save(Beta, Beta_labels, Beta_cluster_num, true_means,
     Betas, Blabels, classifications, means, covariances, cluster_nums,
     file = "path/setting6.RData")
