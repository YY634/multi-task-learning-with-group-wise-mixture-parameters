
#the input of sort_grouping is a dataframe with two variables, "index" and "label".
sort_grouping <- function(df){
  attach(df)
  temp <- df[order(label, index),]
  labels <- sort(unique(temp$label))
  num <- length(labels)
  L <- list()
  for (i in 1:num){
    ind <- which(temp$label == labels[i])
    L[[i]] <- temp$index[ind]
  }
  L <- L[order(sapply(L, function(x) x[1], simplify=TRUE))]
  group_sizes <- rep(0, num)
  for(l in 1:num){
    group_sizes[l] <- length(L[[l]])
  }
  ordered_index <- unlist(L)
  new_labels <- rep(1:num, group_sizes)
  new_df <- as.data.frame(cbind(ordered_index, new_labels))
  colnames(new_df) <- c("index", "label")
  return(new_df)
}


TPR <- function(tr_sparsity, hat_sparsity){
  real_pos <- length(which(tr_sparsity != 0))
  detected_pos <- length(which(tr_sparsity !=0 & hat_sparsity != 0))
  rate <- detected_pos/real_pos
  return(rate)
}
TNR <- function(tr_sparsity, hat_sparsity){
  real_neg <- length(which(tr_sparsity == 0))
  detected_neg <- length(which(tr_sparsity ==0 & hat_sparsity == 0))
  rate <- detected_neg/real_neg
  return(rate)
}
ACC <- function(tr_sparsity, hat_sparsity){
  detected_pos <- length(which(tr_sparsity !=0 & hat_sparsity != 0))
  detected_neg <- length(which(tr_sparsity ==0 & hat_sparsity == 0))
  total <- dim(as.matrix(tr_sparsity))[1] * dim(as.matrix(tr_sparsity))[2]
  rate <- (detected_pos+detected_neg)/total
  return(rate)
}


#setting 1.   
n <- 1800;  d <- 3500;  q <- 700;  p <- 200;
load("path/setting1.RData")
T <- length(classifications)
RES <- matrix(0, T, 5)
for (t in 1:T){
  # to order the labels in classifications[[t]]
  labels_temp <- matrix(0, q, p)
  class_temp <- classifications[[t]]
  for (i in 1:q){
    index2 <- which(class_temp[i,] !=0)
    if (length(index2)==0){}else{
      index1 <- which(class_temp[i,] ==0)
      part1 <- as.data.frame(cbind(index1, 0))
      colnames(part1) <- c("index", "label")
      part2_temp <- as.data.frame(cbind(index2, class_temp[i,index2]))
      colnames(part2_temp) <- c("index", "label")
      part2 <- sort_grouping(part2_temp)
      colnames(part2) <- c("index", "label")
      est <- rbind(part1, part2)
      est <- est[order(est$index),]
      labels_temp[i,] <- as.vector(est$label)
    }
  }
  rm(i, index1, index2, part1, part2, part2_temp, est)
  # to order the labels in Beta_labels (truth)
  true_labels <- matrix(0, q, p)
  for (i in 1:q){
    index4 <- which(Beta_labels[i,] !=0)
    if (length(index4)==0){}else{
      index3 <- which(Beta_labels[i,] ==0)
      part3 <- as.data.frame(cbind(index3, 0))
      colnames(part3) <- c("index", "label")
      part4_temp <- as.data.frame(cbind(index4, Beta_labels[i, index4]))
      colnames(part4_temp) <- c("index", "label")
      part4 <- sort_grouping(part4_temp)
      colnames(part4) <- c("index", "label")
      tr <- rbind(part3, part4)
      tr <- tr[order(tr$index),]
      true_labels[i,] <- as.vector(tr$label)
    }
  }
  rm(i, index3, index4, part3, part4, part4_temp, tr)
  
  eva <- rep(0, 5)
  eva[1] <- sum((Beta-Betas[[t]])^2)/(d*p)
  eva[2] <- ACC(Beta_labels, Blabels[[t]])
  eva[3] <- length(which(true_labels !=  labels_temp))/(q*p) 
  eva[4] <- length(which(Beta_cluster_num != cluster_nums[[t]]))/q
  eva[5] <- sum((true_means - means[[t]])^2)/(d*p)
  RES[t,] <- eva
}

save(RES,  file = "path/eval_setting1.RData")




