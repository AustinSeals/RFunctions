#' Normalize data (in super sample) using Cosbin algorithm
#'
#' @param data data (super sample: mean of each group; if a group of 10 samples has 1000 features, 
#' super sample is a vector with size 1*1000).
#' @return a list contains normalized data and normalization factor.
#' @examples
#' cosbin(data)
cosbin <- function(data) {
  result <-data
  #row.names(result) <- 1:nrow(result)
  
  # Step 1: identify sDEG
  # mycosine(c(1, 0,0), c(1, 1, 1)) = 0.5773503
  # solve 0.5773503 = 0.99 * x - sin(arccos(x)) * sin(arccos(0.99))

  threshold <- 0.728284 #cos(arccos(0.5773503) - arccos(0.98))
  
  sDEG_output<-NULL
  
  while (max(cos_iDEG(result)) >= threshold) {
    temp_sDEG_ind <- which.max(cos_iDEG(result))
    temp_sDEG_group<-cos_iDEG_group(result)[temp_sDEG_ind]
    sDEG_output <- rbind(sDEG_output,c(max(cos_iDEG(result))
                         ,temp_sDEG_group,
                         setdiff(row.names(result),row.names(result[-temp_sDEG_ind,]))))
    result <- result[-temp_sDEG_ind, ]
    result <- totalcount(result)$data
  }
  
  sDEG_output<-sDEG_output[order(sDEG_output[,2]),]
  colnames(sDEG_output)<-c('Cosine Value','Group','Index in Complete Data')
  
  # STEP 2: identify CEG, normalize based on CEG
  # converge: # not change
  while (length(which(cos_iCEG(result) >= 0.998)) != dim(result)[1]) {
    temp_CEG_ind <- which(cos_iCEG(result) >= 0.998)
    result <- result[temp_CEG_ind,]
    result <- totalcount(result)$data
  }
  
  # Final step: normalization
  ind <- row.names(result)
  scalar <- colMeans(data[ind, ] / result)
  for (i in 1:ncol(data)) {
    data[, i] <- data[, i] / scalar[i]
  }
  
  return(list(data = data, norm_factor = scalar, CEG_index = ind, sDEG_output = sDEG_output))
}

#' convert the normalized data (in super sample) to the correct scale (with replicates)
#'
#' @param cosbin_out output of the cosbin function that kindly tells us the normalizing factors 
#' for each sample, and more
#' @param data_all original data, not pre-imputed
#' @param data_complete pre-imputed data by CONTI pre-imputation
#' @return final normalization results
#' @examples
#' cosbin_convert(cosbin_out, data_all, data_complete)
cosbin_convert <- function(cosbin_out, data_all, data_complete) {
  CEG <- data_complete[cosbin_out$CEG_index, ]
  CEG_norm <- totalcount(CEG)
  factor <- CEG_norm$norm_factor
  for (i in 1:dim(data_all)[2]) {
    data_all[, i] <- data_all[, i] /  factor[i]
  }
  return(list(data=data_all,factor=factor))
}

########################## helper functions
#' Low/high-expressed genes are filtered by their L2-norm ranks.
#' @param data Each row is a gene and each column is a sample.
#' @param thres.low The lower bound of percentage of genes to keep for Cosbin
#'     with ranked norm. The value should be between 0 and 1.
#'     The default is 0.05.
#' @param thres.high The higher bound of percentage of genes to keep for Cosbin
#'     with ranked norm. The value should be between 0 and 1.
#'     The default is 1.
data_cleaning <- function(data, thres.low = 0.05, thres.high = 1){
  if (thres.low >= thres.high) {
    stop("thres.low must be smaller than thres.high!")
  }
  
  sigNorm <- apply(data, 1, function(x) norm(matrix(x),"F") )
  Valid <- sigNorm >= quantile(sigNorm, thres.low) &
    sigNorm <= quantile(sigNorm, thres.high)
  data <- data[Valid,]
  return(data)
}

                   
#' Calculate the cosine between vector A and vector B
mycosine <- function(A, B) {
  return(sum(A * B) / sqrt(sum(A ^ 2) * sum(B ^ 2)))
}

#' Calculate the cosine between data and reference vector
cos_ref <- function(data, ref) {
  gene <- dim(data)[1]
  smp <- dim(data)[2]
  
  res <- array(0, dim = gene)
  for (i in 1:gene) {
    res[i] <- mycosine(ref, data[i, ])
  }
  return(res)
}


#' Calculate the cosine between data and CEG reference, e.g. c(1, 1, 1)
cos_iCEG <- function(data) {
  gene <- dim(data)[1]
  smp <- dim(data)[2]
  
  ref <- rep(1, smp)
  
  res <- array(0, dim = gene)
  for (i in 1:gene) {
    res[i] <- mycosine(ref, data[i, ])
  }
  return(res)
}


#' Calculate the cosine between data and closest iDEG reference, e.g. c(1, 0, 0)
cos_iDEG <- function(data) {
  gene <- dim(data)[1]
  smp <- dim(data)[2]
  
  ref <- vector("list", length = smp)
  for (i in 1:smp) {
    ref[[i]] <- rep(0, smp)
    ref[[i]][i] <- 1
  }
  
  res <- array(0, dim = gene)
  
  for (i in 1:gene) {
    temp <- NULL
    for (j in 1:smp) {
      temp <- c(temp, mycosine(ref[[j]], data[i, ]))
    }
    res[i] <- max(temp)
  }
  return(res)
}


#' Calculate the cosine between data and closest iDEG reference, e.g. c(1, 0, 0), and give their group info!
cos_iDEG_group <- function(data) {
  gene <- dim(data)[1]
  smp <- dim(data)[2]
  
  ref <- vector("list", length = smp)
  for (i in 1:smp) {
    ref[[i]] <- rep(0, smp)
    ref[[i]][i] <- 1
  }
  
  res <- array(0, dim = gene)
  
  for (i in 1:gene) {
    temp <- NULL
    for (j in 1:smp) {
      temp <- c(temp, mycosine(ref[[j]], data[i, ]))
    }
    res[i] <- which.max(temp)
  }
  return(res)
}


#' total count normalization
totalcount <- function(data) {
  scalar <- colSums(data) / mean(colSums(data))
  data <- apply(data, 2, function(x)
    x / sum(x)) * mean(colSums(data))
  
  return(list(data = data, norm_factor = scalar))
}


super_sample<-function(input,nRep){
  
  output <- NULL
  start <- 1
  for(rep in nRep){
    output <- cbind(output, rowMeans(input[, start : (start + rep - 1)])) 
    # output will have a column number the same as your group number (3 in this example)
    start <- start + rep
  }
  row.names(output)<-row.names(input)
  return(output)
}


#' Variance check
variance_check<-function(input,nRep,threshold=0.5){
  var_record<-NULL
  
  for (i in (1:dim(input)[1])){
    start <- 1
    var_feature <- NULL
    for(rep in nRep){
      var_feature <- cbind(var_feature, var(as.numeric(input[i, start : (start + rep - 1)])))
      start <- start + rep
    }
    var_record<-rbind(var_record,var_feature)
  }

  colthres<-apply(var_record, 2, function(x) quantile(x,threshold))
  index<-apply(var_record,1,function(row) all(row<colthres))
  index<-row.names(input)[index]
  
  return(index)
}