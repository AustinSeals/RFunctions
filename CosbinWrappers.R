
library(dplyr)
source("RFunctions-main\\Cosbin_functions_advanced.R")

###############

sample_means  =  function(data_ , batch_ , runorder_ , description_  = "Mean Intensity per Sample") {
  
  plot_dat =  cbind.data.frame(SampleID = rownames(data_)  , 
                               Batch = batch_ , 
                               RunOrder  = runorder_  , 
                               Mean  = rowMeans(data_ ,  na.rm = T),
                               Description = description_ ) 
  
  ggplot(plot_dat , aes(x =as.factor(RunOrder) , y = Mean , color  =  Batch )) +
      geom_point() + 
      xlab("Sample Run Order") + 
      ggtitle(description_)
  
  return(plot_dat )

}




############## GET INDEXES OF MISSING VALUES

missing_loc  =  function(data_){
  
  values_to_impute = which(is.na(data_) ,  arr.ind = T)
  return(values_to_impute)
  
}

# define TC function
totalcount <- function(data) {
  scalar <- colSums(data , na.rm = T) / mean(colSums(data , na.rm = T))
  data <- apply(data, 2, function(x)
    x / sum(x , na.rm = T)) * mean(colSums(data , na.rm = T) , na.rm = T)
  
  return(list(data = data, norm_factor = scalar))
}


########## TOTALCOUNT NORMALIZATION
tc_norm =  function(data_){
  # transpose input data 
  data_t = t(data_)

  # set transposed data into another variable 'data_norm_totalcount' 
  data_norm_totalcount = data_t # includes missing 
  data_norm_totalcount[is.na(data_norm_totalcount)]<- 0 # Impute by 0
  
  ## run tc
  totalcount_norm <- totalcount(data_norm_totalcount)
  factor_totalcount <- totalcount_norm$norm_factor # Get the scalers from 0-imputed data
  
  ## use scalers to normalize 'data_'(after_remove) that has NAs
  for (i in 1:dim(data_norm_totalcount)[2]) {
    data_norm_totalcount[, i] <- data_t[, i]/factor_totalcount[i] # Use the scalers on non-imputed data
  }
  
  return(
    list(Data = t(data_norm_totalcount) ,  factors  =  factor_totalcount )
    )
  
  
}


############ NIPALS IMPUTATION
nipals_impute =  function(data_ ,  pcs = 20) {
  
  log_data <- log2(data_+1) # for non-CAM methods
  #Performing imputation using the PCA method with "nipals" algorithm, using the function 
  imputedNip_log <- completeObs(pca(log_data, nPcs= pcs, method="nipals"))
  #Converting the log transformed imputed data back to the original scale by applying the inverse transformation
  imputedNipals <- (2 ^ imputedNip_log)-1
  
  return(imputedNipals)
  
}

########## NIPALS by Batch 
nipals_batch  =  function(data_ ,batch_){
  sample_order  =  rownames(data_)
  genes =  colnames(data_) 
  data =  cbind.data.frame(Batch =  batch_ ,  data_)
  
  batch_list =  split(  data ,  data$Batch)
  imputed_list = vector(mode = "list" ,  length  =  length(batch_list) )	
  
  for (i in 1:length(batch_list)){
    message( paste("Batch:" ,  as.character(i))           ) 
    tmp =  batch_list[[i]]
    tmp  =  tmp[ ,  genes]    
    
    imputed =  nipals_impute(tmp)
    imputed_list[[i]] =  imputed
    
  }
  
  final_df= rlist::list.rbind(imputed_list)
  final_df = final_df[sample_order , ]
  return(final_df)
}

apply_factors  =  function(data_ ,  factors_ ){
  data =  t(data_)
  for (i in 1:dim(data)[2]) {
    data[, i] <- data[, i]/factors_[i] 
  }
  return(t(data))
}



############ COSBIN
# data_ : rows are samples ,  columns are genes ,  with sample ids as rownames
# group : group variables to sort samples before cosbin
cosbin_norm  =  function(data_ , group) {
  
  sample_order  = rownames(data_) # origninal sample order
  genes  = colnames(data_)
  
  data =  cbind.data.frame(Grp  = group , data_  ) # combine xpression data and group data
  
  data = data[order(data$Grp) , ] # sort by group
  nRep =  table(data$Grp)
  
  data_norm_input<-t(data[, genes]   ) # transpose - samples are columns
  
  supersample<-super_sample(data_norm_input,nRep)
  supersample<-data_cleaning(supersample,0.16,0.84)
  
  # Index of high within-batch-variance feature
  iCEG_precheck<-variance_check(data_norm_input,nRep)
  
  cosbin_out<-cosbin(supersample)
  iCEG<-intersect(cosbin_out$CEG_index,iCEG_precheck)
  
  # iCEG totalcount
  iCEG_common<-data_norm_input[iCEG,]
  iCEG_norm <- totalcount(iCEG_common)
  factor_cosbin <- iCEG_norm$norm_factor
  
  data_norm<-data_norm_input
  for (i in 1:dim(data_norm_input)[2]) {
    data_norm[, i] <- data_norm_input[, i]/factor_cosbin[i] 
  }
  
  # Check the results
  test<-data_norm[iCEG,]
  test<-colSums(test) # Identical? G2G!
  print("Test: Are they identical?")
  print(test)
  
  final =  data_norm[ ,  sample_order]
  
  return( 
    list( Data = t(final)  , factors = factor_cosbin[sample_order] )
        )
  
  
}



