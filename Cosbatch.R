
source_file  = " RFunctions-main/Cosbin_functions_advanced.R"
tryCatch(
  source(source_file) , 
  error = function(e) {
    message(paste("Error: Failed to source", source_file))
  }
)

####### cosbatch - cosbin by batch #################

#' Peform cosbin across batches
#'  @description
#'  Cosbin across batches by doing cosbin within each batch then identifying common iCEGs.
#'
#'
#' @param data_	a data matrix. rows are samples ,  columns are genes ,  with sample ids as rownames
#' @param grp_	a vector. group indicator to sort samples before cosbin
#' @param batch_ a vector. the batch id for each sample
#'
#'  @return A list containing:
#' \itemize{
#'   \item \code{Data}: The cosbin normalized data.
#'   \item \code{factors}: The scalind factors for each sample
#' }
cosbatch  =  function(data_ , grp_ ,  batch_){
	sample_order  =  rownames(data_)
	genes =  colnames(data_) 
	data =  cbind.data.frame(Group = grp_ ,  Batch =  batch_ ,  data_) # make sure rownames are preserved
	message(paste("Rownames are Preserved:" , as.character( identical(rownames(data) , rownames(data_)    )   ) ) )
	
	## Within Batch Cosbin - to identify common iCEG genes
	batch_list =  split(  data ,  data$Batch)
	batch_iceg = vector(mode = "list" ,  length  =  length(batch_list) )	

	for (i in 1:length(batch_list)){
	  message( paste("Batch:" ,  as.character(i))           ) 
		tmp =  batch_list[[i]]
		nRep  =  table(tmp$Group) # group counts
		tmp =  t(tmp[order(tmp$Group)  , genes]) # sort by group the transpose
		##
		print(dim(tmp))
		message("Super Sample Calculation")
		supersample =  super_sample(tmp ,  nRep)# function from 'Cosbin_functions_advanced.R'
		supersample  =  data_cleaning(supersample ,  0.16 , 0.84 )# function from 'Cosbin_functions_advanced.R'
		##
		pre_check  =  variance_check(tmp ,  nRep) # function from 'Cosbin_functions_advanced.R'
		##
		cosbin_out  = cosbin(supersample) # function from 'Cosbin_functions_advanced.R'
		iCEG =  intersect(cosbin_out$CEG_index , pre_check)
		batch_iceg[[i]] =  iCEG 
 	}

	## Get commmon ICEGs in full sample set 
	common_iCEG =  Reduce( intersect ,  batch_iceg)
	iCEG_common  =  t(data_[ ,  common_iCEG]) # orignal data
	
	##  Calculate Scalers
	iCEG_norm  =  totalcount(iCEG_common) # run tc on original data with Cosbin derived iCEGs only
	factor_cosbin  = iCEG_norm$norm_factor
	
	## Apply Scalers to samples
	final_cosbin = t(data_)
	for ( i in 1:ncol(final_cosbin)) {
		final_cosbin[ ,i] =  final_cosbin[,i]/factor_cosbin[i]
	}

	return(list(`cosbin_data` = t(final_cosbin) , 
			`factors` =  factor_cosbin , 
			`iCEGs` = common_iCEG )
			) 

}





#' Perform Nipals Imputation.
#'
#' @description
#' Perform nipals imputation on a data matirx. Data is log2 transformed before imputaiton
#' then transformed back to orgininal scale.
#'
#'
#' @param data_ a numeric data matirx
#' @param pcs the number of components to keep for nipals - ppca impuation
#'
#' @return the imputed data matrix
#'
#' @importFrom pcaMethods completeObs pca
#' nipals_impute(data_ , pcs = 20)
#' @export
nipals_impute =  function(data_ ,  pcs = 20) {

  log_data <- log2(data_+1) # for non-CAM methods
  #Performing imputation using the PCA method with "nipals" algorithm, using the function
  imputedNip_log <- pcaMethods::completeObs(pcaMethods::pca(log_data, nPcs= pcs, method="nipals"))
  #Converting the log transformed imputed data back to the original scale by applying the inverse transformation
  imputedNipals <- (2 ^ imputedNip_log)-1

  return(imputedNipals)

}


#' Perform within batch nipals imputation
#'
#' @param data_  the data matrix. samples are rows and genes are columns
#' @param batch_ a vector. the sample batch ids. the length  should be same as nrow(data_)
#'
#' @return the imputed data matrix
#'
#' @importFrom rlist list.rbind
#' @export
#' nipals_batch(data_ , pcs = 20)
#' @export
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








