
source("RFunctions-main\\Cosbin_functions_advanced.R")

####### cosbatch - cosbin by batch #################
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
		supersample =  super_sample(tmp ,  nRep)
		supersample  =  data_cleaning(supersample ,  0.16 , 0.84 )
		##
		pre_check  =  variance_check(tmp ,  nRep) 
		##
		cosbin_out  = cosbin(supersample)
		iCEG =  intersect(cosbin_out$CEG_index , pre_check)
		batch_iceg[[i]] =  iCEG 
 	}

	## Get commmon ICEGs in full sample set 
	common_iCEG =  Reduce( intersect ,  batch_iceg)
	iCEG_common  =  t(data_[ ,  common_iCEG]) # orignal data
	
	##  Calculate Scalers
	iCEG_norm  =  totalcount(iCEG_common) # run tc on original data with iCEGs only
	factor_cosbin  = iCEG_norm$norm_factor
	
	## Apply Scalers to samples
	final_cosbin = t(data_)
	for ( i in 1:ncol(final_cosbin)) {
		final_cosbin[ ,i] =  final_cosbin[,i]/factor_cosbin[i]
	}

	return(list(`cosbined_data` = t(final_cosbin) , 
			`factors` =  factor_cosbin , 
			`iCEGs` = common_iCEG )
			) 





}






