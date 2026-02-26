

#########################
## Load Helper Functions
#######################
scripts_path  = "RFunctions-main/"

source(paste0(scripts_path , "Omics_Meta_Stats.R")   )
source(paste0(scripts_path , "CosbinWrappers.R")   )
source(paste0(scripts_path ,  "Cosbatch.R"))


#######################################
### Full Pipeline 
########################################

prot_pipeline =  function(BatchID  , PhenoType    , expr_data ,  mean_only = T , model = NULL  ) {
  
  
  ## 1. Get <50% missing and <67% features in each batch 
  
  message("1. Calculating which features have <50% and <67% missing in each batch")
  feature_miss_batch  =  batch_miss_filter( cbind.data.frame(BatchID = BatchID ,  expr_data) ,  batch_col = `BatchID` , 
                                            select_cols  =  colnames(expr_data)   )
  message(paste("Proteins with Less then 50% missing in each batch:"  , length(feature_miss_batch$good_vars)   ))
  
  ##
  feature_miss_batch67  =  batch_miss_filter67( cbind.data.frame(BatchID = BatchID ,expr_data) ,  batch_col = `BatchID` , 
                                                select_cols  =  colnames(expr_data)   )
  message(paste("Proteins with Less then 67% missing in each batch:"  , length(feature_miss_batch67$good_vars)   ))
  message("Done with missing analysis")
  
  expr50  =  expr_data[ , feature_miss_batch$good_vars]
  expr67 =  expr_data[ ,  feature_miss_batch67$good_vars]
  
  
  impute_locales50  =  missing_loc(expr50)
  impute_locales67  =  missing_loc(expr67)
  
  ############ START OF PIPELINE ############################################
  
  ## 2.  TC and Impute 50% data 
  
  message("2. Total Count Normaliztion of 50% Data them impute")
  expr50_normed  = tc_norm( expr50  )
  # Impute
  expr50_tc_impute  =  nipals_batch(expr50_normed$Data ,  BatchID)
  print(dim(expr50_tc_impute))
  
  ## 3. Cosbin on 50% data
  message("3. Cosbin Normalize the 50% Data")
  expr50_impute_cosbind = cosbatch( expr50_tc_impute , PhenoType   , BatchID   )
  
  
  ## 4.  Apply 50% TC Factors to 67% Data and Impute 
  
  message("4. Apply 50% TC Factors to 67% Data and Impute ")
  expr67_tc = apply_factors(expr67  , expr50_normed$factors)
  expr67_tc_impute  =  nipals_batch(expr67_tc , BatchID)
  
  ## 5. Apply 50% Cosbin Factors to 67% TC
  
  message("5. Apply 50% Cosbin Factors to 67% TC normed data")
  expr67_tc_impute_cosbin =  apply_factors(expr67_tc_impute ,  expr50_impute_cosbind$factors)
  
  ## 6. Run Combat(sva::Combat) - mean only adjusment
  
  message("6. Run Combat(sva::Combat) - mean only adjusment")
  expr67_combat  = sva:: ComBat(t(expr67_tc_impute_cosbin) ,  
                                BatchID ,
                                mod= model , 
                                par.prior = T , 
                                mean.only = mean_only)
  
  # make matrix non-negative
  min(expr67_combat)->min_combat
  expr67_combat =  expr67_combat + abs(min_combat)
  
  
  ## 7. Final Nipals
  
  message("7. Final Nipals")
  expr67_final_nipals  =  t(expr67_combat)
  expr67_final_nipals[impute_locales67] = NA
  expr67_final_nipals =  nipals_impute(expr67_final_nipals)
  
  
  resObj  = list(vars50 =  feature_miss_batch$good_vars,
                 vars60 =  feature_miss_batch67$good_vars , 
                 tc_factors50  = expr50_normed$factors , 
                 cosbin_factors50  =  expr50_impute_cosbind$factors ,
                 iCEGs = expr50_impute_cosbind$iCEGs,
                 expr50_tc_impute  =  expr50_tc_impute  , 
                 expr50_tc_impute_cosbind = expr50_impute_cosbind , 
                 expr67_tc_impute = expr67_tc_impute,
                 expr67_tc_impute_cosbin = expr67_tc_impute_cosbin,
                 expr67_tc_impute_cosbin_combat = t(expr67_combat),
                 expr67_final_nipals = expr67_final_nipals
                 )  
  
  return(resObj)
  
  
}


#prot_pipeline(data_cleaned$BatchID ,  data_cleaned$PATHCLASS4 ,  data_cleaned[ , 27:4000] ) ->test_pipeline





pipeline_boxplot=  function(matrix , annot  ,  sample_col ,  batch_col  , factors) {
  
  annot$order = 1:nrow(annot)
  color_list = sample_annotation_to_colors(annot  ,  factor_columns = c(batch_col,  factors) , 
                                           numeric_columns = 'order')
  
  proBatch::matrix_to_long(t(matrix) , 
                           sample_annotation = annot , 
                           sample_id_col = sample_col) -> long
  
  longLog =  long
  longLog$Intensity  = log2(longLog$Intensity + 1)
  
  plot_boxplot(longLog, annot,
               batch_col = batch_col, color_scheme = color_list[[batch_col]] , sample_id_col = sample_col  )
  
  
}





pipeline_dendro=  function(matrix , annot  ,  sample_col ,  batch_col  , factors) {
  
  annot$order = 1:nrow(annot)
  color_list = sample_annotation_to_colors(annot  ,  factor_columns = c(batch_col,  factors) , 
                                           numeric_columns = 'order')
  
  plot_hierarchical_clustering( t(matrix) , 
                                sample_annotation = annot ,
                                sample_id_col = sample_col  , 
                                color_list = color_list , 
                                factors_to_plot = c("BatchID" , factors) , 
                                distance = 'euclidean', agglomeration = 'complete', 
                                label_samples = FALSE)
  
  
}



pipeline_PVCA  =  function(matrix , annot  ,  sample_col ,  batch_col  , factors) {
  mat = t(matrix)
  rownames(annot) =  NULL
  plot_PVCA(mat , annot,
            technical_factors = c(batch_col),
            biological_factors = c(factors) ,  
            sample_id_col = sample_col)


}






























