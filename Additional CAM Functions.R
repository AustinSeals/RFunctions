##############################################################################
#               A COLLECTION OF FUNCTION TO ANALYZE CAM RESULTS 
##############################################################################
require(CAM3)
require(debCAM)
require(ggplot2)
require(dplyr)
require(tibble)
require(pcaMethods)

############################################################
# --- Create MDL plot from  CAM object---------------------
# cam.object: cam results object from CAM() or CAM3Run()
# min_size  = Size minimum MDL point on graph. Passed to ggplot
# plot.title  = what to title the plot. passed to ggtile
############################################################# 
mdl.plot  =  function(cam.object , min_size  = 5   , plot.title  =  NULL ){
  K  =  as.numeric(names(cam.object@ASestResult))
  mdls  =  vector()
  i =  1
  for (k  in K){
    i =  i + 1
    mdls[i]=cam.object@ASestResult[[as.character(k)]]@mdl
    
  }
  mdls.df = data.frame(K =  K ,  MDL  =  mdls[K])
  min_mdl  =  min(mdls.df$MDL)
  mdls.df$min_mdl  = ifelse(mdls.df$MDL == min_mdl , min_size ,  1  )
  ########################################################
  ggplot(mdls.df ,  aes(x  = as.factor(K) ,  y  =  MDL ,  size = min_mdl)) +
    geom_point() + 
    ggtitle(plot.title) + 
    scale_size(guide = "none") + 
    xlab("K -  number of sources")
}



############################################################
# --- Extract SGs at a specified cos threshold-------------
# cam.object: cam results object from CAM() or CAM3Run()
# K: the desired K number source to return results for(must be in ca.object)
# cos.thres: the threshold to be a SG
# thres.low/high: for filtering out low and highly expressed features
############################################################# 

get.SG =  function(cam.object ,  K ,  cos.thres = 0.9999999999, thres.low = 0, thres.high = 1 , PCA = FALSE){
  
  # Get Smat
  Smat = Smat(cam.object ,  K  , usingPCA  = PCA)
  MGlist_all <- cotMG(Sest = Smat, cos.thres = cos.thres, thres.low = thres.low, thres.high = thres.high)
  ## Store in   dataframe
  mg.df =  as.data.frame(MGlist_all[[2]])
  mg.df =  rownames_to_column(mg.df ,  "Feature")
  
  return(list( SG.df  = mg.df , SG.list  = MGlist_all[[1]] ))
  
}


#############################################################
# --- Return Top N SGs from each source---------------------
# sg.object: result of get.SG
# top: number of top SGs
###########################################################
balance_sg =  function(sg.object ,  top  = 30 ) {
  
  new_marker = sg.object$SG.df  %>% arrange(source , desc(cos)) %>%
    group_by(source) %>% top_n( top,  cos) %>%
    mutate(. ,  COSINE_RANK_BY_SOURCE  =  base::rank( desc(cos) , ties.method = "min"  ))
  
  return(new_marker)
}


##################################################################################################################################
# ---- Log2 -> standardize -> absolute value the data ----------
# mat: data to transform. Columns are features
# absolute: True(default) of False. Should 3rd transformation(absolute value) be executed. If not,  stop after standardization 
############################################################################################################################
matrix_transform =  function(mat , absolute =  T){
  
  if (absolute == T) {
    
    mat_abs = log2(mat+1)
    mat_abs = scale(mat_abs)
    mat_abs = abs(mat_abs)
    return(mat_abs)
    
  }else {
    mat_scaled = log2(mat+1)
    mat_scaled = scale(mat_scaled)
    return(mat_scaled)
    
  }
  
}






####################################################################
# --- Get Per Sample Average expression of SGs per Source-----------
# dat: data to process.columns are features
# sg_list: a list of vectors contains the SGs for each source
# suffix: Add suffix to new column if needed. Null is default
########################################################################

avg_sgs = function(dat ,  sg_list , suffix  = NULL) {
  results  = data.frame(Sample_ID  =  rownames(dat) )
  
  sgs.to.mean  = sg_list
  
  
  for (i in 1: length(sgs.to.mean) ){
    newcol = rowMeans(dat[ ,  sgs.to.mean[[as.character(i)]] ])
    results[ , (ncol(results) + 1) ] <- newcol
    colnames(results)[ncol(results)] <- paste0("Src ", i , " Avg" , suffix)
  }
  
  return(results[ ,  2:ncol(results)] )
  
}

#####################################################################################################################
# ------- Turn dataset contains Average SGs into long format. Useful for Plotting. Can be used for other data as well
# dat: data to pivot longer
# index.cols:columns that identify unique samples. A vector of names. They will be included in results.
# pivot.longer: columns to pivot. A vector of names

avg_to_long  =  function(dat , index.cols , pivot.cols , name = "Source" ) {
  dat_long  = dat %>% select( c(index.cols ,  pivot.cols) ) %>%
    tidyr::pivot_longer(. , cols = pivot.cols  , names_to = name , values_to = "Average Values of SGs" )
  return(dat_long)
  
  
}



####################################################################
### --- Number of SGs per Source  at various Cosine thresholds------
# cos_candidates: thresholds to display results for
# cam.object: cam results object from CAM() or CAM3Run()
# K: the desired K number source to return results for(must be in ca.object)
# feature_type: A string(Examples: "Proteins" ,"Genes" etc...). Used to title plots. Default is "Genes"
########################################################################

SG.dist =  function(cos_candidates = seq(0.85 , 0.98 , 0.01) ,  
                    cam_object ,
                    K , 
                    feature_type  = "Genes" , 
                    label_size = 5 ){
  
  sg.tallies =  vector(mode = "list" ,  length = length(cos_candidates))
  
  ### get SGs at each threshold
  for( i in 1:length(cos_candidates))   {
    sg =get.SG(cam_object ,  K , cos.thres = cos_candidates[i]  , PCA = F)
    sg.tally  = as.data.frame(table(sg$SG.df$source))
    total_sgs  = nrow(sg$SG.df)
    
    names(sg.tally) = c("Source" , "Total_SGs")
    sg.tally$COS.Threshold =  paste( "Cos>",cos_candidates[i] , ",N=" , total_sgs )
    
    sg.tallies[[i]] =  sg.tally
  }
  ### bar plot
  bar_plot  = bind_rows(sg.tallies) %>% 
    mutate_at(c("Source") , as.factor  ) %>%
    ggplot(. ,  aes(x = Total_SGs ,  y = Source ,  fill = Source )) + 
    geom_bar(stat = "identity") +
    geom_text(
      aes(label = Total_SGs), 
      ## make labels left-aligned
      hjust = 1, nudge_x = -.5 , size =  label_size 
    ) +
    #geom_vline(xintercept = 75)+ 
    #scale_x_continuous(breaks =seq(0,200 , 5) ) + 
    theme(axis.text.x = element_text(angle = 90)) + 
    ggtitle("Total SGs by Source via different Cosine Thresholds" ,  paste("Cosine Range: " ,  min(cos_candidates) , " to " , max(cos_candidates) )) + 
    facet_wrap( ~COS.Threshold )
  
  return(bar_plot)
}



######################################################################
## ---------- Heatmap of Smat with Column Standardized
# dat: Data for heatmap. Features are columns
# font_size   =  font size of column labels
#####################################################################
smat_heatmap  =  function(dat  ,  font_size  =  5){
  
  Heatmap(scale(  dat  ) ,name  = "Source Matrix Heatmap"  ,
          row_order = NULL ,  column_names_gp=gpar(fontsize=font_size)     )
  
}







# Helper function function for predictive screening
screen  =  function(predictors_  , data_ , outcome_ , covars_ = NULL   , family_ = family) {
  names(data_) =  make.names(names(data_))
  terms  =  paste0(c(make.names(predictors_) , make.names(covars_)) ,  collapse = "+")
  
  
  formula_  =  as.formula(paste0(outcome_, "~" ,terms))
  model_  =  glm(formula = formula_ ,  family = family_ ,  data =  data_  )
  model_smry  =  summary(model_)
  
  model_coefs  =  as.data.frame(model_smry$coefficients)
  model_coefs = model_coefs[2,]
  model_coefs = tibble::rownames_to_column(model_coefs ,  "Predictor")
  
  return(model_coefs)
}


############################################################################
# --- Perform univariate predictor screening --------------
# predictors =  a  character list of predictors to model with the outcome. numeric only
# data_  = input data for model
# outcome: outcome of model. character. numeric data only
# covars_  = a list variables to adjust for. Added to model
# family_ =  model family type. Common are gaussian(default) and "binomial"
##############################################################################
predictor_screening =  function(predictors  , data , outcome , covars = NULL   , family = "gaussian") {
  
  res = lapply(predictors, screen,  data_ = data, outcome_= outcome , 
               covars_ = covars  , family_=family ) %>% bind_rows()
  res$Predictor =  predictors
  res$`Pr(>|z|)` =  round(res$`Pr(>|z|)` ,  8 )
  res$FDR = round(p.adjust(res$`Pr(>|z|)` ,  method = "fdr"  ) , 8)
  
  return(res)
}



