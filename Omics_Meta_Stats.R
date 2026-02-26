
library(pastecs)
library(dplyr)
library(tibble)
library(pcaMethods)
library(ggplot2)
library(patchwork)
source('H:/DESKTOP_BACKUP/RFunctions/Missing.R') 

###################################################################################################################################
### Generate Meta Data stats for both samples and Features. Samples shold be rows(identified in row names) and Features are Columns
# sample.col.name  :  Name to apply to the Sample ID column in the results
# feature.col.name :  Name to apply to the Feature ID column in the results
###################################################################################################################################
generate.meta.stats  =  function(data ,  sample.col.name  = "ID" ,  feature.col.name  = "Feature" , pca_miss_thresh = 100){
  
  features.n  =  ncol(data)
  samples.n  =  nrow(data)
  
  feature.meta  =as.data.frame( t(stat.desc(data)))
  feature.meta$Pct.Missing =  round(  (feature.meta$nbr.na /  samples.n) , 3 ) *100
  feature.meta =  feature.meta %>% relocate(Pct.Missing)
  feature.meta =  rownames_to_column(feature.meta  ,  feature.col.name)
  R.Var.Name.Format  =  make.names((feature.meta[,1]))
  feature.meta =  tibble::add_column(feature.meta ,   R.Var.Name.Format ,  .after = feature.col.name)
  feature.meta$index  =  seq(1 , features.n)
  
  low.thresh = length(which(feature.meta$Pct.Missing <=50))
  medium.thresh = length(which(feature.meta$Pct.Missing > 50 &  feature.meta$Pct.Missing <67   ))
  high.thresh = length(which(feature.meta$Pct.Missing >= 67))
  #############
  sample.meta  =as.data.frame(  t(stat.desc( t(data)   ))   )
  sample.meta$Pct.Missing =  round(  (sample.meta$nbr.na /  features.n) , 3 ) *100
  sample.meta =  sample.meta %>% relocate(Pct.Missing)
  sample.meta =  rownames_to_column(sample.meta  ,  sample.col.name)
  ## LOF based on stats 
  lof_input  =  sample.meta %>%
    select(Pct.Missing ,  median , mean , coef.var , min , max) %>%
    select(where(~n_distinct(.) > 1))
  lof_input = scale(lof_input)
  LOF  =  dbscan::glosh(lof_input)
  LOF_outlier  = ifelse(LOF > 0.9  , "YES" , "NO"  )
  ## leverage
  pca_input_features  =  which(feature.meta$Pct.Missing <=pca_miss_thresh & feature.meta$std.dev != 0 )
  pca.fit =  pca(as.matrix(scale(  log2(data[ , pca_input_features ])  )) , nPcs = 20 , method = "ppca"   , scale = "none" ,  center = FALSE   ) # Run ppca
  r2 = pca.fit@R2cum
  
  lev = as.data.frame(pcaMethods::leverage(pca.fit)) # get sample leverage 
  names(lev) = c("Leverage") # formatting results
  lev$index =  seq(1 , nrow(lev))
  lev$Label =  row.names(lev)
  lev$Leverage.Percentile  =  percent_rank(lev$Leverage) # calculate percentile of leverge for each sample
  lev.med =  median(lev$Leverage) # get median
  lev$Outlier.ID =  ifelse(lev$Leverage > 4*lev.med  ,  lev$Label ,  NA) # identify outliers as having leverage more than 4 times the median leverage
  lev$Outlier  =  ifelse(is.na(lev$Outlier.ID) ,  "NO" ,  "YES")
  #####
  sample.meta =  cbind.data.frame(sample.meta ,  lev[ ,  c("index" , "Leverage", "Leverage.Percentile" , "Outlier.ID" , "Outlier")])
  sample.meta =  cbind.data.frame(sample.meta , scores(pca.fit)[ , 1:2]  )
  
 return(list(feature.meta =  feature.meta ,  sample.meta = sample.meta , low.thresh = low.thresh ,  high.thresh = high.thresh ,  medium.thresh = medium.thresh  ))
  
}

###########################################################################################################################################
### Generate informative missing results. For the 
#  data : Samples should be rows(identified in row names) and Features are Columns
#  iv_:  are the indices for the features you wants to test in 'data'. 
# outcome.vector: is a seprate vector containing outcome variable
# outcome.name :  name to help identify  what the outcome variable was in the analysis 
############################################################################################################################################
generate.informative.missing  =  function(iv_ ,  outcome.vector , outcome.name  , missing_cut = 50, pcut = 0.05  ,   data) {
  
  missing.test.df_ =  cbind.data.frame(outcome.vector , data)
  names(missing.test.df_)[1] =  outcome.name
  names(missing.test.df_) =  make.names(names(missing.test.df_))
  missing.test.df_[,1] =  as.factor(missing.test.df_[,1])

  
  missing.reg  = informative.missing.bin( iv = names(missing.test.df_)[iv_ + 1], outcomes =  outcome.name ,  data = missing.test.df_) 
  missing.reg =  missing.reg$result.df
  fdr.pvalue = p.adjust(missing.reg$pvalue ,  "fdr") 
  missing.reg =  add_column(missing.reg , fdr.pvalue ,  .after = "pvalue"  )
  missing.reg$informative.missing =  ifelse(missing.reg$fdr.pvalue <  pcut ,  "YES" , "NO")
  return(missing.reg[ ,  c(1:6 , 13)])
  
}


##############################################################################################################################
### Combind Feature Meta Stats and Missing regression data. Supply meta stats first and missing regrssion results second
################################################################################################################################
combine.feature.meta.data  =  function(feature.meta.df  ,  informative.missing.df) {
  
  return(cbind.data.frame(feature.meta.df ,  informative.missing.df[ , 2: ncol(informative.missing.df) ]   ))
  
  
}




####################################################################################################
##                   Missingness by Batch
#####################################################################################################################

### batch missingness filter
miss = function(.x ) {
  sum(is.na(.x)) / length(.x)
}

meets_thresh  = function(x , cut = 0.5) {
  all(x  <= cut)
  
}


meets_thresh67  = function(x , cut = 0.67) {
  all(x  < cut)
  
}




batch_miss_filter = function( dat_ , batch_col   , select_cols = 1 ) {
  
  miss_smry  = dat_ %>% group_by({{batch_col}}) %>%
    summarise(across( all_of(select_cols) ,  miss)   )
  
  good_vars  = miss_smry[,sapply(miss_smry, meets_thresh )] %>% colnames()
  #good_vars = good_vars %in% select_cols
  
  return(  list( `good_vars` = good_vars ,  `smry` =miss_smry )    )
  
  
}



####
batch_miss_filter67 = function( dat_ , batch_col  , select_cols = 1 ) {
  miss_smry  = dat_ %>% group_by({{batch_col}}) %>%
    summarise(across( all_of(select_cols) ,  miss)   )
  
  good_vars  = miss_smry[,sapply(miss_smry, meets_thresh67 )] %>% colnames()
  #good_vars = good_vars %in% select_cols
  
  return(  list( `good_vars` = good_vars ,  `smry` =miss_smry )    )
  
  
  
  
}







##### Example Work Flow: {do not run}
# 
# lad45.meta.stats =  generate.meta.stats(marker.data ,  sample.col.name = "CaseID" ,  feature.col.name = "Protein")
# lad45.feature.info.missing  =  generate.informative.missing(iv_ = 1:ncol(marker.data)  ,  informative.test.df$Sample.Type ,  outcome.name = "Normal.or.Lesion" ,  data = marker.data)
# full.feature.meta.data  = combine.feature.meta.data(lad45.meta.stats$feature.meta  , lad45.feature.info.missing   )
# 
# # Plot sample missing rate
# sample.miss.hist  = ggplot(lad45.meta.stats$sample.meta, aes(x=Pct.Missing , color= NULL ,  fill  =  NULL)) +
#   geom_histogram(alpha = 0.3) +
#   scale_x_continuous(breaks  =  seq(0 , 100 ,  5 )) +
#   ggtitle("Sample Missing Rate Distribution")
# 
# 
# # Plot sample Mean Expression
# sample.mean.plot  = ggplot(lad45.meta.stats$sample.meta , aes(x =index ,  y = mean , label = CaseID ,  color  =  NULL )) +
#   geom_point() +
#   geom_text(size =  3) +
#   ggtitle("Mean Sample Protein Intensity")+
#   xlab("Sample Index")
# 
# 
# # Plot sample Leverage
# sample.lev.plot  = ggplot(lad45.meta.stats$sample.meta , aes(x =index ,  y = Leverage , label = Outlier.ID ,  color  =  NULL )) +
#   geom_point() +
#   geom_text() +
#   ggtitle("Sample PPCA Leverage" , paste("Only includes features with <=50% missing (" , lad45.meta.stats$low.thresh ,")") )+
#   xlab("Sample Index")
# 
# # feature missing levels
# feature.missing.levels  = data.frame(`Features Missing <=50%` = lad45.meta.stats$low.thresh ,
#            `Features Missing >50 to <67` = lad45.meta.stats$medium.thresh ,
#            `Features Missing >67%` = lad45.meta.stats$high.thresh  )
# names(feature.missing.levels ) = c("Features Missing <=50%" , "Features Missing >50 to <67" , "Features Missing >67%" )
# 
# 
# # patch it all together into one figure
# ( (sample.miss.hist + sample.lev.plot) /  (sample.mean.plot + gridExtra::tableGrob(t(feature.missing.levels) )) ) + plot_annotation(
#   title = paste("LAD45 | Features:" , nrow(lad45.meta.stats$feature.meta)   ,   " | Samples:" ,  nrow(lad45.meta.stats$sample.meta))
# )
# 
# 
##### END example work flow 


