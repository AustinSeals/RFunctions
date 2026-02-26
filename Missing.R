
library(dplyr)
library(ggplot2)
informative.missing.bin  = function(iv ,  outcomes  , missing_cut = 50, pcut = 0.05  , data) {
  
  res =  data.frame(Variable.missing = character() ,  outcome =  character() ,  missing.effect = numeric() ,  z.value  = numeric() ,
                    pvalue  =  numeric(),
                    mean.if.missing.effect  =  numeric(),
                    mean.if.missing.pvalue =  numeric(),
                    model.pvalue= numeric(), 
                    Missing.pct =  numeric() ,  
                    Removal.Candidate = character() , 
                    Informative.Missing =  character())
  
  for (y in 1 :length(outcomes)){
    for (x in  1: length(iv) ){
      message("Testing Feature: ", x, " " ,iv[x], " " ,  round( x/length(iv)*100 , 1   ) ,  "% complete" )
      # prepare data
      data_ =  data[ , c(outcomes[y] ,  iv[x])] # dataset with 2 variables . Outcome and feature
      data_[,1] =  as.factor(data_[,1]) # convert outcome to factor 
      n =  nrow(data_)
      nmiss =  sum(is.na(data_[,2]))
      miss.pct_ =  round((nmiss/n)*100 , 2)
      if ( nmiss  == nrow(data_) | nmiss == 0 ){ # if there is 100% missing or 0%
        # append to 'res' dataframe
        to.append =  data.frame(Variable.missing = iv[x] ,  outcome =  outcomes[y] ,  missing.effect =  NA , z.value = NA , 
                                pvalue =  1,
                                mean.if.missing.effect  =  NA,
                                mean.if.missing.pvalue =  NA,
                                model.pvalue = 1  ,
                                Missing.pct =  miss.pct_ , 
                                Removal.Candidate = NA , 
                                Informative.Missing  = "NO")
        res =  rbind.data.frame(res ,  to.append)
        
        
      }else{
    
      imput_  =  mean(data_[,2] , na.rm =  T) # get median of non missing feature data
      data_$missing_ind =  ifelse(is.na(data_[,2]) ==  T ,  "YES" ,  "NO" ) # create variables indicating if feature value is missing or not
      data_$missing_imput =  ifelse(is.na(data_[,2]) ==  T ,  imput_ ,  data_[,2] ) # create variable is missing feature values replace with median
      # prepare model formula 
      lhs =  paste0(outcomes[y] , "~")
      rhs =  paste0(c("missing_ind" ,  "missing_imput") ,  collapse = "+")
      formula  =  as.formula(paste0(lhs ,  rhs)) # y~missing_ind+missing_imput
      null_formula  = as.formula(paste0(lhs , "1" ))
      # fit linear regression
      glm. = summary(glm(formula , family = "binomial" ,data =  data_))
      glm.mod = glm(formula , family = "binomial" ,data =  data_)
      null_mod  = glm(null_formula ,  family = "binomial" ,  data = data_ )
      # extract info from model coefficients table.  2nd row contains 'missing_ind' estimates
      estimate_  = glm.$coefficients[2,1] 
      z.value_  =  glm.$coefficients[2,3]
      pvalue_  = glm.$coefficients[2,4]
      global.test  = anova(null_mod, glm.mod, test = 'Chisq')
      global.p_ = global.test[2,5]
      # get mean if missing estimates
      estimate_mean_  = tryCatch(glm.$coefficients[3,1] ,  error =  function(err) NA)
      pvalue_mean_ =  tryCatch(glm.$coefficients[3,4],  error =  function(err) NA )
      
      # get missing rate
      n =  nrow(data_)
      nmiss =  sum(is.na(data_[,2]))
      miss.pct_ =  round((nmiss/n)*100 , 2)
      
      # Removal Candidate?
      removal =  ifelse(miss.pct_ >missing_cut & pvalue_ >= pcut , "YES"  , "NO" )
      informative_  =  ifelse(miss.pct_ > missing_cut & pvalue_ < pcut , "YES"  , "NO" )
      # append to 'res' dataframe
      to.append =  data.frame(Variable.missing = iv[x] ,  outcome =  outcomes[y] ,  missing.effect =  estimate_ , z.value = z.value_ , 
                              pvalue =  pvalue_ ,
                              mean.if.missing.effect  =  estimate_mean_,
                              mean.if.missing.pvalue =  pvalue_mean_,
                              model.pvalue = global.p_ ,
                              Missing.pct =  miss.pct_ , 
                              Removal.Candidate = removal ,
                              Informative.Missing  = informative_)
      res =  rbind.data.frame(res ,  to.append)
      
    }
    }
  }
  features.to.use  =  res %>% filter(Removal.Candidate == "NO" | Missing.pct ==  0 ) %>% select(Variable.missing)
  features.below.cut  =  res %>% filter(Missing.pct <=  missing_cut ) %>% select(Variable.missing)
  informative = res %>% filter(Informative.Missing == "YES")
  n.threshold  = res %>% filter( Missing.pct >  missing_cut ) %>% nrow()
  # histo
  hist  = ggplot(res, aes(x= Missing.pct)) +
    geom_histogram(alpha = 0.3) + 
    scale_x_continuous(breaks  =  seq(0 , round(max(res$Missing.pct)) , 5 )) + 
    ggtitle("Feature Missing Rate Distribution") + 
    xlab("Percent Missing")
  
  
  return(list(result.df =  res  , kept.features =  features.to.use[,1] , 
              features.below.cut =  features.below.cut[,1] , 
              informative.missing.features = informative  ,  
              missing.dist = hist,
              features.above.cut  =  n.threshold))
  
  
}




################################################################################
#  Informative missing but std feature before regression
#################################################################################



informative.missing.std  = function(iv ,  outcomes  , missing_cut = 50, pcut = 0.05  , data) {
  
  res =  data.frame(Variable.missing = character() ,  outcome =  character() ,  missing.effect = numeric() ,  z.value  = numeric() ,
                    pvalue  =  numeric(),
                    mean.if.missing.effect  =  numeric(),
                    mean.if.missing.pvalue =  numeric(),
                    model.pvalue= numeric(), 
                    Missing.pct =  numeric() ,  
                    Removal.Candidate = character() , 
                    Informative.Missing =  character())
  
  for (y in 1 :length(outcomes)){
    for (x in  1: length(iv) ){
      message("Testing Feature: ", x, " " ,iv[x] ,  round(x))
      # prepare data
      data_ =  data[ , c(outcomes[y] ,  iv[x])] # dataset with 2 variables . Outcome and feature
      data_[,1] =  as.factor(data_[,1]) # convert outcome to factor 
      n =  nrow(data_)
      nmiss =  sum(is.na(data_[,2]))
      miss.pct_ =  round((nmiss/n)*100 , 2)
      if ( nmiss  == nrow(data_) | nmiss == 0 ){ # if there is 100% missing or 0%
        # append to 'res' dataframe
        to.append =  data.frame(Variable.missing = iv[x] ,  outcome =  outcomes[y] ,  missing.effect =  NA , z.value = NA , 
                                pvalue =  1,
                                mean.if.missing.effect  =  NA,
                                mean.if.missing.pvalue =  NA,
                                model.pvalue = 1  ,
                                Missing.pct =  miss.pct_ , 
                                Removal.Candidate = NA , 
                                Informative.Missing  = "NO")
        res =  rbind.data.frame(res ,  to.append)
        
        
      }else{
        
        imput_  =  mean(data_[,2] , na.rm =  T) # get median of non missing feature data
        data_$missing_ind =  ifelse(is.na(data_[,2]) ==  T ,  "YES" ,  "NO" ) # create variables indicating if feature value is missing or not
        data_$missing_imput =  ifelse(is.na(data_[,2]) ==  T ,  imput_ ,  data_[,2] ) # create variable is missing feature values replace with median
        data_$missing_imput = as.vector(scale(data_$missing_imput)) # scale mean is missing 
        # prepare model formula 
        lhs =  paste0(outcomes[y] , "~")
        rhs =  paste0(c("missing_ind" ,  "missing_imput") ,  collapse = "+")
        formula  =  as.formula(paste0(lhs ,  rhs)) # y~missing_ind+missing_imput
        null_formula  = as.formula(paste0(lhs , "1" ))
        # fit linear regression
        glm. = summary(glm(formula , family = "binomial" ,data =  data_))
        glm.mod = glm(formula , family = "binomial" ,data =  data_)
        null_mod  = glm(null_formula ,  family = "binomial" ,  data = data_ )
        # extract info from model coefficients table.  2nd row contains 'missing_ind' estimates
        estimate_  = glm.$coefficients[2,1] 
        z.value_  =  glm.$coefficients[2,3]
        pvalue_  = glm.$coefficients[2,4]
        global.test  = anova(null_mod, glm.mod, test = 'Chisq')
        global.p_ = global.test[2,5]
        # get mean if missing estimates
        estimate_mean_  = tryCatch(glm.$coefficients[3,1] ,  error =  function(err) NA)
        pvalue_mean_ =  tryCatch(glm.$coefficients[3,4],  error =  function(err) NA )
        
        # get missing rate
        n =  nrow(data_)
        nmiss =  sum(is.na(data_[,2]))
        miss.pct_ =  round((nmiss/n)*100 , 2)
        
        # Removal Candidate?
        removal =  ifelse(miss.pct_ >missing_cut & pvalue_ >= pcut , "YES"  , "NO" )
        informative_  =  ifelse(miss.pct_ > missing_cut & pvalue_ < pcut , "YES"  , "NO" )
        # append to 'res' dataframe
        to.append =  data.frame(Variable.missing = iv[x] ,  outcome =  outcomes[y] ,  missing.effect =  estimate_ , z.value = z.value_ , 
                                pvalue =  pvalue_ ,
                                mean.if.missing.effect  =  estimate_mean_,
                                mean.if.missing.pvalue =  pvalue_mean_,
                                model.pvalue = global.p_ ,
                                Missing.pct =  miss.pct_ , 
                                Removal.Candidate = removal ,
                                Informative.Missing  = informative_)
        res =  rbind.data.frame(res ,  to.append)
        
      }
    }
  }
  features.to.use  =  res %>% filter(Removal.Candidate == "NO" | Missing.pct ==  0 ) %>% select(Variable.missing)
  informative = res %>% filter(Informative.Missing == "YES")
  
  # histo
  hist  = ggplot(res, aes(x= Missing.pct)) +
    geom_histogram(alpha = 0.3) + 
    scale_x_continuous(breaks  =  seq(0 , round(max(res$Missing.pct)) , 5 )) + 
    ggtitle("Feature Missing Rate Distribution") + 
    xlab("Percent Missing")
  
  
  return(list(result.df =  res  , kept.features =  features.to.use[,1] ,  informative.missing.features = informative  ,  missing.dist = hist ))
  
  
}


