
library(dplyr)
library(DT)
library(stringr)
library(tidyr)
library(rlang)

# source("H:/DESKTOP_BACKUP/RFunctions/finalfit_table_wrappers.R", echo=TRUE)



print_vec2 =  function(my_vector){
  for (i in seq_along(my_vector)) {
    cat( i, ":", my_vector[i], "\n")
  }
}




print_table_final_fit_simple  =   function(data_ ,  label_column = "label",  title_ = NULL){
 data_$label_row  =  ifelse(data_[ , label_column] != "" , T , F)
  
 return(  
  datatable( data_  , 
             extensions = "Buttons" ,
             class = c("compact" , "row-border") ,
             options = list(
               columnDefs = list(list(targets = c( label_column), visible = FALSE))  , 
               pageLength = 100 ,  lengthChange = FALSE,
               scrollY = "100vh" , scrollX = TRUE , scrollCollapse = T ,  paging = F ,
               dom = 'Bfrtip', buttons = c("copy") , 
               columnDefs = list(list(className = "dt-center", targets = "_all"))
             ) ,  rownames = F,
             caption = title_  )%>%
    formatStyle(
      columns = names(data_) ,
      valueColumns = label_column , 
      borderTop = styleEqual( c(FALSE , TRUE) , c('none', '2px solid grey'))
    )
  
 )

  
}






print_table_final_fit_ptest2  =   function(data_ ,  label_column = "label", p_column = "p"  , highlite =  c(1)  ,  title_ = NULL){
  
  
  data_[ , "p_numeric"] = readr::parse_number(  
    as.character( data_[ ,  p_column])     
  )  %>% round(. , 3)
  
  
  data_ <- data_ %>%
    mutate(is_na_a = is.na(p_numeric))
  
  data_$label_row  =  ifelse(data_[ , label_column] != "" , T , F)
  
  
  return(  
    datatable( data_  , 
               extensions = "Buttons" ,
               class = c("compact" , "row-border") ,
               options = list(
                 columnDefs = list(list(targets = c(  "label_row" , "p_numeric" , "is_na_a"), visible = FALSE))  , 
                 pageLength = 100 ,  lengthChange = FALSE,
                 scrollY = "100vh" , scrollX = TRUE , scrollCollapse = T ,  paging = F ,
                 dom = 'Bfrtip', buttons = c("copy") , 
                 columnDefs = list(list(className = "dt-center", targets = "_all"))
               ) ,  rownames = F,
               caption = title_  )%>%
      formatStyle(
        columns = names(data_) ,
        valueColumns =  "label_row" , 
        borderTop = styleEqual( c(FALSE , TRUE) , c('none', '2px solid grey'))
      ) %>%
      formatStyle(
        columns = highlite,      # The character/string columns to highlight
        valueColumns = 'p_numeric' , # The column providing the logic (< 0.05)
        backgroundColor = styleInterval( c(.05 , 0.06) , c( "yellow" ,   "orange", "transparent"))
        
        
      )
    
  )
  
  
}







extract_summary_stats <- function(df, cols_to_split) {
  
  df_out <- df %>%
    # 1. Create 'variable_name' by filling down the labels
    # (summary_factorlist leaves blanks in the label column for factor levels)
    mutate(variable_name = na_if(label, "")) %>% 
    fill(variable_name, .direction = "down") %>% 
    
    # 2. Add the stat_type indicator
    mutate(
      stat_type = case_when(
        levels == "Mean (SD)" ~ "Mean",
        levels == "Median (IQR)" ~ "Median",
        TRUE ~ "Count"
      )
    ) %>%
    # Organize columns to the front for clarity
    relocate(variable_name, stat_type, .after = label)
  
  # 3. Extract and split the numeric values for specified columns
  for (col in cols_to_split) {
    val_col  <- paste0(col, "_value") 
    stat_col <- paste0(col, "_stat") 
    
    df_out <- df_out %>%
      mutate(
        # Extract main number
        !!val_col := str_trim(str_replace(!!sym(col), "\\(.*\\)", "")),
        # Extract content within parentheses
        !!stat_col := str_extract(!!sym(col), "(?<=\\().*?(?=\\))"),
        # Convert main value to numeric
        !!val_col := as.numeric(!!sym(val_col))
      )
  }
  
  return(df_out)
}









