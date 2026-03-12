
library(dplyr)
library(DT)
library(stringr)
library(tidyr)
library(rlang)

#' Print a Vector with Indices
#'
#' @description Iterates through a given vector and prints each element alongside its index.
#'
#' @param my_vector A vector of any type to be printed.
#'
#' @return NULL. Prints to the console.
#' @export
#'
#' 
#' print_vec2(c("apple", "banana", "cherry"))
print_vec2 = function(my_vector){
  for (i in seq_along(my_vector)) {
    cat( i, ":", my_vector[i], "\n")
  }
}


#' Create a Simple FinalFit Data Table
#'
#' @description Wraps the \code{DT::datatable} function to create a formatted,
#' interactive HTML table optimized for finalfit output.
#'
#' @param data_ A data frame to be converted into a datatable.
#' @param label_column Character string. The name of the column containing row labels. Defaults to "label_row".
#' @param title_ Character string. Optional caption/title for the table. Defaults to NULL.
#'
#' @return A \code{DT::datatable} object.
#' @import DT
print_table_final_fit_simple  =   function(data_ ,  label_column = "label",  title_ = NULL){
 data_$label_row  =  ifelse(data_[ , label_column] != "" , T , F)
 data_[ ,  label_column] =  str_wrap(  data_[ ,  label_column] ,  width = 20)
   
 return(  
  datatable( data_  , 
             extensions = "Buttons" ,
             class = c("compact" , "row-border") ,
             options = list(
               columnDefs = list(list(targets = c( "label_row"), visible = FALSE))  , 
               pageLength = 100 ,  lengthChange = FALSE,
               scrollY = "100vh" , scrollX = TRUE , scrollCollapse = T ,  paging = F ,
               dom = 'Bfrtip', buttons = c("copy") , 
               columnDefs = list(list(className = "dt-center", targets = "_all"))
             ) ,  rownames = F,
             caption = title_  )%>%
    formatStyle(
      columns = colnames(data_) ,
      valueColumns = "label_row" , 
      borderTop = styleEqual( c(FALSE , TRUE) , c('none', '2px solid grey'))
    )
  
 )

  
}



#' Create an Alternative FinalFit Data Table with P-Value Highlighting
#'
#' @description A variation of \code{print_table_final_fit_ptest} that derives
#' \code{label_row} dynamically from an existing label column.
#'
#' @param data_ A data frame containing finalfit results.
#' @param label_column Character string. The column containing row labels. Defaults to "label".
#' @param p_column Character string. The column containing p-values. Defaults to "p".values must be character
#' @param highlite Numeric vector. The indices of the columns to highlight. Defaults to c(1).
#' @param title_ Character string. Optional caption/title for the table. Defaults to NULL.
#'
#' @return A \code{DT::datatable} object with dynamically generated row labels.
#' @import DT
#' @importFrom dplyr mutate %>%
#' @importFrom readr parse_number
print_table_final_fit_ptest2  =   function(data_ ,  label_column = "label", p_column = "p"  , highlite =  c(1)  ,  title_ = NULL){
  
  
  data_[ , "p_numeric"] = readr::parse_number(  
    as.character( data_[ ,  p_column])     
  )  %>% round(. , 3)
  
  
  data_ <- data_ %>%
    mutate(is_na_a = is.na(p_numeric))
  
  data_$label_row  =  ifelse(data_[ , label_column] != "" , T , F)
  
  data_[ ,  label_column] =  str_wrap(  data_[ ,  label_column] ,  width = 20)
  
  
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





#' Extract FinalFit Summary Statistics
#'
#' @description Extracts and splits numeric values and statistical summary strings
#' from finalfit formatted summary columns (e.g., separating mean and standard deviation).
#'
#' @param df A data frame produced by \code{finalfit::summary_factorlist}.
#' @param cols_to_split Character vector. Names of columns from which to extract
#' numeric values and statistics (e.g., removing text inside parentheses).
#'
#' @return A mutated \code{tibble} or \code{data.frame} containing the extracted statistics
#' alongside a new \code{variable_name} and \code{stat_type}.
#' @importFrom dplyr mutate fill relocate case_when na_if %>%
#' @importFrom stringr str_trim str_replace str_extract
#' @importFrom rlang := sym
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










