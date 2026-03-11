#' Create a Dummy-Coded Matrix from a Dataframe
#'
#' This function transforms a dataframe containing a mix of numeric and categorical 
#' variables into a fully numeric matrix. It is specifically designed to prepare 
#' data for matrix factorization, regularized regression, or neural networks.
#'
#' @param df A data frame to be transformed.
#' @param full_rank Logical. If \code{TRUE}, the function drops the first level of each 
#'   factor to avoid the "dummy variable trap" (perfect multi-collinearity). 
#'   Set to \code{FALSE} for standard one-hot encoding. Default is \code{FALSE}.
#' @param handle_na Character. How to handle missing values. Options are \code{"na.pass"} 
#'   (leaves NAs in the matrix) or \code{"na.omit"} (removes rows with NAs). 
#'   Default is \code{"na.pass"}.
#'
#' @return A numeric matrix where all categorical variables (factors/characters) 
#'   have been expanded into binary (0/1) columns.
#' 
#' @details 
#' The function automatically converts character columns to factors before encoding. 
#' It utilizes \code{model.matrix} to ensure that the mapping of levels to columns 
#' is consistent and handles complex factor interactions if they exist in the data.
#'
#' @examples
#' \code{
#' data_mtx <- dummy_code_df(iris, full_rank = FALSE)
#' head(data_mtx)
#' }
#'
#' @importFrom stats model.matrix as.formula
dummy_code_df <- function(df, full_rank = FALSE, handle_na = "na.pass") {
  
  # 1. Validation: Ensure input is a data frame
  if (!is.data.frame(df)) {
    stop("Input 'df' must be a data frame.")
  }
  
  # 2. Pre-processing: Convert characters to factors
  # This ensures model.matrix treats them as categorical variables
  df_factored <- as.data.frame(unclass(df), stringsAsFactors = TRUE)
  
  # 3. Define the formula
  # '~ .' means include all columns. 
  # '+ 0' or '- 1' removes the global intercept, which is usually preferred 
  # for matrix factorization so that all levels of the first factor are represented.
  if (full_rank) {
    formula_choice <- stats::as.formula("~ .")
  } else {
    formula_choice <- stats::as.formula("~ . + 0")
  }
  
  # 4. Generate the model matrix
  # We use do.call to dynamically handle the na.action parameter
  mm <- stats::model.matrix(
    object = formula_choice, 
    data = df_factored, 
    na.action = handle_na
  )
  
  # 5. Clean up the '(Intercept)' column if full_rank was TRUE
  # (model.matrix adds an intercept by default when using '~ .')
  if ("(Intercept)" %in% colnames(mm)) {
    mm <- mm[, colnames(mm) != "(Intercept)", drop = FALSE]
  }
  
  return(mm)
}