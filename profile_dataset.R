
#' Profile the distribution and type of all variables in a dataset
#'
#' @param df A dataframe to analyze.
#' @return A tidy dataframe summarizing each column's characteristics.
profile_dataset <- function(df) {
  
  results_list <- list()
  
  for (col_name in names(df)) {
    column_data <- df[[col_name]]
    
    # 1. Basic Metadata
    n_total <- length(column_data)
    n_miss  <- sum(is.na(column_data))
    n_unique <- length(unique(na.omit(column_data)))
    data_type <- class(column_data)[1]
    
    # 2. Heuristic Classification
    # This helps decide if a numeric var should actually be a factor (e.g., "Sex" coded 1/2)
    classification <- "Continuous"
    if (n_unique == 1) {
      classification <- "Constant (Zero Variance)"
    } else if (n_unique == 2) {
      classification <- "Binary"
    } else if (is.character(column_data) || is.factor(column_data)) {
      classification <- "Categorical"
    } else if (n_unique < 10 && is.numeric(column_data)) {
      classification <- "Categorical (Numeric Coded?)"
    } else if (n_unique == n_total) {
      classification <- "Likely ID/Unique Key"
    }
    
    # 3. Descriptive Stats (Conditional on type)
    summary_val <- ""
    if (is.numeric(column_data)) {
      mu <- round(mean(column_data, na.rm = TRUE), 2)
      sigma <- round(sd(column_data, na.rm = TRUE), 2)
      summary_val <- paste0("Mean: ", mu, " (SD: ", sigma, ")")
    } else {
      # For categorical, show the top level and its frequency
      top_level <- names(sort(table(column_data), decreasing = TRUE))[1]
      count <- sort(table(column_data), decreasing = TRUE)[1]
      pct <- round((count / n_total) * 100, 1)
      summary_val <- paste0("Mode: '", top_level, "' (", pct, "%)")
    }
    
    results_list[[col_name]] <- data.frame(
      Variable = col_name,
      Class = data_type,
      Inferred_Type = classification,
      Unique_Values = n_unique,
      Missing_N = n_miss,
      Missing_Pct = round((n_miss / n_total) * 100, 2),
      Summary_Stats = summary_val,
      stringsAsFactors = FALSE
    )
  }
  
  final_profile <- do.call(rbind, results_list)
  rownames(final_profile) <- NULL
  return(final_profile)
}





#' Automatically convert likely categorical variables to factors
#'
#' @param df A dataframe.
#' @param max_levels The maximum number of unique values a numeric variable 
#'                   can have to be considered a factor (default is 10).
#' @param exclude A character vector of column names to skip (e.g., ID columns).
#' @return A dataframe with updated column classes.
auto_factorize <- function(df, max_levels = 10, exclude = NULL) {
  
  # 1. Run our profile function to get the "Inferred_Type"
  profile <- profile_dataset(df)
  
  # 2. Identify variables that should be factors
  # Criteria: Is already a character OR is Binary OR is Numeric with few levels
  vars_to_factor <- profile$Variable[
    profile$Inferred_Type %in% c("Binary", "Categorical (Numeric Coded?)", "Categorical")
  ]
  
  # 3. Filter out excluded variables and those with zero variance (constants)
  vars_to_factor <- setdiff(vars_to_factor, exclude)
  vars_to_factor <- setdiff(vars_to_factor, profile$Variable[profile$Inferred_Type == "Constant (Zero Variance)"])
  
  # 4. Filter by the user-defined max_levels threshold for numeric variables
  # This prevents accidentally factorizing things like 'Age' if it only has 12 unique values.
  final_vars <- vars_to_factor[sapply(vars_to_factor, function(v) {
    length(unique(na.omit(df[[v]]))) <= max_levels
  })]
  
  if (length(final_vars) == 0) {
    message("No variables met the criteria for factor conversion.")
    return(df)
  }
  
  # 5. Perform the conversion
  message(paste("Converting the following to factors:", paste(final_vars, collapse = ", ")))
  
  df[final_vars] <- lapply(df[final_vars], as.factor)
  
  return(df)
}
