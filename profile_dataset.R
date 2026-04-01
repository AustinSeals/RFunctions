
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