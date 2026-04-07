
#' Run Univariate Models adjusted for a fixed set of covariates
#'
#' @param df A dataframe containing the data.
#' @param dep_var A string for the dependent variable.
#' @param predictors A character vector of the main predictors of interest.
#' @param covariates A character vector of variables to adjust for in every model.
#' @return A tidy dataframe of the adjusted estimates for the main predictors.
run_univariate_adjusted_models <- function(df, dep_var, predictors, covariates = NULL) {
  
  results_list <- list()
  
  # Create the string for covariates (e.g., "+ age + sex")
  covar_string <- ""
  if (!is.null(covariates) && length(covariates) > 0) {
    covar_string <- paste(" +", paste(covariates, collapse = " + "))
  }
  
  for (pred in predictors) {
    
    # Construct formula: dep ~ predictor + adj1 + adj2...
    formula_str <- paste(dep_var, "~", pred, covar_string)
    
    model <- tryCatch({
      glm(as.formula(formula_str), data = df, family = binomial(link = "logit"))
    }, error = function(e) {
      warning(paste("Model failed for predictor:", pred))
      return(NULL)
    })
    
    if (is.null(model)) next
    
    # Extract stats
    coef_summary <- summary(model)$coefficients
    ci <- suppressMessages(confint.default(model))
    
    # --- The Selection Logic ---
    # We only want the rows corresponding to the 'pred', not the covariates or intercept.
    # We find rows that start with the predictor name.
    all_rows <- rownames(coef_summary)
    pred_rows_idx <- which(grepl(paste0("^", pred), all_rows))
    
    if (length(pred_rows_idx) == 0) next
    
    # Extract only the predictor-specific rows
    coef_subset <- coef_summary[pred_rows_idx, , drop = FALSE]
    ci_subset <- ci[pred_rows_idx, , drop = FALSE]
    
    # Logic for Reference Level
    ref_level <- "N/A (Continuous)"
    if (is.factor(df[[pred]]) || is.character(df[[pred]])) {
      ref_level <- levels(as.factor(df[[pred]]))[1]
    }
    
    # Clean labels
    levels_detected <- sub(paste0("^", pred), "", rownames(coef_subset))
    display_names <- ifelse(
      levels_detected == "", 
      pred, 
      paste0(pred, " (", levels_detected, " vs ", ref_level, ")")
    )
    
    # Build dataframe for this predictor
    pred_df <- data.frame(
      Variable = display_names,
      Predictor_Base = pred,
      Comparison_Level = ifelse(levels_detected == "", "Continuous", levels_detected),
      Reference_Level = ref_level,
      Beta = coef_subset[, "Estimate"],
      Odds_Ratio = exp(coef_subset[, "Estimate"]),
      Lower_CI_95 = exp(ci_subset[, 1]),
      Upper_CI_95 = exp(ci_subset[, 2]),
      P_Value = coef_subset[, "Pr(>|z|)"],
      Adjusted_For = paste(covariates, collapse = ", "),
      stringsAsFactors = FALSE
    )
    
    results_list[[pred]] <- pred_df
  }
  
  final_df <- do.call(rbind, results_list)
  rownames(final_df) <- NULL
  return(final_df)
}









#' Perform Multivariable Logistic Regression
#'
#' @param df A dataframe containing the data.
#' @param dep_var A string representing the dependent variable name.
#' @param predictors A character vector of all predictor variable names to include in the model.
#' @return A tidy dataframe with adjusted Beta, OR, CIs, and Reference Level info.
run_multivariable_model <- function(df, dep_var, predictors) {
  
  # Construct the full formula: dep_var ~ pred1 + pred2 + ...
  formula_str <- paste(dep_var, "~", paste(predictors, collapse = " + "))
  
  # Fit the single multivariable model
  model <- tryCatch({
    glm(as.formula(formula_str), data = df, family = binomial(link = "logit"))
  }, error = function(e) {
    stop(paste("The multivariable model failed to converge. Check for perfect separation or collinearity.\nError:", e$message))
  })
  
  # Extract coefficient summary and Wald CIs
  coef_summary <- summary(model)$coefficients
  ci <- suppressMessages(confint.default(model))
  
  # Remove the Intercept
  coef_summary <- coef_summary[-1, , drop = FALSE]
  ci <- ci[-1, , drop = FALSE]
  
  # Create a mapping of which predictor produced which row in the output
  # (Important because categorical variables create multiple rows)
  raw_terms <- rownames(coef_summary)
  
  results_list <- list()
  
  for (i in seq_along(raw_terms)) {
    term_name <- raw_terms[i]
    
    # Identify which base predictor this term belongs to
    # We find the predictor name that is a prefix of the term_name
    base_pred <- predictors[sapply(predictors, function(p) grepl(paste0("^", p), term_name))][1]
    
    # Logic for Reference Level
    ref_level <- "N/A (Continuous)"
    if (is.factor(df[[base_pred]]) || is.character(df[[base_pred]])) {
      ref_level <- levels(as.factor(df[[base_pred]]))[1]
    }
    
    # Clean labels: "cyl (6 vs 4)"
    level_detected <- sub(paste0("^", base_pred), "", term_name)
    
    display_name <- ifelse(
      level_detected == "", 
      base_pred, 
      paste0(base_pred, " (", level_detected, " vs ", ref_level, ")")
    )
    
    results_list[[i]] <- data.frame(
      Variable = display_name,
      Predictor_Base = base_pred,
      Comparison_Level = ifelse(level_detected == "", "Continuous", level_detected),
      Reference_Level = ref_level,
      Beta = coef_summary[i, "Estimate"],
      Odds_Ratio = exp(coef_summary[i, "Estimate"]),
      Lower_CI_95 = exp(ci[i, 1]),
      Upper_CI_95 = exp(ci[i, 2]),
      P_Value = coef_summary[i, "Pr(>|z|)"],
      stringsAsFactors = FALSE
    )
  }
  
  final_df <- do.call(rbind, results_list)
  return(final_df)
}








