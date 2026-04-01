library(survival)
#' Run Adjusted Univariate Cox Models with PH Assumption Check
#'
#' @param df A dataframe containing the data.
#' @param time_var A string for the time-to-event variable.
#' @param event_var A string for the status/event variable (1=event, 0=censored).
#' @param predictors A character vector of the main predictors of interest.
#' @param covariates A character vector of variables to adjust for.
#' @return A tidy dataframe including Hazard Ratios and PH Assumption status.
run_cox_adjusted_models <- function(df, time_var, event_var, predictors, covariates = NULL) {
  
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("The 'survival' package is required.")
  }
  
  results_list <- list()
  covar_string <- if (!is.null(covariates)) paste(" +", paste(covariates, collapse = " + ")) else ""
  
  for (pred in predictors) {
    
    formula_str <- paste0("survival::Surv(", time_var, ", ", event_var, ") ~ ", pred, covar_string)
    
    model <- tryCatch({
      survival::coxph(as.formula(formula_str), data = df)
    }, error = function(e) {
      warning(paste("Cox model failed for predictor:", pred))
      return(NULL)
    })
    
    if (is.null(model)) next
    
    # --- PH Assumption Test (Schoenfeld residuals) ---
    # cox.zph returns a test for each term + a global test.
    ph_test <- survival::cox.zph(model)
    ph_p_values <- ph_test$table[, "p"]
    
    # Extract model stats
    summ <- summary(model)
    coef_table <- summ$coefficients
    ci_table <- summ$conf.int
    
    # Isolate rows belonging to the main 'pred'
    all_rows <- rownames(coef_table)
    pred_rows_idx <- which(grepl(paste0("^", pred), all_rows))
    
    if (length(pred_rows_idx) == 0) next
    
    # Reference level logic
    ref_level <- "N/A (Continuous)"
    if (is.factor(df[[pred]]) || is.character(df[[pred]])) {
      ref_level <- levels(as.factor(df[[pred]]))[1]
    }
    
    # Map PH test results back to the specific levels of our predictor
    # Note: ph_test$table rows match the model coefficients order
    ph_results <- ifelse(ph_p_values[pred_rows_idx] >= 0.05, "Yes", "No")
    
    # Clean labels
    levels_detected <- sub(paste0("^", pred), "", all_rows[pred_rows_idx])
    display_names <- ifelse(levels_detected == "", pred, paste0(pred, " (", levels_detected, " vs ", ref_level, ")"))
    
    # Build dataframe
    pred_df <- data.frame(
      Variable = display_names,
      Predictor_Base = pred,
      Comparison_Level = ifelse(levels_detected == "", "Continuous", levels_detected),
      Reference_Level = ref_level,
      Hazard_Ratio = ci_table[pred_rows_idx, "exp(coef)"],
      Lower_CI_95 = ci_table[pred_rows_idx, "lower .95"],
      Upper_CI_95 = ci_table[pred_rows_idx, "upper .95"],
      P_Value = coef_table[pred_rows_idx, "Pr(>|z|)"],
      PH_Assumption_Met = ph_results,  # NEW COLUMN
      PH_Test_P = ph_p_values[pred_rows_idx], # Added for transparency
      Adjusted_For = paste(covariates, collapse = ", "),
      stringsAsFactors = FALSE
    )
    
    results_list[[pred]] <- pred_df
  }
  
  final_df <- do.call(rbind, results_list)
  rownames(final_df) <- NULL
  return(final_df)
}



#' Perform Multivariable Cox Regression with PH Assumption Check
#'
#' @param df A dataframe containing the data.
#' @param time_var A string for the time-to-event variable.
#' @param event_var A string for the status/event variable (1=event, 0=censored).
#' @param predictors A character vector of all predictor variable names to include.
#' @return A tidy dataframe with adjusted HRs and PH Assumption status.
run_multivariable_cox_ph <- function(df, time_var, event_var, predictors) {
  
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("The 'survival' package is required.")
  }
  
  # 1. Build and Fit the Multivariable Model
  formula_str <- paste0("survival::Surv(", time_var, ", ", event_var, ") ~ ", 
                        paste(predictors, collapse = " + "))
  
  model <- tryCatch({
    survival::coxph(as.formula(formula_str), data = df)
  }, error = function(e) {
    stop(paste("The multivariable Cox model failed.\nError:", e$message))
  })
  
  # 2. Run the PH Assumption Test (Schoenfeld residuals)
  # This calculates p-values for every term in the model
  ph_test <- survival::cox.zph(model)
  ph_table <- ph_test$table
  
  # 3. Extract standard model statistics
  summ <- summary(model)
  coef_table <- summ$coefficients
  ci_table <- summ$conf.int
  raw_terms <- rownames(coef_table)
  
  results_list <- list()
  
  for (i in seq_along(raw_terms)) {
    term_name <- raw_terms[i]
    
    # Identify base predictor
    base_pred <- predictors[sapply(predictors, function(p) grepl(paste0("^", p), term_name))][1]
    
    # Reference level logic
    ref_level <- "N/A (Continuous)"
    if (is.factor(df[[base_pred]]) || is.character(df[[base_pred]])) {
      ref_level <- levels(as.factor(df[[base_pred]]))[1]
    }
    
    # Clean labels
    level_detected <- sub(paste0("^", base_pred), "", term_name)
    display_name <- ifelse(
      level_detected == "", 
      base_pred, 
      paste0(base_pred, " (", level_detected, " vs ", ref_level, ")")
    )
    
    # Extract PH p-value for THIS specific term
    # Note: the row names in ph_table match the model coefficients
    ph_p <- ph_table[term_name, "p"]
    
    results_list[[i]] <- data.frame(
      Variable = display_name,
      Predictor_Base = base_pred,
      Comparison_Level = ifelse(level_detected == "", "Continuous", level_detected),
      Reference_Level = ref_level,
      Beta = coef_table[i, "coef"],
      Hazard_Ratio = ci_table[i, "exp(coef)"],
      Lower_CI_95 = ci_table[i, "lower .95"],
      Upper_CI_95 = ci_table[i, "upper .95"],
      P_Value = coef_table[i, "Pr(>|z|)"],
      PH_Assumption_Met = ifelse(ph_p >= 0.05, "Yes", "No"),
      PH_Test_P = ph_p,
      stringsAsFactors = FALSE
    )
  }
  
  final_df <- do.call(rbind, results_list)
  
  # Add the Global Model PH Test as an attribute to the dataframe
  # This is useful for checking if the model as a whole is valid
  global_p <- ph_table["GLOBAL", "p"]
  attr(final_df, "Global_PH_P") <- global_p
  
  message(paste("Global PH Assumption Check P-value:", round(global_p, 4)))
  
  rownames(final_df) <- NULL
  return(final_df)
}








#' Generate a Forest Plot for Hazard Ratios
#'
#' @param results_df The dataframe output from run_cox_adjusted_models or run_multivariable_cox.
#' @param title A string for the plot title.
#' @return A ggplot object.
plot_cox_forest <- function(results_df, title = "Hazard Ratios (95% CI)") {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("The 'ggplot2' package is required.")
  }
  
  # Ensure the data is from a Cox model
  if (!"Hazard_Ratio" %in% colnames(results_df)) {
    stop("This function requires a 'Hazard_Ratio' column. Use run_multivariable_cox first.")
  }
  
  # 1. Reverse the variable order for top-to-bottom reading
  results_df$Variable <- factor(results_df$Variable, levels = rev(unique(results_df$Variable)))
  
  # 2. Build the Plot
  p <- ggplot2::ggplot(results_df, ggplot2::aes(x = Hazard_Ratio, y = Variable)) +
    # Vertical reference line at HR = 1 (Null Hypothesis)
    ggplot2::geom_vline(xintercept = 1, linetype = "dashed", color = "darkred", linewidth = 0.8) +
    # Confidence Interval bars
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = Lower_CI_95, xmax = Upper_CI_95, color = Predictor_Base), 
                            height = 0.2, linewidth = 1, show.legend = FALSE) +
    # Point estimate (the HR)
    ggplot2::geom_point(ggplot2::aes(color = Predictor_Base), size = 4, show.legend = FALSE) +
    # Logarithmic scale for Hazard Ratios
    ggplot2::scale_x_log10(breaks = c(0.1, 0.25, 0.5, 1, 2, 4, 10)) +
    # Labs and Annotation
    ggplot2::labs(
      title = title,
      subtitle = NULL,
      x = "Hazard Ratio (Log Scale)",
      y = NULL
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 16),
      plot.subtitle = ggplot2::element_text(size = 10, color = "grey40", face = "italic"),
      axis.text.y = ggplot2::element_text(size = 11, face = "bold", color = "black"),
      axis.title.x = ggplot2::element_text(size = 12, margin = ggplot2::margin(t = 15)),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(color = "grey95")
    )
  
  return(p)
}


