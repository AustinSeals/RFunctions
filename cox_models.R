library(survival)

#' Run Univariate Cox Models adjusted for a fixed set of covariates
#'
#' @param df A dataframe containing the data.
#' @param time_var A string for the time-to-event variable.
#' @param event_var A string for the status/event variable (1=event, 0=censored).
#' @param predictors A character vector of the main predictors of interest.
#' @param covariates A character vector of variables to adjust for in every model.
#' @return A tidy dataframe of the adjusted Hazard Ratios for the main predictors.
run_cox_adjusted_models <- function(df, time_var, event_var, predictors, covariates = NULL) {
  
  # Ensure the survival package is available
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("The 'survival' package is required. Please install it.")
  }
  
  results_list <- list()
  
  # Create the string for covariates (e.g., "+ age + sex")
  covar_string <- ""
  if (!is.null(covariates) && length(covariates) > 0) {
    covar_string <- paste(" +", paste(covariates, collapse = " + "))
  }
  
  for (pred in predictors) {
    
    # Construct Cox formula: Surv(time, event) ~ predictor + adj1 + adj2...
    formula_str <- paste0("survival::Surv(", time_var, ", ", event_var, ") ~ ", pred, covar_string)
    
    model <- tryCatch({
      survival::coxph(as.formula(formula_str), data = df)
    }, error = function(e) {
      warning(paste("Cox model failed for predictor:", pred))
      return(NULL)
    })
    
    if (is.null(model)) next
    
    # Extract coefficients and Confidence Intervals
    # Note: Cox models in R use summary(model)$conf.int for HR and CIs directly
    summ <- summary(model)
    coef_table <- summ$coefficients  # Contains Beta, HR, se, z, p
    ci_table <- summ$conf.int        # Contains HR, exp(-beta), lower .95, upper .95
    
    # --- The Selection Logic ---
    # Isolate only the rows belonging to the main 'pred'
    all_rows <- rownames(coef_table)
    pred_rows_idx <- which(grepl(paste0("^", pred), all_rows))
    
    if (length(pred_rows_idx) == 0) next
    
    # Logic for Reference Level
    ref_level <- "N/A (Continuous)"
    if (is.factor(df[[pred]]) || is.character(df[[pred]])) {
      ref_level <- levels(as.factor(df[[pred]]))[1]
    }
    
    # Clean labels
    levels_detected <- sub(paste0("^", pred), "", all_rows[pred_rows_idx])
    display_names <- ifelse(
      levels_detected == "", 
      pred, 
      paste0(pred, " (", levels_detected, " vs ", ref_level, ")")
    )
    
    # Build dataframe
    pred_df <- data.frame(
      Variable = display_names,
      Predictor_Base = pred,
      Comparison_Level = ifelse(levels_detected == "", "Continuous", levels_detected),
      Reference_Level = ref_level,
      Beta = coef_table[pred_rows_idx, "coef"],
      Hazard_Ratio = ci_table[pred_rows_idx, "exp(coef)"],
      Lower_CI_95 = ci_table[pred_rows_idx, "lower .95"],
      Upper_CI_95 = ci_table[pred_rows_idx, "upper .95"],
      P_Value = coef_table[pred_rows_idx, "Pr(>|z|)"],
      Adjusted_For = paste(covariates, collapse = ", "),
      stringsAsFactors = FALSE
    )
    
    results_list[[pred]] <- pred_df
  }
  
  final_df <- do.call(rbind, results_list)
  rownames(final_df) <- NULL
  return(final_df)
}


#' Perform Multivariable Cox Proportional Hazards Regression
#'
#' @param df A dataframe containing the data.
#' @param time_var A string for the time-to-event variable.
#' @param event_var A string for the status/event variable (1=event, 0=censored).
#' @param predictors A character vector of all predictor variable names to include.
#' @return A tidy dataframe with adjusted Beta, HR, CIs, and Reference Level info.
run_multivariable_cox <- function(df, time_var, event_var, predictors) {
  
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("The 'survival' package is required. Please install it.")
  }
  
  # Construct the full formula: Surv(time, event) ~ pred1 + pred2 + ...
  formula_str <- paste0("survival::Surv(", time_var, ", ", event_var, ") ~ ", 
                        paste(predictors, collapse = " + "))
  
  # Fit the single multivariable Cox model
  model <- tryCatch({
    survival::coxph(as.formula(formula_str), data = df)
  }, error = function(e) {
    stop(paste("The multivariable Cox model failed. Check for collinearity or sparse events.\nError:", e$message))
  })
  
  # Extract summary statistics
  summ <- summary(model)
  coef_table <- summ$coefficients
  ci_table <- summ$conf.int
  
  # Get term names (e.g., "age", "sexFemale")
  raw_terms <- rownames(coef_table)
  results_list <- list()
  
  for (i in seq_along(raw_terms)) {
    term_name <- raw_terms[i]
    
    # Identify which base predictor this term belongs to
    base_pred <- predictors[sapply(predictors, function(p) grepl(paste0("^", p), term_name))][1]
    
    # Logic for Reference Level
    ref_level <- "N/A (Continuous)"
    if (is.factor(df[[base_pred]]) || is.character(df[[base_pred]])) {
      ref_level <- levels(as.factor(df[[base_pred]]))[1]
    }
    
    # Clean labels: "sex (Female vs Male)"
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
      Beta = coef_table[i, "coef"],
      Hazard_Ratio = ci_table[i, "exp(coef)"],
      Lower_CI_95 = ci_table[i, "lower .95"],
      Upper_CI_95 = ci_table[i, "upper .95"],
      P_Value = coef_table[i, "Pr(>|z|)"],
      stringsAsFactors = FALSE
    )
  }
  
  final_df <- do.call(rbind, results_list)
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
      subtitle = "<- Decreased Hazard (Protective) | Increased Hazard (Risk) ->",
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


