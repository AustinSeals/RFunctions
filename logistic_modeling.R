

#' Perform Univariate Logistic Regression across multiple predictors
#'
#' @param df A dataframe containing the data.
#' @param dep_var A string representing the dependent variable name (must be binary 0/1).
#' @param predictors A character vector of predictor variable names.
#' @return A tidy dataframe with Beta, OR, CIs, and P-values.
run_univariate_models <- function(df, dep_var, predictors) {
  
  # Initialize an empty list to store results
  results_list <- list()
  
  for (pred in predictors) {
    
    # Create the formula: dep_var ~ predictor
    formula_str <- paste(dep_var, "~", pred)
    
    # Fit the logistic regression model
    model <- tryCatch({
      glm(as.formula(formula_str), data = df, family = binomial(link = "logit"))
    }, error = function(e) {
      warning(paste("Model failed for predictor:", pred, "\nReason:", e$message))
      return(NULL)
    })
    
    # Skip to next iteration if model failed
    if (is.null(model)) next
    
    # Extract coefficient summary
    coef_summary <- summary(model)$coefficients
    
    # Calculate Confidence Intervals (using Wald CIs for speed/stability in loops)
    ci <- suppressMessages(confint.default(model))
    
    # Drop the intercept row (usually row 1)
    coef_summary <- coef_summary[-1, , drop = FALSE]
    ci <- ci[-1, , drop = FALSE]
    
    # If a predictor has only NA/dropped levels, skip
    if (nrow(coef_summary) == 0) next
    
    # Extract values
    betas <- coef_summary[, "Estimate"]
    p_values <- coef_summary[, "Pr(>|z|)"]
    odds_ratios <- exp(betas)
    lower_ci <- exp(ci[, 1])
    upper_ci <- exp(ci[, 2])
    
    # Clean up variable labels (handling categorical factors)
    raw_terms <- rownames(coef_summary)
    
    # R pastes the variable name and the factor level together (e.g., "ColorRed").
    # We strip the variable name to isolate the level.
    levels <- sub(paste0("^", pred), "", raw_terms)
    
    # If the level is empty, it means the variable was continuous.
    # Otherwise, it's categorical and we format it nicely.
    display_names <- ifelse(levels == "", 
                            pred, 
                            paste0(pred ,"(", levels, ")"))
    
    # Combine into a dataframe for this specific predictor
    pred_df <- data.frame(
      Variable = display_names,
      Predictor_Base = pred,
      Level = ifelse(levels == "", "Continuous", levels),
      Beta = betas,
      Odds_Ratio = odds_ratios,
      Lower_CI_95 = lower_ci,
      Upper_CI_95 = upper_ci,
      P_Value = p_values,
      stringsAsFactors = FALSE
    )
    
    # Append to our list
    results_list[[pred]] <- pred_df
  }
  
  # Combine all predictor dataframes into one final dataframe
  final_df <- do.call(rbind, results_list)
  rownames(final_df) <- NULL # Clean up rownames
  
  return(final_df)
}



#' Perform Univariate Logistic Regression WITH Reference Level Labels
#'
#' @param df A dataframe containing the data.
#' @param dep_var A string representing the dependent variable name.
#' @param predictors A character vector of predictor variable names.
#' @return A tidy dataframe including Beta, OR, CIs, and Reference Level info.
run_univariate_models2 <- function(df, dep_var, predictors) {
  
  results_list <- list()
  
  for (pred in predictors) {
    # Create the formula
    formula_str <- paste(dep_var, "~", pred)
    
    # Fit the model
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
    
    # Remove Intercept
    coef_summary <- coef_summary[-1, , drop = FALSE]
    ci <- ci[-1, , drop = FALSE]
    
    if (nrow(coef_summary) == 0) next
    
    # Logic for Reference Level
    # If the variable is a factor or character, find the first level
    ref_level <- "N/A (Continuous)"
    if (is.factor(df[[pred]]) || is.character(df[[pred]])) {
      ref_level <- levels(as.factor(df[[pred]]))[1]
    }
    
    # Get raw term names (e.g., "cyl6", "cyl8")
    raw_terms <- rownames(coef_summary)
    
    # Clean labels: "cyl (Level: 6 vs Ref: 4)"
    levels_detected <- sub(paste0("^", pred), "", raw_terms)
    
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
      Beta = coef_summary[, "Estimate"],
      Odds_Ratio = exp(coef_summary[, "Estimate"]),
      Lower_CI_95 = exp(ci[, 1]),
      Upper_CI_95 = exp(ci[, 2]),
      P_Value = coef_summary[, "Pr(>|z|)"],
      stringsAsFactors = FALSE
    )
    
    results_list[[pred]] <- pred_df
  }
  
  final_df <- do.call(rbind, results_list)
  rownames(final_df) <- NULL
  return(final_df)
}





#' Generate a Forest Plot from Univariate Logistic Regression Results
#'
#' @param results_df The dataframe output from run_univariate_models().
#' @param title A string for the plot title.
#' @return A ggplot object.
plot_forest <- function(results_df, title = "Univariate Odds Ratios (95% CI)") {
  
  # Check if ggplot2 is installed
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("The 'ggplot2' package is required to generate this plot. Please install it.")
  }
  
  # Lock in the factor levels in reverse order so the first variable plots at the top
  results_df$Variable <- factor(results_df$Variable, levels = rev(unique(results_df$Variable)))
  
  # Create the plot
  p <- ggplot2::ggplot(results_df, ggplot2::aes(x = Odds_Ratio, y = Variable)) +
    # Add the line of no effect (Null Hypothesis: OR = 1)
    ggplot2::geom_vline(xintercept = 1, linetype = "dashed", color = "darkred", linewidth = 0.8, alpha = 0.7) +
    # Add horizontal error bars for the Confidence Intervals
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = Lower_CI_95, xmax = Upper_CI_95, color = Predictor_Base), 
                            height = 0.25, linewidth = 0.8, show.legend = FALSE) +
    # Add the point estimate for the Odds Ratio
    ggplot2::geom_point(ggplot2::aes(color = Predictor_Base), size = 3, show.legend = FALSE) +
    # Transform the X-axis to a log scale (standard practice for Odds Ratios)
    ggplot2::scale_x_log10() +
    # Labeling and theming
    ggplot2::labs(
      title = title,
      x = "Odds Ratio (Log Scale)",
      y = NULL
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 14, hjust = 0.5),
      axis.text.y = ggplot2::element_text(size = 11, face = "bold", color = "black"),
      axis.title.x = ggplot2::element_text(size = 12, margin = ggplot2::margin(t = 10)),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(color = "grey90", linetype = "dotted")
    )
  
  return(p)
}
