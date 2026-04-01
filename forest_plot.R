


#' Column-Agnostic Forest Plot Function
#'
#' @param df The results dataframe.
#' @param label_col The column containing the variable names/labels (y-axis).
#' @param estimate_col The column containing the point estimates (x-axis).
#' @param low_col The column containing the lower confidence limit.
#' @param high_col The column containing the upper confidence limit.
#' @param group_col The column used to color-code variables (e.g., Predictor_Base).
#' @param x_label The label for the x-axis.
#' @param title The plot title.
#' @param log_scale Logical; if TRUE, uses a log10 scale (standard for OR/HR).
#' @return A ggplot object.
plot_forest_agnostic <- function(df, 
                                 label_col = "Variable", 
                                 estimate_col = "Hazard_Ratio", 
                                 low_col = "Lower_CI_95", 
                                 high_col = "Upper_CI_95", 
                                 group_col = "Predictor_Base",
                                 x_label = "Estimate (95% CI)",
                                 title = "Forest Plot",
                                 log_scale = TRUE) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Please install 'ggplot2'.")
  }
  
  # Ensure the y-axis labels maintain the order they appear in the dataframe
  df[[label_col]] <- factor(df[[label_col]], levels = rev(unique(df[[label_col]])))
  
  # Build the base plot using .data[[]] for tidy evaluation of string arguments
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[estimate_col]], y = .data[[label_col]])) +
    # Reference line at 1.0 (or 0 for linear models)
    ggplot2::geom_vline(xintercept = ifelse(log_scale, 1, 0), 
                        linetype = "dashed", color = "darkred", linewidth = 0.8) +
    # Confidence Intervals
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = .data[[low_col]], 
                                         xmax = .data[[high_col]], 
                                         color = .data[[group_col]]), 
                            height = 0.2, linewidth = 1, show.legend = FALSE) +
    # Point Estimates
    ggplot2::geom_point(ggplot2::aes(color = .data[[group_col]]), 
                        size = 3.5, show.legend = FALSE) +
    ggplot2::labs(title = title, x = x_label, y = NULL) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 14),
      axis.text.y = ggplot2::element_text(size = 10, face = "bold", color = "black"),
      panel.grid.minor = ggplot2::element_blank()
    )
  
  # Apply log scale if requested
  if (log_scale) {
    p <- p + ggplot2::scale_x_log10()
  }
  
  return(p)
}