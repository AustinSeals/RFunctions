
#' Fully Customizable Forest Plot Function
#'
#' @param df The results dataframe.
#' @param label_col The column for the y-axis labels.
#' @param estimate_col The column for the x-axis point estimates.
#' @param low_col The column for the lower CI.
#' @param high_col The column for the upper CI.
#' @param group_col The column used for color-coding.
#' @param x_title Custom string for the X-axis title.
#' @param y_title Custom string for the Y-axis title (defaults to NULL).
#' @param title Custom string for the plot title.
#' @param log_scale Logical; set to TRUE for OR/HR, FALSE for Beta coefficients.
#' @return A ggplot object.
plot_forest_custom <- function(df, 
                               label_col = "Variable", 
                               estimate_col = "Hazard_Ratio", 
                               low_col = "Lower_CI_95", 
                               high_col = "Upper_CI_95", 
                               group_col = "Predictor_Base",
                               x_title = "Estimate (95% CI)",
                               y_title = NULL,
                               title = "Forest Plot",
                               log_scale = TRUE) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Please install 'ggplot2'.")
  }

  # Maintain original data order on the Y-axis
  df[[label_col]] <- factor(df[[label_col]], levels = rev(unique(df[[label_col]])))

  # Map variables to the plot
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[estimate_col]], y = .data[[label_col]])) +
    # Vertical reference line (1 for ratios, 0 for linear differences)
    ggplot2::geom_vline(xintercept = ifelse(log_scale, 1, 0), 
                        linetype = "dashed", color = "darkred", linewidth = 0.8) +
    # Horizontal Error Bars
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = .data[[low_col]], 
                                        xmax = .data[[high_col]], 
                                        color = .data[[group_col]]), 
                            height = 0.2, linewidth = 1, show.legend = FALSE) +
    # Point Estimates
    ggplot2::geom_point(ggplot2::aes(color = .data[[group_col]]), 
                        size = 3.5, show.legend = FALSE) +
    # Use the custom labels provided by the user
    ggplot2::labs(
      title = title, 
      x = x_title, 
      y = y_title
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 14, margin = ggplot2::margin(b = 10)),
      axis.text.y = ggplot2::element_text(size = 10, face = "bold", color = "black"),
      axis.title.x = ggplot2::element_text(size = 11, margin = ggplot2::margin(t = 10)),
      axis.title.y = ggplot2::element_text(size = 11, margin = ggplot2::margin(r = 10)),
      panel.grid.minor = ggplot2::element_blank()
    )

  if (log_scale) {
    # Adding common log breaks for better readability
    p <- p + ggplot2::scale_x_log10(breaks = c(0.1, 0.2, 0.5, 1, 2, 5, 10))
  }

  return(p)
}
