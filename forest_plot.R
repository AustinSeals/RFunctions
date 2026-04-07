#' Final Agnostic Forest Plot with Optional Coloring
#'
#' @param df The results dataframe.
#' @param color_var The column name to use for coloring (e.g., "Predictor_Base"). 
#'                  Set to NULL for a single-color plot.
#' @param label_col The column for y-axis labels.
#' @param estimate_col The column for x-axis point estimates.
#' @param low_col The column for the lower CI.
#' @param high_col The column for the upper CI.
#' @param x_title Custom string for the X-axis title.
#' @param y_title Custom string for the Y-axis title.
#' @param title Custom string for the plot title.
#' @param log_scale Logical; TRUE for OR/HR, FALSE for linear models.
#' @param x_breaks breaks for the log scale x-axis . Defaault c(0.1, 0.5, 1, 2, 5, 10) 
#' @return A ggplot object.
plot_forest_final <- function(df, 
                              color_var = NULL,
                              label_col = "Variable",
                              estimate_col = "Hazard_Ratio",
                              low_col = "Lower_CI_95",
                              high_col = "Upper_CI_95",
                              x_title = "Estimate (95% CI)",
                              y_title = NULL,
                              title = "Forest Plot",
                              log_scale = TRUE
                              x_breaks  =  c(0.1, 0.5, 1, 2, 5, 10) 
                             ) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Please install 'ggplot2'.")
  }

  # Ensure Y-axis order matches the dataframe
  df[[label_col]] <- factor(df[[label_col]], levels = rev(unique(df[[label_col]])))

  # Initialize the base plot
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[estimate_col]], y = .data[[label_col]])) +
    # Vertical reference line
    ggplot2::geom_vline(xintercept = ifelse(log_scale, 1, 0), 
                        linetype = "dashed", color = "darkred", linewidth = 0.8)

  # Logic for Conditional Coloring
  if (!is.null(color_var) && color_var %in% colnames(df)) {
    # If color_var is provided and exists, map it to color
    p <- p + 
      ggplot2::geom_errorbarh(ggplot2::aes(xmin = .data[[low_col]], 
                                          xmax = .data[[high_col]], 
                                          color = .data[[color_var]]), 
                              height = 0.2, linewidth = 1) +
      ggplot2::geom_point(ggplot2::aes(color = .data[[color_var]]), size = 3.5) +
      ggplot2::labs(color = color_var)
  } else {
    # If no color_var is specified, use a default single color
    p <- p + 
      ggplot2::geom_errorbarh(ggplot2::aes(xmin = .data[[low_col]], 
                                          xmax = .data[[high_col]]), 
                              color = "midnightblue", height = 0.2, linewidth = 1) +
      ggplot2::geom_point(color = "midnightblue", size = 3.5)
  }

  # Add labels and theme
  p <- p + 
    ggplot2::labs(title = title, x = x_title, y = y_title) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 14),
      axis.text.y = ggplot2::element_text(size = 10, face = "bold", color = "black"),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = if(is.null(color_var)) "none" else "right"
    )

  if (log_scale) {
    p <- p + ggplot2::scale_x_log10(breaks = x_breaks)
  }

  return(p)
}
