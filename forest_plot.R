#' Final Agnostic Forest Plot (Legend Suppressed)
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
                              log_scale = TRUE) {
  
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

  # Logic for Conditional Coloring (Legend specifically suppressed in geoms)
  if (!is.null(color_var) && color_var %in% colnames(df)) {
    p <- p + 
      ggplot2::geom_errorbarh(ggplot2::aes(xmin = .data[[low_col]], 
                                          xmax = .data[[high_col]], 
                                          color = .data[[color_var]]), 
                              height = 0.2, linewidth = 1, show.legend = FALSE) +
      ggplot2::geom_point(ggplot2::aes(color = .data[[color_var]]), 
                          size = 3.5, show.legend = FALSE)
  } else {
    p <- p + 
      ggplot2::geom_errorbarh(ggplot2::aes(xmin = .data[[low_col]], 
                                          xmax = .data[[high_col]]), 
                              color = "midnightblue", height = 0.2, linewidth = 1) +
      ggplot2::geom_point(color = "midnightblue", size = 3.5)
  }

  # Add labels and global theme settings
  p <- p + 
    ggplot2::labs(title = title, x = x_title, y = y_title) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 14),
      axis.text.y = ggplot2::element_text(size = 10, face = "bold", color = "black"),
      panel.grid.minor = ggplot2::element_blank(),
      # Absolute suppression of the legend
      legend.position = "none" 
    )

  if (log_scale) {
    p <- p + ggplot2::scale_x_log10(breaks = c(0.1, 0.5, 1, 2, 5, 10))
  }

  return(p)
}








#' Sorted Agnostic Forest Plot (Categorized, Color-Coded, and Wrapped)
#'
#' @param df The results dataframe.
#' @param label_col The column for y-axis labels.
#' @param wrap_width The character width at which to wrap y-axis text. Default is 30.
#' @param estimate_col The column for x-axis point estimates.
#' @param low_col The column for the lower CI.
#' @param high_col The column for the upper CI.
#' @param x_title Custom string for the X-axis title.
#' @param y_title Custom string for the Y-axis title.
#' @param title Custom string for the plot title.
#' @param log_scale Logical; TRUE for OR/HR (Null = 1), FALSE for linear (Null = 0).
#' @return A ggplot object.
plot_forest_categorized <- function(df, 
                              label_col = "Variable",
                              wrap_width = 100,
                              estimate_col = "Hazard_Ratio",
                              low_col = "Lower_CI_95",
                              high_col = "Upper_CI_95",
                              x_title = "Estimate (95% CI)",
                              y_title = NULL,
                              title = "Forest Plot",
                              log_scale = TRUE) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Please install 'ggplot2'.")
  if (!requireNamespace("stringr", quietly = TRUE)) stop("Please install 'stringr'.")
  
  # 1. Determine Null Value
  null_val <- ifelse(log_scale, 1, 0)
  
  # 2. Categorize Effect based on Confidence Intervals
  df$Effect_Type <- mapply(function(low, high) {
    if (is.na(low) || is.na(high)) return("Non-Predictor")
    if (low > null_val) return("Positive Predictor")
    if (high < null_val) return("Negative Predictor")
    return("Non-Predictor")
  }, df[[low_col]], df[[high_col]])
  
  # 3. Create a Sorting Rank: Positive (3) -> Negative (2) -> Non (1)
  df$Sort_Rank <- ifelse(df$Effect_Type == "Positive Predictor", 3,
                         ifelse(df$Effect_Type == "Negative Predictor", 2, 1))
  
  # 4. Sort and Wrap Labels
  # Sort first so we can extract the correct color order
  df <- df[order(df$Sort_Rank, df[[estimate_col]]), ]
  
  # Wrap the text labels
  df[[label_col]] <- stringr::str_wrap(df[[label_col]], width = wrap_width)
  
  # Re-set factor levels based on the newly wrapped and sorted labels
  df[[label_col]] <- factor(df[[label_col]], levels = unique(df[[label_col]]))
  
  # 5. Define Internal Color Palette
  impact_colors <- c(
    "Positive Predictor" = "#2E8B57", # SeaGreen
    "Negative Predictor" = "#CD5C5C", # IndianRed
    "Non-Predictor"      = "#708090"  # SlateGray
  )
  
  # 6. Extract colors for the Y-axis
  y_axis_colors <- impact_colors[df$Effect_Type]
  
  # 7. Build the Plot
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[estimate_col]], 
                                        y = .data[[label_col]], 
                                        color = Effect_Type)) +
    ggplot2::geom_vline(xintercept = null_val, linetype = "dashed", 
                        color = "darkred", linewidth = 0.8) +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = .data[[low_col]], 
                                         xmax = .data[[high_col]]), 
                            height = 0.2, linewidth = 1, show.legend = FALSE) +
    ggplot2::geom_point(size = 3.5, show.legend = FALSE) +
    ggplot2::scale_color_manual(values = impact_colors) +
    ggplot2::labs(
      title = title, 
      subtitle = "Green: Positive | Red: Negative | Gray: Non-Significant",
      x = x_title, 
      y = y_title
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 14),
      plot.subtitle = ggplot2::element_text(size = 9, color = "grey40", face = "italic"),
      # Applying the color vector to the axis text
      axis.text.y = ggplot2::element_text(size = 10, face = "bold", color = y_axis_colors, lineheight = 0.9),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "none"
    )
  
  if (log_scale) {
    p <- p + ggplot2::scale_x_log10(breaks = c(0.1, 0.5, 1, 2, 5, 10))
  }
  
  return(p)
}


