

#' Check Linearity of Continuous Predictors using Martingale Residuals
#'
#' @param model A coxph or glm object.
#' @param var_name The string name of the continuous variable.
plot_linearity_check <- function(model, var_name, df) {
  # For Cox models, Martingale residuals vs the predictor should be linear
  res <- residuals(model, type = "martingale")
  
  ggplot2::ggplot(df, ggplot2::aes_string(x = var_name, y = "res")) +
    ggplot2::geom_point(alpha = 0.3) +
    ggplot2::geom_smooth(method = "loess", color = "blue") +
    ggplot2::labs(title = paste("Linearity Check:", var_name),
                  y = "Martingale Residuals",
                  subtitle = "The blue line should be roughly straight/horizontal") +
    ggplot2::theme_minimal()
}