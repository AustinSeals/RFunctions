#' Create a finalfit-style summary table for repeated measures data
#'
#' @description
#' Generates a beautifully formatted summary table comparing explanatory variables 
#' across the levels of a dependent categorical variable, specifically designed for 
#' clustered or repeated-measures data. It mimics the output matrix of 
#' \code{finalfit::summary_factorlist()} while calculating p-values that correctly 
#' account for within-subject correlation using mixed-effects models.
#'
#' @details
#' The function automatically detects whether an explanatory variable is continuous 
#' or categorical and fits the appropriate mixed-effects model to generate a Type III Wald 
#' Chi-Square p-value via \code{car::Anova()}:
#' \itemize{
#'   \item \strong{Continuous (mean):} Fits a linear mixed model (\code{lmer}).
#'   \item \strong{Continuous (median):} Fits a rank-transformed linear mixed model 
#'         (\code{lmer(rank(x) ~ ...)}) to handle skewness and unbalanced/missing repeated measures robustly.
#'   \item \strong{Categorical:} Fits a binomial generalized linear mixed model (\code{glmer}).
#' }
#' 
#' @param data A data frame or tibble containing the dataset.
#' @param dependent Character vector of length 1. The name of the categorical grouping 
#'   variable (e.g., Treatment Arm).
#' @param explanatory Character vector. The names of the variables to be summarized 
#'   and tested across the dependent variable.
#' @param random_effect Character vector of length 1. The name of the clustering 
#'   variable (e.g., Patient ID, Subject ID) to be used as a random intercept.
#' @param cont_summary Character string. Either \code{"mean"} to summarize continuous 
#'   variables as Mean (SD), or \code{"median"} for Median [IQR]. Default is \code{"mean"}.
#' @param digits Integer. The number of decimal places to display for numeric 
#'   summaries and percentages. Default is \code{1}.
#' @param col_pct Logical. If \code{TRUE}, categorical summaries display column 
#'   percentages. If \code{FALSE}, displays row percentages. Default is \code{TRUE}.
#' @param total_col Logical. If \code{TRUE}, includes an overall total column 
#'   combining all groups. Default is \code{TRUE}.
#' @param add_missing Logical. If \code{TRUE}, appends a row showing the count of 
#'   \code{NA} (missing) values for each explanatory variable. Default is \code{TRUE}.
#'
#' @return A \code{tibble} containing the formatted summary statistics and p-values, 
#'   structured with columns for the variable label, level, dependent group levels, 
#'   total (optional), and the p-value.
#'
#' @importFrom lme4 lmer glmer
#' @importFrom car Anova
#' @importFrom dplyr filter group_by summarise mutate select bind_rows bind_cols if_else all_of
#' @importFrom tidyr pivot_wider
#' @importFrom rlang sym
#' 
#' @export
#'
#' @examples
#' \dontrun{
#' set.seed(42)
#' mock_data <- data.frame(
#'   PatientID = rep(1:100, each = 3),
#'   Timepoint = rep(c("Baseline", "Month 1", "Month 6"), 100),
#'   Group = rep(sample(c("Drug A", "Placebo"), 100, replace = TRUE), each = 3),
#'   HeartRate = rnorm(300, mean = 75, sd = 10),
#'   Relapse = sample(c("Yes", "No"), 300, replace = TRUE)
#' )
#'
#' # Generate table with medians and column percentages
#' summary_table <- summary_mixed_factorlist(
#'   data = mock_data,
#'   dependent = "Group",
#'   explanatory = c("HeartRate", "Relapse"),
#'   random_effect = "PatientID",
#'   cont_summary = "median",
#'   col_pct = TRUE
#' )
#' 
#' print(summary_table)
#' }


library(lme4)
library(car)
library(dplyr)
library(tidyr)
library(purrr)
###############################################################################
summary_mixed_factorlist <- function(data, dependent, explanatory, random_effect, 
                                     cont_summary = "mean", digits = 1, 
                                     col_pct = TRUE, total_col = TRUE, add_missing = TRUE) {
  
  if (!cont_summary %in% c("mean", "median")) {
    stop("cont_summary must be either 'mean' or 'median'")
  }
  
  dep_sym <- sym(dependent)
  dep_levels <- levels(as.factor(data[[dependent]]))
  results_list <- list()
  
  for (exp_var in explanatory) {
    is_numeric <- is.numeric(data[[exp_var]]) && !is.factor(data[[exp_var]]) && length(unique(data[[exp_var]])) > 2
    p_val_str <- "Error"
    
    # --- STEP 1: HYPOTHESIS TESTING ---
    model_data <- data %>% filter(!is.na(!!sym(exp_var)), !is.na(!!dep_sym))
    
    tryCatch({
      if (is_numeric) {
        if (cont_summary == "mean") {
          # Parametric: Standard Linear Mixed Model
          formula_str <- paste0("`", exp_var, "` ~ `", dependent, "` + (1 | `", random_effect, "`)")
          mod <- lme4::lmer(as.formula(formula_str), data = model_data, REML = FALSE)
          p_val <- car::Anova(mod, type = "III")$`Pr(>Chisq)`[1]
        } else {
          # Non-Parametric: Rank-Transformed Linear Mixed Model
          # This handles skewed data AND missing/unbalanced repeated measures perfectly.
          formula_str <- paste0("rank(`", exp_var, "`) ~ `", dependent, "` + (1 | `", random_effect, "`)")
          mod <- lme4::lmer(as.formula(formula_str), data = model_data, REML = FALSE)
          p_val <- car::Anova(mod, type = "II")$`Pr(>Chisq)`[1]
        }
      } else {
        # Categorical: Generalized Linear Mixed Model (Binomial)
        formula_str <- paste0("as.factor(`", exp_var, "`) ~ `", dependent, "` + (1 | `", random_effect, "`)")
        mod <- lme4::glmer(as.formula(formula_str), data = model_data, family = binomial)
        p_val <- car::Anova(mod, type = "II")$`Pr(>Chisq)`[1]
      }
      p_val_str <- if(p_val < 0.001) "<0.001" else sprintf("%.3f", p_val)
    }, error = function(e) {
      p_val_str <<- "N/A"
    })
    
    # --- STEP 2: SUMMARY STATISTICS STRATIFIED BY DEPENDENT ---
    if (is_numeric) {
      base_summary <- data %>%
        filter(!is.na(!!sym(exp_var)), !is.na(!!dep_sym)) %>%
        group_by(!!dep_sym)
      
      if (cont_summary == "mean") {
        summary_df <- base_summary %>%
          summarise(val = sprintf(paste0("%.", digits, "f (%.", digits, "f)"), mean(!!sym(exp_var)), sd(!!sym(exp_var))), .groups = 'drop') %>%
          pivot_wider(names_from = !!dep_sym, values_from = val)
        
        if (total_col) {
          total_val <- data %>% filter(!is.na(!!sym(exp_var))) %>%
            summarise(Total = sprintf(paste0("%.", digits, "f (%.", digits, "f)"), mean(!!sym(exp_var)), sd(!!sym(exp_var))))
          summary_df <- bind_cols(summary_df, total_val)
        }
        summary_df <- summary_df %>% mutate(label = exp_var, levels = "Mean (SD)", p = p_val_str)
        
      } else {
        summary_df <- base_summary %>%
          summarise(val = sprintf(paste0("%.", digits, "f [%.", digits, "f, %.", digits, "f]"), 
                                  median(!!sym(exp_var)), quantile(!!sym(exp_var), 0.25), quantile(!!sym(exp_var), 0.75)), .groups = 'drop') %>%
          pivot_wider(names_from = !!dep_sym, values_from = val)
        
        if (total_col) {
          total_val <- data %>% filter(!is.na(!!sym(exp_var))) %>%
            summarise(Total = sprintf(paste0("%.", digits, "f [%.", digits, "f, %.", digits, "f]"), 
                                      median(!!sym(exp_var)), quantile(!!sym(exp_var), 0.25), quantile(!!sym(exp_var), 0.75)))
          summary_df <- bind_cols(summary_df, total_val)
        }
        summary_df <- summary_df %>% mutate(label = exp_var, levels = "Median [IQR]", p = p_val_str)
      }
      
      if (add_missing) {
        missing_counts <- data %>%
          group_by(!!dep_sym) %>%
          summarise(n_miss = sum(is.na(!!sym(exp_var))), .groups = 'drop') %>%
          pivot_wider(names_from = !!dep_sym, values_from = n_miss, names_prefix = "miss_")
        
        missing_row <- data.frame(label = "", levels = "Missing")
        for (lev in dep_levels) {
          missing_row[[lev]] <- as.character(missing_counts[[paste0("miss_", lev)]])
        }
        if (total_col) {
          missing_row[["Total"]] <- as.character(sum(is.na(data[[exp_var]])))
        }
        missing_row[["p"]] <- ""
        summary_df <- bind_rows(summary_df, missing_row)
      }
      
    } else {
      cat_data <- data %>% filter(!is.na(!!sym(exp_var)), !is.na(!!dep_sym))
      counts <- table(cat_data[[exp_var]], cat_data[[dependent]])
      totals_row <- rowSums(counts)
      
      pct_matrix <- if (col_pct) {
        prop.table(counts, margin = 2) * 100 
      } else {
        prop.table(cbind(counts, Total = totals_row)[, 1:length(dep_levels)], margin = 1) * 100 
      }
      
      levels_exp <- rownames(counts)
      summary_rows <- list()
      
      for (i in seq_along(levels_exp)) {
        row_vals <- list(label = if(i == 1) exp_var else "", levels = levels_exp[i])
        for (lev in dep_levels) {
          row_vals[[lev]] <- sprintf(paste0("%d (%.", digits, "f%%)"), counts[levels_exp[i], lev], pct_matrix[levels_exp[i], lev])
        }
        if (total_col) {
          tot_pct <- (totals_row[i] / sum(totals_row)) * 100
          row_vals[["Total"]] <- sprintf(paste0("%d (%.", digits, "f%%)"), totals_row[i], tot_pct)
        }
        row_vals[["p"]] <- if(i == 1) p_val_str else ""
        summary_rows[[i]] <- as_tibble(row_vals)
      }
      summary_df <- bind_rows(summary_rows)
      
      if (add_missing) {
        missing_counts <- data %>%
          group_by(!!dep_sym) %>%
          summarise(n_miss = sum(is.na(!!sym(exp_var))), .groups = 'drop') %>%
          pivot_wider(names_from = !!dep_sym, values_from = n_miss, names_prefix = "miss_")
        
        missing_row <- data.frame(label = "", levels = "Missing")
        for (lev in dep_levels) {
          missing_row[[lev]] <- as.character(missing_counts[[paste0("miss_", lev)]])
        }
        if (total_col) {
          missing_row[["Total"]] <- as.character(sum(is.na(data[[exp_var]])))
        }
        missing_row[["p"]] <- ""
        summary_df <- bind_rows(summary_df, missing_row)
      }
    }
    
    results_list[[exp_var]] <- summary_df
  }
  
  final_table <- bind_rows(results_list)
  target_cols <- c("label", "levels", dep_levels)
  if (total_col) target_cols <- c(target_cols, "Total")
  target_cols <- c(target_cols, "p")
  
  final_table <- final_table %>% select(all_of(target_cols))
  colnames(final_table)[1] <- paste0("Dependent: ", dependent)
  
  return(final_table)
}


###################################################################################

set.seed(42)
mock_data <- data.frame(
  PatientID = rep(1:100, each = 3),
  Timepoint = rep(c("Baseline", "Month 1", "Month 6"), 100),
  TreatmentGroup = rep(sample(c("Drug A", "Placebo"), 100, replace = TRUE), each = 3),
  BloodPressure = rnorm(300, mean = 120, sd = 15),
  SymptomPresent = sample(c("Yes", "No" , "Maybe"), 300, replace = TRUE)
)

# Inject random NAs to test feature updates
mock_data$BloodPressure[sample(1:300, 15)] <- NA
mock_data$SymptomPresent[sample(1:300, 22)] <- NA

explanatory_vars <- c("BloodPressure", "SymptomPresent")

# Execute with Column Percentages, Totals, and Missing tracking
table_col_pct <- summary_mixed_factorlist(
  data = mock_data, 
  dependent = "TreatmentGroup", 
  explanatory = explanatory_vars, 
  random_effect = "PatientID",
  cont_summary = "median",
  col_pct = TRUE,      # Column percentages ( % of Treatment Group )
  total_col = TRUE,    # Extra 'Total' column
  add_missing = TRUE   # Extra 'Missing' row per variable block
)

print(table_col_pct)


###########################################################################################
exp_var = "BloodPressure"
dependent = "TreatmentGroup"
random_effect = "PatientID"
model_data = mock_data

formula_str <- paste0("rank(`", exp_var, "`) ~ `", dependent, "` + (1 | `", random_effect, "`)")
mod <- lme4::lmer(as.formula(formula_str), data = model_data, REML = FALSE)

mod
car::Anova(mod , type = 'III')
