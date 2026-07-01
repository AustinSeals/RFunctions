library(dplyr)
library(fastDummies)
library(tibble)
library(stringi)

simple_stats = function(dat, dummy_cut = 8, grouper, sort_categories = FALSE) {
  
  dat = dat[!is.na(dat[, grouper]), ]
  
  ### Get number of distinct values per column
  cols = names(dat)
  grouper_n = unique(dat[, grouper])
  dist = sapply(dat, function(x) n_distinct(x))
  
  ### Get each column's data type
  types = sapply(dat, function(x) class(x)[1]) 
  message("Detected Data Types:")
  print(types)
  
  dat = dat %>% mutate_if(is.factor, as.character)
  
  ### Create temp dummy column
  dat$temp123 = dat[, grouper]
  
  ### Select which types to dummy
  to_dummy = cols[(types %in% c("numeric", "integer", "logical", "factor", "character")) & (dist %in% c(2:dummy_cut))]
  to_dummy = to_dummy[to_dummy != grouper]
  to_dummy = c(to_dummy, "temp123")
  
  ### Select others
  others = cols[(types %in% c("numeric", "integer")) & (dist > dummy_cut)]
  others = others[others != grouper]
  
  ### Dummify explicitly keeping all levels 
  dummy_data = fastDummies::dummy_cols(
    dat[, c(to_dummy, others)], 
    select_columns = to_dummy, 
    remove_selected_columns = TRUE,
    remove_first_dummy = FALSE
  )
  
  dummy_names = fastDummies::dummy_cols(
    dat[, to_dummy], 
    select_columns = to_dummy, 
    remove_selected_columns = TRUE,
    remove_first_dummy = FALSE
  ) %>% names(.)
  
  # Create a map linking dummy column names back to their original base variable
  base_var_map = character(length(dummy_names))
  names(base_var_map) = dummy_names
  for (v in to_dummy) {
    matches = startsWith(dummy_names, paste0(v, "_"))
    base_var_map[matches] = v
  }
  
  newData = data.frame(Group = dat[, grouper], dummy_data, check.names = FALSE)
  names(newData)[1] = "Group"
  
  # NEW: Create a combined dataset for summarization that includes an "Overall" cohort
  newData_overall = newData
  newData_overall$Group = "Overall"
  newData_combined = bind_rows(newData_overall, newData)
  
  # Lock "Overall" as the first factor level so it appears as the first column
  newData_combined$Group = factor(newData_combined$Group, levels = c("Overall", sort(unique(as.character(newData$Group)))))
  
  message("(1 of 4) Data Wrangling Done")
  
  ### Summarize Data (Calculates overall means & proportions)
  smry = newData_combined %>% 
    group_by(Group) %>%
    summarise_if(is.numeric, mean, na.rm = TRUE) %>%
    tibble::column_to_rownames(var = "Group") %>%
    mutate_all(~ round(., 3)) %>%
    t(.) %>%
    as.data.frame(.) %>%
    tibble::rownames_to_column(var = "Variable")
  
  # Adjusted label names to reflect they apply to the overall column too
  smry = tibble::add_column(smry, Measure = ifelse(smry$Variable %in% dummy_names, "Percent", "Mean"), .after = "Variable")
  
  # Calculate missing ONLY on original data to avoid double counting
  missing = sapply(newData[, smry$Variable], function(x) sum(is.na(x)))
  smry$Missing = missing
  
  message("(2 of 4) Summary Stats Completed")
  
  ### Run Association Test (Strictly run on original newData without the Overall group)
  pv = c()
  total_n = c()
  
  for (i in 1:nrow(smry)) {
    form = as.formula(paste0("`", smry$Variable[i], "` ~ Group"))
    
    if (smry$Measure[i] == "Mean") {
      # ANOVA for parametric continuous variables
      p = summary(aov(form, data = newData))[[1]][["Pr(>F)"]][1]
      n_ = sum(!is.na(newData[, smry$Variable[i]]))
    } else {
      p = tryCatch(
        fisher.test(newData$Group, newData[, smry$Variable[i]])$p.value, 
        error = function(e) { return(1) }
      )
      n_ = sum(newData[, smry$Variable[i]], na.rm = TRUE)
    }
    
    pv = c(pv, p)
    total_n = c(total_n, n_)
  }
  
  smry$P.value = round(pv, 4)
  smry$Sig = ifelse(smry$P.value < 0.05, 1, 0)
  smry = tibble::add_column(smry, Total_N = total_n, .after = "Variable")
  
  message("(3 of 4) Association Test Completed")
  
  ### Format percent data dynamically picking up "Overall" and the distinct groups
  vars = names(smry)[!names(smry) %in% c("Variable", "Measure", "Missing", "P.value", "Sig", "Total_N")]
  
  for(v in vars) {
    smry[, v] = ifelse(smry$Measure == "Percent", round(smry[, v] * 100, 1), smry[, v])  
  }
  
  ### Add N
  N_df = newData_combined %>% 
    group_by(Group) %>% 
    summarise(Group.N = n()) %>% 
    tibble::column_to_rownames(var = "Group") %>% 
    t(.) %>%
    as.data.frame(.) %>%
    tibble::rownames_to_column(var = "Measure") %>%
    mutate_all(as.character)
  
  ## ADD SD Formatting  
  for (i in 1:nrow(smry)) {
    if (smry$Measure[i] == "Mean") {
      
      # Calculate SD natively including the Overall group
      sd_tab = newData_combined %>% 
        group_by(Group) %>% 
        summarise(sd = round(sd(!!sym(smry$Variable[i]), na.rm = TRUE), 1)) %>%
        mutate(sd_str = paste0(" (", sd, ")")) %>%
        select(Group, sd_str) %>%
        tibble::column_to_rownames(var = "Group") %>%
        t(.) %>% 
        as.vector(.)
      
      current_values = smry[i, vars] %>% as.vector(.)
      smry[i, vars] = paste0(as.character(round(as.numeric(current_values), 1)), sd_tab)
      
    } else {
      current_values = smry[i, vars] %>% as.vector(.)
      smry[i, vars] = paste0(as.character(current_values), "%")
    }
  }
  
  final_smry = bind_rows(N_df, smry)
  final_smry = final_smry[, c("Variable", "Total_N", "Measure", vars, "P.value", "Sig")]
  
  message("(4 of 4) Formatting Add group sizes/counts ")
  
  ### remove temp
  final_smry = final_smry %>% filter(!grepl("temp123", Variable))
  
  # Sorting Logic
  if (sort_categories) {
    message("Sorting categories by frequency...")
    
    final_smry = final_smry %>%
      mutate(
        Original_Row = row_number(),
        Base_Variable = ifelse(
          Measure == "Percent" & !is.na(Variable), 
          base_var_map[Variable], 
          coalesce(Variable, "Top_Row")
        ),
        Numeric_N = as.numeric(Total_N) 
      ) %>%
      group_by(Base_Variable) %>%
      mutate(
        Category_Count = n(),
        Var_Order = min(Original_Row),
        Sort_Value = ifelse(Measure == "Percent" & Category_Count > 2, -Numeric_N, Original_Row)
      ) %>%
      ungroup() %>%
      arrange(Var_Order, Sort_Value) %>%
      select(-Original_Row, -Base_Variable, -Numeric_N, -Category_Count, -Var_Order, -Sort_Value)
  }
  
  final_smry$Variable = ifelse(
    final_smry$Measure == "Percent",
    stringi::stri_replace_last_fixed(final_smry$Variable, "_", ": "), 
    final_smry$Variable
  )
  
  message("Done")
  
  return(final_smry)
}


#######################################################################


simple_stats_np = function(dat, dummy_cut = 8, grouper, sort_categories = FALSE) {
  
  dat = dat[!is.na(dat[, grouper]), ]
  
  ### Get number of distinct values per column
  cols = names(dat)
  grouper_n = unique(dat[, grouper])
  dist = sapply(dat, function(x) n_distinct(x))
  
  ### Get each column's data type
  types = sapply(dat, function(x) class(x)[1]) 
  message("Detected Data Types:")
  print(types)
  
  dat = dat %>% mutate_if(is.factor, as.character)
  
  ### Create temp dummy column
  dat$temp123 = dat[, grouper]
  
  ### Select which types to dummy
  to_dummy = cols[(types %in% c("numeric", "integer", "logical", "factor", "character")) & (dist %in% c(2:dummy_cut))]
  to_dummy = to_dummy[to_dummy != grouper]
  to_dummy = c(to_dummy, "temp123")
  
  ### Select others
  others = cols[(types %in% c("numeric", "integer")) & (dist > dummy_cut)]
  others = others[others != grouper]
  
  ### Dummify explicitly keeping all levels 
  dummy_data = fastDummies::dummy_cols(
    dat[, c(to_dummy, others)], 
    select_columns = to_dummy, 
    remove_selected_columns = TRUE,
    remove_first_dummy = FALSE
  )
  
  dummy_names = fastDummies::dummy_cols(
    dat[, to_dummy], 
    select_columns = to_dummy, 
    remove_selected_columns = TRUE,
    remove_first_dummy = FALSE
  ) %>% names(.)
  
  # Create a map linking dummy column names back to their original base variable
  base_var_map = character(length(dummy_names))
  names(base_var_map) = dummy_names
  for (v in to_dummy) {
    matches = startsWith(dummy_names, paste0(v, "_"))
    base_var_map[matches] = v
  }
  
  newData = data.frame(Group = dat[, grouper], dummy_data, check.names = FALSE)
  names(newData)[1] = "Group"
  
  # NEW: Create a combined dataset for summarization that includes an "Overall" cohort
  newData_overall = newData
  newData_overall$Group = "Overall"
  newData_combined = bind_rows(newData_overall, newData)
  
  # Lock "Overall" as the first factor level so it appears as the first column
  newData_combined$Group = factor(newData_combined$Group, levels = c("Overall", sort(unique(as.character(newData$Group)))))
  
  message("(1 of 4) Data Wrangling Done")
  
  ### Summarize Data (Patched to calculate mean for categories, median for continuous)
  smry = newData_combined %>% 
    group_by(Group) %>%
    summarise(
      across(all_of(dummy_names), ~ mean(.x, na.rm = TRUE)),
      across(all_of(others), ~ median(.x, na.rm = TRUE))
    ) %>%
    tibble::column_to_rownames(var = "Group") %>%
    mutate_all(~ round(., 3)) %>%
    t(.) %>%
    as.data.frame(.) %>%
    tibble::rownames_to_column(var = "Variable")
  
  # Adjusted label names to reflect they apply to the overall column too
  smry = tibble::add_column(smry, Measure = ifelse(smry$Variable %in% dummy_names, "Percent", "Median"), .after = "Variable")
  
  # Calculate missing ONLY on original data to avoid double counting
  missing = sapply(newData[, smry$Variable], function(x) sum(is.na(x)))
  smry$Missing = missing
  
  message("(2 of 4) Summary Stats Completed")
  
  ### Run Association Test (Strictly run on original newData without the Overall group)
  pv = c()
  total_n = c()
  
  for (i in 1:nrow(smry)) {
    form = as.formula(paste0("`", smry$Variable[i], "` ~ Group"))
    
    if (smry$Measure[i] == "Median") {
      p = kruskal.test(form, data = newData)$p.value
      n_ = sum(!is.na(newData[, smry$Variable[i]]))
    } else {
      p = tryCatch(
        fisher.test(newData$Group, newData[, smry$Variable[i]])$p.value, 
        error = function(e) { return(1) }
      )
      n_ = sum(newData[, smry$Variable[i]], na.rm = TRUE)
    }
    
    pv = c(pv, p)
    total_n = c(total_n, n_)
  }
  
  smry$P.value = round(pv, 4)
  smry$Sig = ifelse(smry$P.value < 0.05, 1, 0)
  smry = tibble::add_column(smry, Total_N = total_n, .after = "Variable")
  
  message("(3 of 4) Association Test Completed")
  
  ### Format percent data dynamically picking up "Overall" and the distinct groups
  vars = names(smry)[!names(smry) %in% c("Variable", "Measure", "Missing", "P.value", "Sig", "Total_N")]
  
  for(v in vars) {
    smry[, v] = ifelse(smry$Measure == "Percent", round(smry[, v] * 100, 1), smry[, v])  
  }
  
  ### Add N
  N_df = newData_combined %>% 
    group_by(Group) %>% 
    summarise(Group.N = n()) %>% 
    tibble::column_to_rownames(var = "Group") %>% 
    t(.) %>%
    as.data.frame(.) %>%
    tibble::rownames_to_column(var = "Measure") %>%
    mutate_all(as.character)
  
  ## ADD IQR Formatting  
  for (i in 1:nrow(smry)) {
    if (smry$Measure[i] == "Median") {
      
      # Calculate Q1 and Q3 natively including the Overall group
      iqr_tab = newData_combined %>% 
        group_by(Group) %>% 
        summarise(
          q1 = round(quantile(!!sym(smry$Variable[i]), 0.25, na.rm = TRUE), 1),
          q3 = round(quantile(!!sym(smry$Variable[i]), 0.75, na.rm = TRUE), 1)
        ) %>%
        mutate(iqr_str = paste0(" (", q1, " - ", q3, ")")) %>%
        select(Group, iqr_str) %>%
        tibble::column_to_rownames(var = "Group") %>%
        t(.) %>% 
        as.vector(.)
      
      current_values = smry[i, vars] %>% as.vector(.)
      smry[i, vars] = paste0(as.character(round(as.numeric(current_values), 1)), iqr_tab)
      
    } else {
      current_values = smry[i, vars] %>% as.vector(.)
      smry[i, vars] = paste0(as.character(current_values), "%")
    }
  }
  
  final_smry = bind_rows(N_df, smry)
  final_smry = final_smry[, c("Variable", "Total_N", "Measure", vars, "P.value", "Sig")]
  
  message("(4 of 4) Formatting Add group sizes/counts ")
  
  ### remove temp
  final_smry = final_smry %>% filter(!grepl("temp123", Variable))
  
  # Sorting Logic
  if (sort_categories) {
    message("Sorting categories by frequency...")
    
    final_smry = final_smry %>%
      mutate(
        Original_Row = row_number(),
        Base_Variable = ifelse(
          Measure == "Percent" & !is.na(Variable), 
          base_var_map[Variable], 
          coalesce(Variable, "Top_Row")
        ),
        Numeric_N = as.numeric(Total_N) 
      ) %>%
      group_by(Base_Variable) %>%
      mutate(
        Category_Count = n(),
        Var_Order = min(Original_Row),
        Sort_Value = ifelse(Measure == "Percent" & Category_Count > 2, -Numeric_N, Original_Row)
      ) %>%
      ungroup() %>%
      arrange(Var_Order, Sort_Value) %>%
      select(-Original_Row, -Base_Variable, -Numeric_N, -Category_Count, -Var_Order, -Sort_Value)
  }
  
  final_smry$Variable = ifelse(
    final_smry$Measure == "Percent",
    stringi::stri_replace_last_fixed(final_smry$Variable, "_", ": "), 
    final_smry$Variable
  )
  
  message("Done")
  
  return(final_smry)
}

