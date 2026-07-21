

### EXAMPLE CAM WORKFLOW FOR A TWO STEP PIPELINE #############################


# Step 1: Find optimal K and save the fitted CAM object
#opt_res <- find_optimal_K(data = expr_matrix, K_range = 2:8)

# Step 2: Extract results using the optimal K (or override with any K you want)
#final_results <- extract_CAM3_results(
#  cam_obj = opt_res$cam_object,
#  data    = expr_matrix,
#  K       = opt_res$optimal_K
#)

# Inspect outputs:
#head(final_results$A_matrix)
#head(final_results$SG_averages)





#' Find Optimal K for CAM3 using Minimum Description Length (MDL)
#' 
#' @param data A gene expression matrix (Rows = Genes, Columns = Observations/Samples).
#' @param K_range A numeric vector specifying the range of K to test (e.g., 2:10).
#' @param dim.rdc Reduced data dimension; should be not less than maximum candidate K.
#' @param plot_MDL Logical; if TRUE, plots the MDL curve and highlights the minimum.
#' @param ... Additional arguments passed to the `CAM3()` function.
#' 
#' @return A list containing:
#'   - \code{optimal_K}: The integer K that minimizes the MDL.
#'   - \code{MDL_values}: A named vector of MDL values for each tested K.
#'   - \code{input_parameters}: A list of all CAM parameters used (excluding the data matrix).
#'   - \code{cam_object}: The complete CAM object.
#'   - \code{params}: A list of the parameters passed to CAM.
find_optimal_K <- function(data, K_range = 2:8, dim.rdc = 10 ,  plot_MDL = TRUE, ...) {
  
  # 1. Capture the input parameters (saving K_range and all arguments passed via '...')
  additional_params <- list(...)
  cam_params <- c(list(K = K_range), additional_params)
  
  # 2. Check for the package namespace
  pkg <- NULL
  if (requireNamespace("CAM3", quietly = TRUE)) {
    pkg <- "CAM3"
  } else if (requireNamespace("debCAM", quietly = TRUE)) {
    pkg <- "debCAM"
  } else {
    stop("Please install 'debCAM' from Bioconductor or 'CAM3' from GitHub.")
  }
  
  # 3. Run CAM across the range of K
  cam_fun <- get("CAM3Run", envir = asNamespace(pkg))
  message(sprintf("Running %s::CAM for K = %s to %s...", 
                  pkg, min(K_range), max(K_range)))
  
  # Execute CAM using the data, K_range, and the captured '...' parameters
  cam_res <- cam_fun(data, K = K_range, ...)
  
  # 4. Extract the MDL values
  mdl_vals  =  vector()
  i =  1
  for (k  in K_range){
    mdl_vals[i]=cam_res@ASestResult[[as.character(k)]]@mdl
    i =  i + 1
  }
  
  # Ensure the vector is named properly by K
  if (is.null(names(mdl_vals))) {
    names(mdl_vals) <- as.character(K_range)
  }
  
  # 5. Identify the Optimal K (minimum MDL)
  opt_K_idx <- which.min(mdl_vals)
  opt_K <- as.numeric(names(mdl_vals)[opt_K_idx])
  message(sprintf("Optimal K found: %d", opt_K))
  
  # 6. Plot the MDL curve
  if (plot_MDL) {
    plot(x = as.numeric(names(mdl_vals)), y = mdl_vals, 
         type = "b", pch = 19, col = "darkgray", lwd = 2,
         xlab = "Number of Sources (K)", 
         ylab = "Minimum Description Length (MDL)",
         main = "Optimal K Selection via MDL Curve",
         xaxt = "n")
    axis(1, at = as.numeric(names(mdl_vals)))
    
    # Highlight the optimal K in red
    points(opt_K, mdl_vals[opt_K_idx], col = "red", pch = 19, cex = 1.5)
    legend("topright", legend = paste("Optimal K =", opt_K), 
           col = "red", pch = 19, bty = "n")
  }
  
  # 7. Return combined results including the parameter list
  return(list(
    optimal_K = opt_K,
    MDL_values = mdl_vals,
    params = cam_params,
    cam_object = cam_res
  ))
}





#' Run CAM3 and Extract Matrices, Marker Genes, and SG Averages
#' 
#' @param data A gene expression matrix (Rows = Genes, Columns = Observations/Samples).
#' @param K An integer specifying the number of subpopulations/sources.
#' @param dim.rdc Reduced data dimension; should be not less than maximum candidate K.
#' @param cos.thres cosine similarity threshold for marker gene selection. MGs with cos similarity greater than this threshold are kept .default is 0.95
#' @param ... Additional arguments passed to the `CAM3()` function.
#' 
#' @return A list containing:
#'   - \code{A_matrix}: Estimated proportion matrix.
#'   - \code{S_matrix}: Estimated signature matrix.
#'   - \code{SG_list}: List of Specific Genes (marker genes) for each source at the specified cos threshold. .
#'   - \code{SG_averages}: Matrix of (Observations x Sources) with the average SG expression.
#'   - \code{params}: A list of the parameters passed to CAM.
run_CAM3_wrapper <- function(data, K, dim.rdc = 10, cos.thres = 0.95,  ...) {
  # 1. Check for the package (looks for CAM3 first, then debCAM)
  pkg <- NULL
  if (requireNamespace("CAM3", quietly = TRUE)) {
    pkg <- "CAM3"
  } else if (requireNamespace("debCAM", quietly = TRUE)) {
    pkg <- "debCAM"
  } else {
    stop("Please install 'debCAM' from Bioconductor or 'CAM3' from GitHub to run this wrapper.")
  }
  
  
  
  # Capture the input parameters (saving K_range and all arguments passed via '...')
  additional_params <- list(...)
  cam_params <- c(list(K = K , dim.rdc = dim.rdc , cos.thres =  cos.thres), additional_params)
  
  
  # 2. Run the CAM algorithm
  cam_fun <- get("CAM3Run", envir = asNamespace(pkg))
  message(sprintf("Running %s::CAM with K = %d...", pkg, K))
  cam_res <- cam_fun(data, K = K, ...)
  
  # 3. Extract A and S matrices for the specified K
  # Uses standard accessors Amat()/Smat() or falls back to S4 slot extraction
  A_mat <- tryCatch({
    get("Amat", envir = asNamespace(pkg))(cam_res, K)
  }, error = function(e) {
    cam_res@ASestResult[[as.character(K)]]@Aest
  })
  
  S_mat <- tryCatch({
    get("Smat", envir = asNamespace(pkg))(cam_res, K)
  }, error = function(e) {
    cam_res@ASestResult[[as.character(K)]]@Sest
  })
  
  # 4. Extract SG (Specific Genes / Marker Genes) list
  SG_list = cotMG(Sest = S_mat, cos.thres = cos.thres, thres.low = 0, thres.high = 1)[[1]]
  
  # 5. Calculate the average of all SGs per Source per observation in the original data
  num_obs <- ncol(data)
  num_sources <- length(SG_list)
  
  # Ensure the sources have readable names
  source_names <- names(SG_list)
  if (is.null(source_names)) {
    source_names <- paste0("Source_", seq_len(num_sources))
    names(SG_list) <- source_names
  }
  
  # Initialize the SG_averages matrix (Rows = Samples/Observations, Cols = Sources)
  SG_averages <- matrix(NA, nrow = num_obs, ncol = num_sources)
  rownames(SG_averages) <- colnames(data)
  colnames(SG_averages) <- source_names
  
  for (k in seq_along(SG_list)) {
    genes_k <- SG_list[[k]]
    
    # Safely match gene indices or names to the original data
    if (is.character(genes_k)) {
      if (is.null(rownames(data))) {
        stop("The dataset must have row names because the SG list contains gene names.")
      }
      valid_genes <- intersect(genes_k, rownames(data))
    } else {
      valid_genes <- intersect(genes_k, seq_len(nrow(data)))
    }
    
    # Calculate column means across the extracted specific genes
    if (length(valid_genes) > 1) {
      SG_averages[, k] <- colMeans(data[valid_genes, , drop = FALSE], na.rm = TRUE)
    } else if (length(valid_genes) == 1) {
      SG_averages[, k] <- data[valid_genes, ]
    } else {
      warning(sprintf("No valid Specific Genes found in the original data for %s", source_names[k]))
      SG_averages[, k] <- NA
    }
  }
  
  # 6. Return combined results
  return(list(
    A_matrix = A_mat,
    S_matrix = S_mat,
    SG_list = SG_list,
    SG_averages = SG_averages,
    params =  cam_params
  ))
}









#' Extract Matrices, Marker Genes, and SG Averages
#' 
#' @param cam_obj An existing CAM object (e.g., from find_optimal_K or debCAM::CAM).
#' @param data A gene expression matrix (Rows = Genes, Columns = Observations/Samples).
#' @param K An integer specifying the number of subpopulations/sources.
#' @param cos.thres cosine similarity threshold for marker gene selection. MGs with cos similarity greater than this threshold are kept .default is 0.95
#'
#' 
#' @return A list containing:
#'   - \code{A_matrix}: Estimated proportion matrix.
#'   - \code{S_matrix}: Estimated signature matrix.
#'   - \code{SG_list}: List of Specific Genes (marker genes) for each source at the specified cos threshold. 
#'   - \code{SG_averages}: Matrix of (Observations x Sources) with the average SG expression.
extract_CAM3_results <- function( cam_obj ,  data, K, cos.thres = 0.95,) {
  # 1. Check for the package (looks for CAM3 first, then debCAM)
  pkg <- NULL
  if (requireNamespace("CAM3", quietly = TRUE)) {
    pkg <- "CAM3"
  } else if (requireNamespace("debCAM", quietly = TRUE)) {
    pkg <- "debCAM"
  } else {
    stop("Please install 'debCAM' from Bioconductor or 'CAM3' from GitHub to run this wrapper.")
  }
  
  cam_res = cam_obj
  
  # 2. Extract A and S matrices for the specified K
  # Uses standard accessors Amat()/Smat() or falls back to S4 slot extraction
  A_mat <- tryCatch({
    get("Amat", envir = asNamespace(pkg))(cam_res, K)
  }, error = function(e) {
    cam_res@ASestResult[[as.character(K)]]@Aest
  })
  
  S_mat <- tryCatch({
    get("Smat", envir = asNamespace(pkg))(cam_res, K)
  }, error = function(e) {
    cam_res@ASestResult[[as.character(K)]]@Sest
  })
  
  # 3. Extract SG (Specific Genes / Marker Genes) list
  SG_list = cotMG(Sest = S_mat, cos.thres = cos.thres, thres.low = 0, thres.high = 1)[[1]]
  
  # 4. Calculate the average of all SGs per Source per observation in the original data
  num_obs <- ncol(data)
  num_sources <- length(SG_list)
  
  # Ensure the sources have readable names
  source_names <- names(SG_list)
  if (is.null(source_names)) {
    source_names <- paste0("Source_", seq_len(num_sources))
    names(SG_list) <- source_names
  }
  
  # Initialize the SG_averages matrix (Rows = Samples/Observations, Cols = Sources)
  SG_averages <- matrix(NA, nrow = num_obs, ncol = num_sources)
  rownames(SG_averages) <- colnames(data)
  colnames(SG_averages) <- source_names
  
  for (k in seq_along(SG_list)) {
    genes_k <- SG_list[[k]]
    
    # Safely match gene indices or names to the original data
    if (is.character(genes_k)) {
      if (is.null(rownames(data))) {
        stop("The dataset must have row names because the SG list contains gene names.")
      }
      valid_genes <- intersect(genes_k, rownames(data))
    } else {
      valid_genes <- intersect(genes_k, seq_len(nrow(data)))
    }
    
    # Calculate column means across the extracted specific genes
    if (length(valid_genes) > 1) {
      SG_averages[, k] <- colMeans(data[valid_genes, , drop = FALSE], na.rm = TRUE)
    } else if (length(valid_genes) == 1) {
      SG_averages[, k] <- data[valid_genes, ]
    } else {
      warning(sprintf("No valid Specific Genes found in the original data for %s", source_names[k]))
      SG_averages[, k] <- NA
    }
  }
  
  # 5. Return combined results
  return(list(
    A_matrix = A_mat,
    S_matrix = S_mat,
    SG_list = SG_list,
    SG_averages = SG_averages,
  ))
}




