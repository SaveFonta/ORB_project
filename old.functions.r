adj_univariate_old <- function(mi_results, delta = 0.5, select_type = "zscore", theta_col = "log_RR2", se_col = "se_RR2") {
  
  draws <- mi_results$draws           
  data_ordered <- mi_results$data_ordered
  unrep_idx <- mi_results$unrep_idx
  M <- nrow(draws)
  
  theta_MA_m <- numeric(M)
  var_MA_m <- numeric(M)
  log_weights <- numeric(M)
  valid <- rep(FALSE, M)
  
  for (m in 1:M) {
    # Create the m-th complete dataset
    current_data <- data_ordered
    current_data[[theta_col]][unrep_idx] <- draws[m, ]
    
    #now the df is complete

    # Run Meta-MA on complete 
    res_m <- tryCatch(
      rma(yi = current_data[[theta_col]], sei = current_data[[se_col]], method = "REML"),
      error = function(e) NULL
    )
    if (is.null(res_m)) next
    
    valid[m] <- TRUE
    
    
    theta_MA_m[m] <- as.numeric(res_m$beta)
    var_MA_m[m] <- res_m$se^2
    
    # Calculate Weight (Eq. 15 / Eq. 17)
    imputed_theta <- draws[m, ]
    if (select_type == "zscore") {
      z_vals <- imputed_theta / data_ordered[[se_col]][unrep_idx]
      log_weights[m] <- -delta * sum(z_vals)
    } else {
      log_weights[m] <- -delta * sum(imputed_theta)
    }
  }
  
  #take only valid iterations
  theta_MA_m  <- theta_MA_m[valid]
  var_MA_m    <- var_MA_m[valid]
  log_weights <- log_weights[valid]


  # Normalize weights (log-sum-exp)
  max_log <- max(log_weights)
  w_norm <- exp(log_weights - max_log) / sum(exp(log_weights - max_log))
  
  # Adjusted (Eq. 16)
  theta_adj <- sum(w_norm * theta_MA_m)
  
  # Rubin's Rules for Variance (Eq. 18)
  var_within <- sum(w_norm * var_MA_m)
  var_between <- sum(w_norm * (theta_MA_m - theta_adj)^2)
  total_var <- var_within + var_between
  
return(data.frame(
    Approach   = paste("Selection on", select_type),
    Estimate   = theta_adj,                             # log scale
    SE         = sqrt(total_var),
    CI_Lower   = theta_adj - 1.96 * sqrt(total_var),
    CI_Upper   = theta_adj + 1.96 * sqrt(total_var)
))
}

adj_bivariate <- function(mi_results, delta = 0.5, select_type = "zscore") {

  draws             <- mi_results$draws
  alternated.vector <- mi_results$alternated.vector
  unrep_idx         <- mi_results$unrep_idx
  V_full            <- mi_results$V_full
  M                 <- nrow(draws)

  # extract se missing to compute the z values
  sei_unrep <- alternated.vector$sei[unrep_idx]  

  log_weights <- if (select_type == "zscore") {
    -delta * as.numeric(draws %*% (1 / sei_unrep))
  } else {
    -delta * rowSums(draws)
  }

  # --- (2) Fit rma.mv on each completed dataset ---
  fits <- lapply(seq_len(M), function(m) {
    d <- alternated.vector
    d$yi[unrep_idx] <- draws[m, ]
    tryCatch(
      rma.mv(yi, V = V_full, mods = ~ outcome - 1,
             random = ~ outcome | study_id, struct = "UN",
             data = d, method = "REML",
             control = list(rel.tol = 1e-8, maxiter = 1000)),
      error = function(e) NULL
    )
  })

  # --- (3) Drop failed fits ---
  valid       <- !sapply(fits, is.null)
  fits        <- fits[valid]
  log_weights <- log_weights[valid]

  # --- (4) Extract both outcomes at once with sapply ---
  theta_MA_m <- t(sapply(fits, function(r) as.numeric(r$beta)))  # M x 2
  var_MA_m   <- t(sapply(fits, function(r) diag(r$vb)))           # M x 2

  # Normalize weights
  w_norm <- { lw <- log_weights - max(log_weights); exp(lw) / sum(exp(lw)) }

  # Adjusted estimates (vectorized across both outcomes)
  theta_adj <- colSums(w_norm * theta_MA_m)               # length-2 vector
  var_within  <- colSums(w_norm * var_MA_m)
  var_between <- colSums(w_norm * sweep(theta_MA_m, 2, theta_adj, "-")^2)
  total_var   <- var_within + var_between

  data.frame(
    Outcome  = c("O1", "O2"),
    Approach = paste("Bivariate Selection on", select_type),
    Estimate = theta_adj,
    SE       = sqrt(total_var),
    CI_Lower = theta_adj - 1.96 * sqrt(total_var),
    CI_Upper = theta_adj + 1.96 * sqrt(total_var)
  )
}