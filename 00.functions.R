
library(MASS)
library(Matrix)
library(metafor)
















# ------------------------------------------
# Univariate
# -----------------------------------------


# Univariate Imputation under MAR
run_univariate_imputation <- function(data, theta_col = "log_RR2", se_col = "se_RR2", m = 1000, new_version = FALSE, method.re = "REML") {
  
  # Find reported and unreported indexes
  rep_idx_start <- which(!is.na(data[[theta_col]]))
  unrep_idx_start <- which(is.na(data[[theta_col]]))
  
  # Order correctly (first  rep, then unrep)
  data_rep <- data[rep_idx_start, ]
  data_unrep <- data[unrep_idx_start, ]
  data_ordered <- rbind(data_rep, data_unrep)
  
  n_rep <- nrow(data_rep)
  n_unrep <- nrow(data_unrep)
  K <- nrow(data_ordered)
  
  # Indexes after reordering
  rep_idx <- 1:n_rep
  unrep_idx <- (n_rep + 1):K
  
  # Fit REML on the reported data (Naive Estimate)
  res_naive <- rma(yi = data_rep[[theta_col]], sei = data_rep[[se_col]], method = method.re)
  theta_MA <- as.numeric(res_naive$beta)
  tau2_hat <- as.numeric(res_naive$tau2)
  SE2_mean <- res_naive$se^2



  J_K <- matrix(1, K, K) 

  #Total cov matrix
  Sigma <- diag(data_ordered[[se_col]]^2) + (SE2_mean * J_K) #paper version
  
  if (new_version) Sigma <- diag(data_ordered[[se_col]]^2 + tau2_hat) + (SE2_mean * J_K) 


  # Partitioning (Eq. 9)
  Sigma_RR <- Sigma[rep_idx, rep_idx]
  Sigma_UU <- Sigma[unrep_idx, unrep_idx]
  Sigma_UR <- Sigma[unrep_idx, rep_idx]
  Sigma_RU <- Sigma[rep_idx, unrep_idx]
  
  # Conditionals (Eq. 11) 
  theta_R <- data_rep[[theta_col]]
  mu_cond <- as.numeric(theta_MA + Sigma_UR %*% solve(Sigma_RR) %*% (theta_R - theta_MA))
  Sigma_cond <- Sigma_UU - Sigma_UR %*% solve(Sigma_RR) %*% Sigma_RU
  
  # Generate M Imputations (we get a matrix m x n_unrep)
  imputed_draws <- mvrnorm(n = m, mu = mu_cond, Sigma = Sigma_cond)

 #probably mvnorm create a vector if we only have one unreported study and can create problems 
  if (n_unrep == 1) {
  imputed_draws <- matrix(imputed_draws, ncol = 1)
}


  return(list(
    draws = imputed_draws,
    unrep_idx = unrep_idx,
    res_naive = res_naive,
    data_ordered = data_ordered,
    theta_col    = theta_col,  
    se_col       = se_col     
  ))
}

#  Importance Sampling 
adj_univariate <- function(mi_results, delta = 0.5, select_type = "zscore", method.re = "REML") {

  draws        <- mi_results$draws
  data_ordered <- mi_results$data_ordered
  unrep_idx    <- mi_results$unrep_idx
  M            <- nrow(draws)

  theta_col    <- mi_results$theta_col  #  from mi_results
  se_col       <- mi_results$se_col     #  from mi_results

  # extract the se useful for the z scores
  sei_unrep <- data_ordered[[se_col]][unrep_idx]   

if (select_type == "zscore") {
    # divides every row in 'draws' by the 'sei_unrep' vector
    z_matrix <- sweep(draws, 2, sei_unrep, FUN = "/")
    log_weights <- -delta * rowSums(z_matrix)
  } else {
    log_weights <- -delta * rowSums(draws)
  }

  # rma on each complete dataset 
  fits <- lapply(seq_len(M), function(m) {
    d <- data_ordered
    d[[theta_col]][unrep_idx] <- draws[m, ] #complete df
    tryCatch(
      rma(yi = d[[theta_col]], sei = d[[se_col]], method = method.re),
      error = function(e) NULL
    )
  })

  # only take not failed 
  valid       <- !sapply(fits, is.null)
  fits        <- fits[valid]
  log_weights <- log_weights[valid]

  # extract estimates
  theta_MA_m <- sapply(fits, function(r) as.numeric(r$beta))
  var_MA_m   <- sapply(fits, function(r) r$se^2)

  # Normalize weights (log-sum-exp)
  max_log <- max(log_weights)
  w_norm <- exp(log_weights - max_log) / sum(exp(log_weights - max_log))


  # Adjusted (Eq. 16)
  theta_adj  <- sum(w_norm * theta_MA_m)

  # Rubin's Rules (Eq. 18)
  var_within <- sum(w_norm * var_MA_m)
  var_between <- sum(w_norm * (theta_MA_m - theta_adj)^2)
  total_var <- var_within + var_between

  data.frame(
    Approach = paste("Selection on", select_type),
    Estimate = theta_adj,
    SE       = sqrt(total_var),
    CI_Lower = theta_adj - 1.96 * sqrt(total_var),
    CI_Upper = theta_adj + 1.96 * sqrt(total_var)
  )
}
















































# ---------------------------
## Bivariate functions 
# ----------------------------

run_bivariate_imputation <- function(data, theta_cols, se_cols, rho_w = 0.7, m = 1000, new_version = FALSE, method.re = "REML") {
  
  K <- nrow(data)
  
  # study-major vector 1x2K
  y_vec <- as.numeric(t(data[, theta_cols]))
  se_vec <- as.numeric(t(data[, se_cols]))
  
  rep_idx <- which(!is.na(y_vec))
  unrep_idx <- which(is.na(y_vec))
  
  # builde the V's (within study cov)
  V_list <- list()
  for (i in 1:K) { # one for each study 
    v1 <- data[[se_cols[1]]][i]^2
    v2 <- data[[se_cols[2]]][i]^2
    cov_12 <- rho_w * sqrt(v1) * sqrt(v2) 
    V_list[[i]] <- matrix(c(v1, cov_12, cov_12, v2), 2, 2)
  }

  V_full <- as.matrix(bdiag(V_list))
  
  # Create Long Data 
  # --> df that has row 1 = study 1 outcome 1, row 2 = study 1 outcome 2, row 3 = study 2 outcome 1 ....

  alternated.vector <- data.frame(
    study_id = rep(data$study_id, each = 2),
    outcome  = factor(rep(c("O1", "O2"), times = K)),
    yi       = y_vec,
    sei      = se_vec,
    obs_id   = 1:(2*K)
  )
  
  alternated.vector_rep <- alternated.vector[rep_idx, ]
  V_rep <- V_full[rep_idx, rep_idx]
  
  # Fit Bivariate REML  (Eq. 12 & 13)
  res_naive <- rma.mv(yi, V = V_rep, 
                      mods = ~ outcome - 1, 
                      random = ~ outcome | study_id, struct = "UN",  # true effects varies for each study id (heterogeneity)
                      data = alternated.vector_rep, method = method.re, 
                      control=list(rel.tol=1e-8, maxiter=1000))
  
  # Construct Total Covariance Matrix Sigma
  Cov_theta_MA <- vcov(res_naive)
  J_K <- matrix(1, K, K)
  
  #  Eq. 13
  Sigma <- V_full + kronecker(J_K, Cov_theta_MA)

  if (new_version) { 
  tau2_1 <- res_naive$tau2[1]
  tau2_2 <- res_naive$tau2[2]
  rho_b  <- res_naive$rho
  
  Psi <- matrix(c(
    tau2_1, rho_b * sqrt(tau2_1) * sqrt(tau2_2),
    rho_b * sqrt(tau2_1) * sqrt(tau2_2), tau2_2
  ), nrow = 2, ncol = 2)
    I_K <- diag(K)
  
  Sigma <- V_full + kronecker(I_K, Psi) + kronecker(J_K, Cov_theta_MA)
  }

  # partitioning 
  Sigma_RR <- Sigma[rep_idx, rep_idx]
  Sigma_UU <- Sigma[unrep_idx, unrep_idx]
  Sigma_UR <- Sigma[unrep_idx, rep_idx]
  Sigma_RU <- Sigma[rep_idx, unrep_idx]
  
  #  Conditionals (Eq. 11)
  theta_MA_vec <- rep(coef(res_naive), times = K) 
  theta_R_MA <- theta_MA_vec[rep_idx]
  theta_U_MA <- theta_MA_vec[unrep_idx]
  theta_R <- y_vec[rep_idx]
  
  mu_cond <- theta_U_MA + Sigma_UR %*% solve(Sigma_RR) %*% (theta_R - theta_R_MA)
  Sigma_cond <- Sigma_UU - Sigma_UR %*% solve(Sigma_RR) %*% Sigma_RU
  
  # Generate M Imputations
  imputed_draws <- mvrnorm(n = m, mu = as.numeric(mu_cond), Sigma = Sigma_cond)
  if (length(unrep_idx) == 1) {
    imputed_draws <- matrix(imputed_draws, ncol = 1)
  }
  
  return(list(
    draws = imputed_draws,
    unrep_idx = unrep_idx,
    res_naive = res_naive,
    alternated.vector = alternated.vector,
    V_full = V_full,
    theta_cols = theta_cols
  ))
}



adj_bivariate <- function(mi_results, delta = 0.5, select_type = "zscore", method.re = "REML") {

  draws <- mi_results$draws
  alternated.vector <- mi_results$alternated.vector
  unrep_idx <- mi_results$unrep_idx
  V_full <- mi_results$V_full
  M <- nrow(draws)

  # extract the se useful for the z scores
  sei_unrep <- alternated.vector$sei[unrep_idx]  


      # NOTE: here the paper is not clear. 
    # I decided to have one weight for each imputation m. otherwise we'd have w_j^(m) but this mean two weight for each 
    # imputation in the topiramate case. (this choice doesnt affect simulation study)


  if (select_type == "zscore") {
    # divides every row in 'draws' by the 'sei_unrep' vector
    z_matrix <- sweep(draws, 2, sei_unrep, FUN = "/")
    log_weights <- -delta * rowSums(z_matrix)
  } else {
    log_weights <- -delta * rowSums(draws)
  }

  # loop for each draw
  fits <- lapply(seq_len(M), function(m) {
    d <- alternated.vector

    # create complete dataset
    d$yi[unrep_idx] <- draws[m, ]

    #fit the results
    tryCatch( 
      rma.mv(yi, V = V_full, mods = ~ outcome - 1,
             random = ~ outcome | study_id, struct = "UN",
             data = d, method = method.re,
             control = list(rel.tol = 1e-8, maxiter = 1000)),
      error = function(e) NULL
    )
  })

  # extract valid
  valid       <- !sapply(fits, is.null)
  fits        <- fits[valid]
  log_weights <- log_weights[valid]

  # Extract 
  theta_MA_m <- t(sapply(fits, function(r) as.numeric(r$beta)))  # M x 2
  var_MA_m   <- t(sapply(fits, function(r) diag(r$vb)))           # M x 2

  # Normalize weights
  max_log <- max(log_weights)
  w_norm <- exp(log_weights - max_log) / sum(exp(log_weights - max_log))
  
  # Adjusted estimates 
  theta_adj <- colSums(w_norm * theta_MA_m)               # length-2 vector

  # Rubin rule
  var_within  <- colSums(w_norm * var_MA_m)
  var_between <- colSums(w_norm * sweep(theta_MA_m, 2, theta_adj, "-")^2)
  total_var   <- var_within + var_between

  results_df <- data.frame(
    Outcome  = c("O1", "O2"),
    Approach = paste("Bivariate Selection on", select_type),
    Estimate = theta_adj,
    SE       = sqrt(total_var),
    CI_Lower = theta_adj - 1.96 * sqrt(total_var),
    CI_Upper = theta_adj + 1.96 * sqrt(total_var)
  )

    rownames(results_df) <- NULL

    return(results_df)
}









































# --------------------------------
# Simulate data
# -------------------------------

# here we generate the full data 

# K = number of studies
# theta = vector containing the true values for the two outcomes
# tau2 = vector containing the true tau
# rho_b is the between study correlation
# rho_w withing study correlation 
# n_arm is fixed to 50 in the whole simulation 


generate_bivariate_ma <- function(K = 12, theta = c(0.4, 0.4), tau2 = c(0.06, 0.06), 
                                  rho_b = 0.4, rho_w = 0.4, n_arm = 50) {
  
  # Between-study covariance matrix (Psi)
  cov_b <- rho_b * sqrt(tau2[1]) * sqrt(tau2[2])
  Psi <- matrix(c(tau2[1], cov_b, cov_b, tau2[2]), 2, 2)
  
  # True study-specific effects (bivariate normal)
  theta_i <- mvrnorm(n = K, mu = theta, Sigma = Psi)
  
  
  #  matrix V for Wishart distribution
  V_scale <- (1 / ((n_arm - 1) * n_arm)) * matrix(c(1, rho_w, rho_w, 1), 2, 2)
  df_wishart <- 2 * (n_arm - 1)
  
  # Initialize vectors for obs effects and standard errors
  y_obs <- matrix(NA, nrow = K, ncol = 2)
  se_obs <- matrix(NA, nrow = K, ncol = 2)
  
  for (i in 1:K) {
    #  within-study covariance matrix for study i using Wishart (Sigma_i)
    Sigma_i <- rWishart(n = 1, df = df_wishart, Sigma = V_scale)[,,1]
    
    # Draw observed effects (bivariate)
    y_obs[i, ] <- mvrnorm(n = 1, mu = theta_i[i, ], Sigma = Sigma_i)
    
    # Extract standard errors
    se_obs[i, 1] <- sqrt(Sigma_i[1, 1])
    se_obs[i, 2] <- sqrt(Sigma_i[2, 2])
  }
  
  data <- data.frame(
    study_id = 1:K,
    O1_yi = y_obs[, 1], O1_sei = se_obs[, 1],
    O2_yi = y_obs[, 2], O2_sei = se_obs[, 2]
  )

  attr(data, "theta1") <- theta[1]
  attr(data, "tau2_1") <- tau2[1]
  attr(data, "n_arm") <- n_arm

  return(data)
}


















expit <- function(x) exp(x) / (1 + exp(x))
logit <- function(p) log(p / (1 - p))

# Now we need to choose:
# p1 = prob of unreporting
# delta value used in the simulation 
# selection type 



impose_orb <- function(data, p1 = 0.2, delta_sim = 0.5, select_type = "zscore") {
  
  #extract stuff from attributes
  theta1 <- attr(data, "theta1")
  tau2_1 <- attr(data, "tau2_1")
  n_arm <- attr(data, "n_arm")




  # Determine the selection variable s_1 and its expected value E(s_1)
  if (select_type == "zscore") {
    s1 <- data$O1_yi / data$O1_sei #vector of real value 
    E_s1 <- theta1 / sqrt((2/n_arm) + tau2_1) #expected score to define teoreticcal alpha
  } else if (select_type == "effect") {
    s1 <- data$O1_yi 
    E_s1 <- theta1 # if we do selection on effect, the exp value is directly the true value
  } else {
    stop("select_type must be 'zscore' or 'effect'")
  }
  
  target_reported <- 1 - p1
  
  # Avoid logit of 1 or 0 if p1 is extreme (this should never happen in the code)
  if(target_reported >= 1) target_reported <- 0.999 
  if(target_reported <= 0) target_reported <- 0.001
  
  alpha_1 <- logit(target_reported) - (delta_sim * E_s1)
  
  # Calculate reporting probabilities for each outocme
  prob_report <- expit(alpha_1 + delta_sim * s1)
  
  # Simulate reporting or not (1 = reported, 0 = missing)
  reported_flag <- rbinom(n = nrow(data), size = 1, prob = prob_report)
  
  # Apply missingness to Outcome 1
  data_orb <- data
  data_orb$O1_yi[reported_flag == 0] <- NA
 
 #QUESTION-->  should I also unreport the SE?? I think it is just one addditional step of inputing SEs...
 #data_orb$O1_sei[reported_flag == 0] <- NA 
  
  return(data_orb)
}













