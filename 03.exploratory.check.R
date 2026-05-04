# ---------------------------------------------------------
# Some exploratory check to see unbiasdness and how it behaves under varying M
# ---------------------------------------------------------
source("00.functions.R")
library(parallel)

RNGkind("L'Ecuyer-CMRG")
set.seed(1)

# Simulation parameters
n_sim <- 1000 
true_theta <- 0.4
n_cores <- 31

M <- 1000 # imputations

scenarios <- expand.grid(
  K     = c(25, 6),
  tau2  = c(0.06, 0.36),
  delta = c(0, 0.5, 1)
)

# function to subset to m = c(50,200,1000) from the generated df of M = 1000 imputations 
subset_mi <- function(mi, m_use) {
    # it get passed a null df in the case of a failure in rma naive
  if (is.null(mi)) return(NULL)
  
  mi_sub <- mi
  mi_sub$draws <- mi_sub$draws[1:m_use, , drop = FALSE]
  return(mi_sub)
}

# Wrapper for univ 
safe_adj_uni <- function(mi, m_use, delta, method.re) {
  if (is.null(mi)) return(c(est = NA_real_, ess = NA_real_))

  tryCatch({
    mi_subset <- subset_mi(mi, m_use)
    res <- adj_univariate(mi_subset, delta = delta, method.re = method.re, track.ess = TRUE)
    return(c(est = res$Estimate[1], ess = res$ess[1]))
  }, error = function(e) {
    return(c(est = NA_real_, ess = NA_real_))
  })
}

# Wrapper for biv 
safe_adj_biv <- function(mi, m_use, delta) {
  if (is.null(mi)) return(c(est = NA_real_, ess = NA_real_, failed.prop = NA_real_))

  tryCatch({
    mi_subset <- subset_mi(mi, m_use)
    res <- adj_bivariate(mi_subset, delta = delta)
    return(c(est = res$Estimate[1], ess = res$ess[1], failed.prop = res$failed.proportion[1]))
  }, error = function(e) {
    return(c(est = NA_real_, ess = NA_real_, failed.prop = NA_real_))
  })
}

run_comprehensive_scenario <- function(K, tau2, delta, n_sim, true_theta, n_cores, M) {
  start <- Sys.time()
  
  results_list <- mclapply(1:n_sim, function(i) {
    tryCatch({
      # Generate full data
      full_data <- generate_bivariate_ma(K = K, theta = c(true_theta, true_theta), 
                                         rho_w = 0.4, tau2 = c(tau2, tau2))
      
      # --- Baseline (No ORB) - REML ---
      res_full_reml <- rma(yi = O1_yi, sei = O1_sei, data = full_data, method = "REML", 
                           control=list(stepadj=0.1, rel.tol=1e-8, maxiter=1000))
      est_full_reml <- as.numeric(res_full_reml$beta)

      # --- Baseline (No ORB) - PM ---
      res_full_pm <- rma(yi = O1_yi, sei = O1_sei, data = full_data, method = "PM", 
                         control=list(stepadj=0.1, rel.tol=1e-8, maxiter=1000))
      est_full_pm <- as.numeric(res_full_pm$beta)

      
      # Impose ORB 
      obs_data <- impose_orb(full_data, p1 = 0.4, delta_sim = delta, select_type = "zscore", orb.se = TRUE)
      
      # If no studies are missing, skip  
      if (all(!is.na(obs_data$O1_yi))) {
        return(data.frame(
        # specific simulation parameters
        K              = K,
        tau2           = tau2,
        delta          = delta,

        
        # Naive results
        full.reml      = est_full_reml, 
        full.pm        = est_full_pm,
        naive.uni.reml = NA, 
        naive.uni.pm   = NA,
        naive.uni_new_reml = NA,
        naive.uni_new_pm = NA, 
        naive.biv      = NA, 
        naive.biv_new = NA,
        
        # Estimates (50, 200, 1000)
        uni.reml_50       = NA,
        uni.reml_200      = NA,
        uni.reml_1000     = NA,
        
        uni.pm_50         = NA,
        uni.pm_200        = NA,
        uni.pm_1000       = NA,
        
        uni_new.reml_50   = NA,
        uni_new.reml_200  = NA,
        uni_new.reml_1000 = NA,
        
        uni_new.pm_50     = NA,
        uni_new.pm_200    = NA,
        uni_new.pm_1000   = NA,
        
        biv_50            = NA,
        biv_200           = NA,
        biv_1000          = NA,
        
        biv_new_50        = NA,
        biv_new_200       = NA,
        biv_new_1000      = NA,
        
        # ess 
        ess_uni.reml_50      = NA,
        ess_uni.reml_200     = NA,
        ess_uni.reml_1000     = NA,

        ess_uni.pm_50         = NA,
        ess_uni.pm_200        = NA,
        ess_uni.pm_1000       = NA,

        ess_uni_new.reml_50   = NA,
        ess_uni_new.reml_200  = NA,
        ess_uni_new.reml_1000 = NA,

        ess_uni_new.pm_50     = NA,
        ess_uni_new.pm_200    = NA,
        ess_uni_new.pm_1000   = NA,
        
        ess_biv_50        = NA,
        ess_biv_200       = NA,
        ess_biv_1000      = NA,
        
        ess_biv_new_50    = NA,
        ess_biv_new_200   = NA,
        ess_biv_new_1000  = NA,
        
        # Failure 
        failed_biv_1000     = NA,
        failed_biv_new_1000 = NA
      ))
      }

      # --- Impute SEs  ---
      obs_data <- impute_missing_se(obs_data, target_theta_col = "O1_yi", target_se_col = "O1_sei", n_col = "n_total")

      # ---------------------------------------------------------
      # 1. Univariate (Old) - REML
      # ---------------------------------------------------------
      mi_uni_reml <- run_univariate_imputation(obs_data, theta_col = "O1_yi", se_col = "O1_sei", new_version = FALSE, method.re = "REML", m = M)
      naive.uni.reml <- as.numeric(mi_uni_reml$res_naive$beta) 
      
      uni_reml_50   <- safe_adj_uni(mi_uni_reml, 50, delta, "REML")
      uni_reml_200  <- safe_adj_uni(mi_uni_reml, 200, delta, "REML")
      uni_reml_1000 <- safe_adj_uni(mi_uni_reml, 1000, delta, "REML")

      # ---------------------------------------------------------
      # 2. Univariate (Old) - PM
      # ---------------------------------------------------------
      mi_uni_pm <- run_univariate_imputation(obs_data, theta_col = "O1_yi", se_col = "O1_sei", new_version = FALSE, method.re = "PM", m = M)
      naive.uni.pm <- as.numeric(mi_uni_pm$res_naive$beta) 
      
      uni_pm_50   <- safe_adj_uni(mi_uni_pm, 50, delta, "PM")
      uni_pm_200  <- safe_adj_uni(mi_uni_pm, 200, delta, "PM")
      uni_pm_1000 <- safe_adj_uni(mi_uni_pm, 1000, delta, "PM")

      # ---------------------------------------------------------
      # 3. Univariate (New) - REML
      # ---------------------------------------------------------
      mi_uni_new_reml <- run_univariate_imputation(obs_data, theta_col = "O1_yi", se_col = "O1_sei", new_version = TRUE, method.re = "REML", m = M)
      naive.uni_new_reml <- as.numeric(mi_uni_new_reml$res_naive$beta) 

      uni_new_reml_50   <- safe_adj_uni(mi_uni_new_reml, 50, delta, "REML")
      uni_new_reml_200  <- safe_adj_uni(mi_uni_new_reml, 200, delta, "REML")
      uni_new_reml_1000 <- safe_adj_uni(mi_uni_new_reml, 1000, delta, "REML")

      # ---------------------------------------------------------
      # 4. Univariate (New) - PM
      # ---------------------------------------------------------
      mi_uni_new_pm <- run_univariate_imputation(obs_data, theta_col = "O1_yi", se_col = "O1_sei", new_version = TRUE, method.re = "PM", m = M)
      naive.uni_new_pm <- as.numeric(mi_uni_new_pm$res_naive$beta) 

      uni_new_pm_50   <- safe_adj_uni(mi_uni_new_pm, 50, delta, "PM")
      uni_new_pm_200  <- safe_adj_uni(mi_uni_new_pm, 200, delta, "PM")
      uni_new_pm_1000 <- safe_adj_uni(mi_uni_new_pm, 1000, delta, "PM")

      # ---------------------------------------------------------
      # 5. Bivariate (Old)
      # ---------------------------------------------------------
      mi_biv <- run_bivariate_imputation(obs_data, theta_cols = c("O1_yi", "O2_yi"), se_cols = c("O1_sei", "O2_sei"), new_version = FALSE, rho_w = 0.4, m = M)
      naive.biv <- as.numeric(mi_biv$res_naive$beta[1]) 
      
      biv_50   <- safe_adj_biv(mi_biv, 50, delta)
      biv_200  <- safe_adj_biv(mi_biv, 200, delta)
      biv_1000 <- safe_adj_biv(mi_biv, 1000, delta)

      # ---------------------------------------------------------
      # 6. Bivariate (New)
      # ---------------------------------------------------------
      mi_biv_new <- run_bivariate_imputation(obs_data, theta_cols = c("O1_yi", "O2_yi"), se_cols = c("O1_sei", "O2_sei"), new_version = TRUE, rho_w = 0.4, m = M)
      naive.biv_new <- as.numeric(mi_biv_new$res_naive$beta[1]) 
      

      biv_new_50   <- safe_adj_biv(mi_biv_new, 50, delta)
      biv_new_200  <- safe_adj_biv(mi_biv_new, 200, delta)
      biv_new_1000 <- safe_adj_biv(mi_biv_new, 1000, delta)

      # ---------------------------------------------------------
      # return just a row df
      # ---------------------------------------------------------
        return(data.frame(
        # specific simulation parameters
        K              = K,
        tau2           = tau2,
        delta          = delta,

        
        # Naive results
        full.reml      = est_full_reml, 
        full.pm        = est_full_pm,
        naive.uni.reml = naive.uni.reml, 
        naive.uni.pm   = naive.uni.pm,
        naive.uni_new_reml = naive.uni_new_reml,
        naive.uni_new_pm = naive.uni_new_pm, 
        naive.biv      = naive.biv, 
        naive.biv_new = naive.biv_new,
        
        # Estimates (50, 200, 1000)
        uni.reml_50       = uni_reml_50["est"],
        uni.reml_200      = uni_reml_200["est"],
        uni.reml_1000     = uni_reml_1000["est"],
        
        uni.pm_50         = uni_pm_50["est"],
        uni.pm_200        = uni_pm_200["est"],
        uni.pm_1000       = uni_pm_1000["est"],
        
        uni_new.reml_50   = uni_new_reml_50["est"],
        uni_new.reml_200  = uni_new_reml_200["est"],
        uni_new.reml_1000 = uni_new_reml_1000["est"],
        
        uni_new.pm_50     = uni_new_pm_50["est"],
        uni_new.pm_200    = uni_new_pm_200["est"],
        uni_new.pm_1000   = uni_new_pm_1000["est"],
        
        biv_50            = biv_50["est"],
        biv_200           = biv_200["est"],
        biv_1000          = biv_1000["est"],
        
        biv_new_50        = biv_new_50["est"],
        biv_new_200       = biv_new_200["est"],
        biv_new_1000      = biv_new_1000["est"],
        
        # ess 
        ess_uni.reml_50      = uni_reml_50["ess"],
        ess_uni.reml_200     = uni_reml_200["ess"],
        ess_uni.reml_1000     = uni_reml_1000["ess"],

        ess_uni.pm_50         = uni_pm_50["ess"],
        ess_uni.pm_200        = uni_pm_200["ess"],
        ess_uni.pm_1000       = uni_pm_1000["ess"],

        ess_uni_new.reml_50   = uni_new_reml_50["ess"],
        ess_uni_new.reml_200  = uni_new_reml_200["ess"],
        ess_uni_new.reml_1000 = uni_new_reml_1000["ess"],
        
        ess_uni_new.pm_50     = uni_new_pm_50["ess"],
        ess_uni_new.pm_200    = uni_new_pm_200["ess"],
        ess_uni_new.pm_1000   = uni_new_pm_1000["ess"],
        
        ess_biv_50        = biv_50["ess"],
        ess_biv_200       = biv_200["ess"],
        ess_biv_1000      = biv_1000["ess"],
        
        ess_biv_new_50    = biv_new_50["ess"],
        ess_biv_new_200   = biv_new_200["ess"],
        ess_biv_new_1000  = biv_new_1000["ess"],
        
        # Failure 
        failed_biv_1000     = biv_1000["failed.prop"],
        failed_biv_new_1000 = biv_new_1000["failed.prop"]
      ))
      
    }, error = function(e) {
      # Return NAs if the entire iteration fails
        return(data.frame(
        # specific simulation parameters
        K              = K,
        tau2           = tau2,
        delta          = delta,

        
        # Naive results
        full.reml      = NA, 
        full.pm        = NA,
        naive.uni.reml = NA, 
        naive.uni.pm   = NA,
        naive.uni_new_reml = NA,
        naive.uni_new_pm = NA, 
        naive.biv      = NA, 
        naive.biv_new = NA,
        
        # Estimates (50, 200, 1000)
        uni.reml_50       = NA,
        uni.reml_200      = NA,
        uni.reml_1000     = NA,
        
        uni.pm_50         = NA,
        uni.pm_200        = NA,
        uni.pm_1000       = NA,
        
        uni_new.reml_50   = NA,
        uni_new.reml_200  = NA,
        uni_new.reml_1000 = NA,
        
        uni_new.pm_50     = NA,
        uni_new.pm_200    = NA,
        uni_new.pm_1000   = NA,
        
        biv_50            = NA,
        biv_200           = NA,
        biv_1000          = NA,
        
        biv_new_50        = NA,
        biv_new_200       = NA,
        biv_new_1000      = NA,
        
        # ess 
        ess_uni.reml_50      = NA,
        ess_uni.reml_200     = NA,
        ess_uni.reml_1000     = NA,

        ess_uni.pm_50         = NA,
        ess_uni.pm_200        = NA,
        ess_uni.pm_1000       = NA,

        ess_uni_new.reml_50   = NA,
        ess_uni_new.reml_200  = NA,
        ess_uni_new.reml_1000 = NA,

        ess_uni_new.pm_50     = NA,
        ess_uni_new.pm_200    = NA,
        ess_uni_new.pm_1000   = NA,
        
        ess_biv_50        = NA,
        ess_biv_200       = NA,
        ess_biv_1000      = NA,
        
        ess_biv_new_50    = NA,
        ess_biv_new_200   = NA,
        ess_biv_new_1000  = NA,
        
        # Failure 
        failed_biv_1000     = NA,
        failed_biv_new_1000 = NA
      ))
    })
  }, mc.cores = n_cores, mc.preschedule = TRUE, mc.set.seed = TRUE)
  
  end <- Sys.time()
  
  # Combine the replicates into a single data frame
  results_df        <- do.call(rbind, results_list)
  results_df$time_min <- as.numeric(difftime(end, start, units = "mins"))
  
  return(results_df)
}


# NOTE: I realized that if the impose_ORB creates no missingess, I hadle it easy at the beginning, but 'full.reml' and 'full.pm' still will have values while the others will all be NULL.
# instead if impose_ORB creates too much missigness, i.e. it only leave one study or zero, then 'mi_uni_reml$res_naive$beta' will throw an error and it will fall back to the case where all columns are set to NA.
# basically the wrapper 'safe_adj' is not that useful since it will crash anyway. 

# ---------------------------------------------------------
# Loop Over the Grid
# ---------------------------------------------------------
all_results <- vector("list", nrow(scenarios))

for (s in 1:nrow(scenarios)) {
  cat(sprintf("\nRunning scenario %d/%d  (K = %d, tau2 = %.2f, delta = %.1f)...\n",
              s, nrow(scenarios), scenarios$K[s], scenarios$tau2[s], scenarios$delta[s]))
  
  set.seed(s) 
  all_results[[s]] <- run_comprehensive_scenario(K          = scenarios$K[s],
                                                 tau2       = scenarios$tau2[s],
                                                 delta      = scenarios$delta[s],
                                                 n_sim      = n_sim,
                                                 true_theta = true_theta,
                                                 n_cores    = n_cores,
                                                 M          = M) 
}

# Bind all scenarios together
all_results_df <- do.call(rbind, all_results)
saveRDS(all_results_df, file = "data/exploratory.check.rds")

cat("'all_results_df' saved in 'data/exploratory.check.rds'")




# ---------------------------------------------------------
# Summary Table
# ---------------------------------------------------------
summary_tbl <- do.call(rbind, lapply(1:nrow(scenarios), function(s) {
  # Subset to the current scenario
  sub  <- subset(all_results_df, K == scenarios$K[s] & tau2 == scenarios$tau2[s] & delta == scenarios$delta[s])
  
  # Bias calculations
  cols_to_avg <- c("full.reml", "full.pm", "naive.uni.reml", "naive.uni.pm", 
                   "naive.uni_new_reml", "naive.uni_new_pm", "naive.biv", "naive.biv_new", 
                   "uni.reml_50", "uni.reml_200", "uni.reml_1000",
                   "uni.pm_50", "uni.pm_200", "uni.pm_1000",
                   "uni_new.reml_50", "uni_new.reml_200", "uni_new.reml_1000",
                   "uni_new.pm_50", "uni_new.pm_200", "uni_new.pm_1000",
                   "biv_50", "biv_200", "biv_1000",
                   "biv_new_50", "biv_new_200", "biv_new_1000")
  
  bias <- colMeans(sub[, cols_to_avg], na.rm = TRUE) - true_theta
  
  # avg ess
  avg_ess_uni.reml_1000       <- mean(sub$ess_uni.reml_1000, na.rm = TRUE)
  avg_ess_uni.reml_50         <- mean(sub$ess_uni.reml_50, na.rm = TRUE)
  avg_ess_uni.reml_200        <- mean(sub$ess_uni.reml_200, na.rm = TRUE)
  avg_ess_uni.pm_1000         <- mean(sub$ess_uni.pm_1000, na.rm = TRUE)
  avg_ess_uni.pm_50           <- mean(sub$ess_uni.pm_50, na.rm = TRUE)
  avg_ess_uni.pm_200          <- mean(sub$ess_uni.pm_200, na.rm = TRUE)
  avg_ess_uni_new.reml_1000   <- mean(sub$ess_uni_new.reml_1000, na.rm = TRUE)
  avg_ess_uni_new.reml_50     <- mean(sub$ess_uni_new.reml_50, na.rm = TRUE)
  avg_ess_uni_new.reml_200    <- mean(sub$ess_uni_new.reml_200, na.rm = TRUE)
  avg_ess_uni_new.pm_1000     <- mean(sub$ess_uni_new.pm_1000, na.rm = TRUE)
  avg_ess_uni_new.pm_50       <- mean(sub$ess_uni_new.pm_50, na.rm = TRUE)
  avg_ess_uni_new.pm_200      <- mean(sub$ess_uni_new.pm_200, na.rm = TRUE)
  
  avg_ess_biv_50              <- mean(sub$ess_biv_50, na.rm = TRUE)
  avg_ess_biv_200             <- mean(sub$ess_biv_200, na.rm = TRUE)
  avg_ess_biv_1000            <- mean(sub$ess_biv_1000, na.rm = TRUE)
  
  avg_ess_biv_new_50          <- mean(sub$ess_biv_new_50, na.rm = TRUE)
  avg_ess_biv_new_200         <- mean(sub$ess_biv_new_200, na.rm = TRUE)
  avg_ess_biv_new_1000        <- mean(sub$ess_biv_new_1000, na.rm = TRUE)
  
  avg_failed_biv_1000         <- mean(sub$failed_biv_1000, na.rm = TRUE)
  avg_failed_biv_new_1000     <- mean(sub$failed_biv_new_1000, na.rm = TRUE)
  
  # build one row each scenario
  data.frame(K         = scenarios$K[s],
             tau2      = scenarios$tau2[s],
             delta     = scenarios$delta[s],
             t(round(bias, 5)),
             avg_ESS_Uni_REML_M50       = round(avg_ess_uni.reml_50, 1),
             avg_ESS_Uni_REML_M200      = round(avg_ess_uni.reml_200, 1),
             avg_ESS_Uni_REML_M1000     = round(avg_ess_uni.reml_1000, 1),
             avg_ESS_Uni_PM_M50         = round(avg_ess_uni.pm_50, 1),
             avg_ESS_Uni_PM_M200        = round(avg_ess_uni.pm_200, 1),
             avg_ESS_Uni_PM_M1000       = round(avg_ess_uni.pm_1000, 1),
             avg_ESS_Uni_New_REML_M50   = round(avg_ess_uni_new.reml_50, 1),
             avg_ESS_Uni_New_REML_M200  = round(avg_ess_uni_new.reml_200, 1),
             avg_ESS_Uni_New_REML_M1000 = round(avg_ess_uni_new.reml_1000, 1),
             avg_ESS_Uni_New_PM_M50     = round(avg_ess_uni_new.pm_50, 1),
             avg_ESS_Uni_New_PM_M200    = round(avg_ess_uni_new.pm_200, 1),
             avg_ESS_Uni_New_PM_M1000   = round(avg_ess_uni_new.pm_1000, 1),
             
             avg_ESS_Biv_M50            = round(avg_ess_biv_50, 1),
             avg_ESS_Biv_M200           = round(avg_ess_biv_200, 1),
             avg_ESS_Biv_M1000          = round(avg_ess_biv_1000, 1),
             
             avg_ESS_Biv_New_M50        = round(avg_ess_biv_new_50, 1),
             avg_ESS_Biv_New_M200       = round(avg_ess_biv_new_200, 1),
             avg_ESS_Biv_New_M1000      = round(avg_ess_biv_new_1000, 1),
             
             avg_Fail_Biv_M1000         = round(avg_failed_biv_1000, 4),
             avg_Fail_Biv_New_M1000     = round(avg_failed_biv_new_1000, 4),
             
             time_min  = round(sub$time_min[1], 2))
}))


print(summary_tbl, row.names = FALSE)
saveRDS(summary_tbl, file = "data/unbiasedness.comprehensive.rds")
cat ("SIMULATION completed!")


