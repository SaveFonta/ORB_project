# ---------------------------------------------------------
# FINAL script to send to cluster
# ---------------------------------------------------------

source("00.functions.R")
library(parallel)

RNGkind("L'Ecuyer-CMRG")
set.seed(1)

# ---------------------------------------------------------
# Global Parameters
# ---------------------------------------------------------
n_sim         <- 1900   # (maybe also 1000 can work? Dunno)
M_imputations <- 1000   
n_cores       <- 25    

# ---------------------------------------------------------
# Grid, can easily change to evaluate more scenarios, this is taken from Table 3
# ---------------------------------------------------------
scenarios_grid <- expand.grid(
  K           = c(6, 12, 25),                      
  p1          = c(0, 0.2, 0.4),       # Should we actually include p1 = 0?? I wouldn't
  tau2_val    = c(0, 0.02, 0.06, 0.36),     
  theta_1     = c(0, 0.4, 0.8), 
  theta_2     = c(0, 0.4, 0.8),          
  rho_b       = c(0, 0.4, 1),                  
  rho_w       = c(0, 0.4),                     
  delta_sim   = seq(0, 1, by = 0.2),       
  delta_est   = seq(0, 1, by = 0.2),  
  select_type = c("zscore", "effect"),
  stringsAsFactors = FALSE
)

# ---------------------------------------------------------
# Wrappers I wrote to handle two cases of errors: 
# 1) -> failure of convergence of rma inside imputation function (either for numerical issues or more probably cause there are less than 2 reported studies)
# 2) -> failure of convergence of rma inside adjust function for numerical issues 

# both ways it invalidates the whole draw also for the pther methods. Maybe can change this logic? Not sure
# ---------------------------------------------------------

safe_adj_uni <- function(mi, delta, sel_type) {

  # it get passed a null df in the case of a failure in rma naive
  if (is.null(mi)) return(c(est = NA, ci_l = NA, ci_u = NA, ess = NA))

  tryCatch({
    res <- adj_univariate(mi, delta = delta, select_type = sel_type, method.re = "REML", track.ess = TRUE)
    return(c(est = res$Estimate[1], ci_l = res$CI_Lower[1], ci_u = res$CI_Upper[1], ess = res$ess[1]))
  }, error = function(e) c(est = NA, ci_l = NA, ci_u = NA, ess = NA))
}

safe_adj_biv <- function(mi, delta, sel_type) {
  if (is.null(mi)) return(c(est = NA, ci_l = NA, ci_u = NA, ess = NA, fail = NA))
  
  tryCatch({
    res <- adj_bivariate(mi, delta = delta, select_type = sel_type)
    return(c(est = res$Estimate[1], ci_l = res$CI_Lower[1], ci_u = res$CI_Upper[1], 
             ess = res$ess[1], fail = res$failed.proportion[1]))
  }, error = function(e) c(est = NA, ci_l = NA, ci_u = NA, ess = NA, fail = NA))
}
# ---------------------------------------------------------
# MAIN FUNCTION
# ---------------------------------------------------------
function_to_create <- function(scenario_idx) { 
  s           <- scenarios_grid[scenario_idx, ]
  true_theta  <- s$theta_1
  
  #Extract specifics of this scenarios

  # params for generating 
  K <- s$K
  theta_1 <- s$theta_1
  theta_2 <-  s$theta_2
  tau2_val <- s$tau2_val
  rho_b <- s$rho_b
  rho_w <- s$rho_w
  delta_sim <- s$delta_sim

  #params for orb
  p1 <- s$p1
  select_type <- s$select_type

  # params for the estimation 
  delta_est <- s$delta_est
  # rho_w
  # select_type
  

  # Run n_sim replicates for this specific row of the grid
  replicates <- lapply(1:n_sim, function(i) {
    tryCatch({
      # Generate complete data
      full_data <- generate_bivariate_ma(K = K, theta = c(theta_1, theta_2), 
                                         tau2 = c(tau2_val, tau2_val),
                                         rho_b = rho_b, rho_w = rho_w)
      
      # Baseline, full model (basically p_1 = 0)
      res_full <- rma(yi = O1_yi, sei = O1_sei, data = full_data, method = "REML")
      est_full <- as.numeric(res_full$beta)
      se_full  <- as.numeric(res_full$se)
      
      # Impose ORB
      obs_data <- impose_orb(full_data, p1 = p1, delta_sim = delta_sim, 
                             select_type = select_type, orb.se = TRUE)
      
      # by chance, ORB sometimes doesnt create any missing outcomes (more common when K = 6), in that case we finish here this simulation returning the full model....
      # I don't know another way to handle this, cause we could even use some kind of repeat block in the case of no missingness.
      # But this would bias a bit the results I think ? Since we would condition the distribution of the draw only on those draws that produces missingness? 

      if (all(!is.na(obs_data$O1_yi))) return(NULL)

      
      # Impute Se
      obs_data <- impute_missing_se(obs_data, target_theta_col = "O1_yi", target_se_col = "O1_sei", n_col = "n_total")
      
      # fit the functions
      mi_uni  <- run_univariate_imputation(obs_data, theta_col = "O1_yi", se_col = "O1_sei", 
                                           new_version = TRUE, method.re = "REML", m = M_imputations)
      adj_uni <- safe_adj_uni(mi_uni, delta = delta_est, sel_type = select_type)
      
      mi_biv  <- run_bivariate_imputation(obs_data, theta_cols = c("O1_yi", "O2_yi"), 
                                          se_cols = c("O1_sei", "O2_sei"), 
                                          new_version = TRUE, rho_w = rho_w, m = M_imputations)
      adj_biv <- safe_adj_biv(mi_biv, delta = delta_est, sel_type = select_type)

      naive_uni_est <- if (is.null(mi_uni)) NA_real_ else as.numeric(mi_uni$res_naive$beta)
      naive_biv_est <- if (is.null(mi_biv)) NA_real_ else as.numeric(mi_biv$res_naive$beta[1])
      
      return(data.frame(full = est_full, 
                        naive_uni = naive_uni_est,
                        naive_biv = naive_biv_est,
                        uni = adj_uni["est"], 
                        biv = adj_biv["est"],
                        uni_ci_l = adj_uni["ci_l"], 
                        uni_ci_u = adj_uni["ci_u"],
                        biv_ci_l = adj_biv["ci_l"], 
                        biv_ci_u = adj_biv["ci_u"],
                        u_ess = adj_uni["ess"], 
                        b_ess = adj_biv["ess"], 
                        b_f = adj_biv["fail"]))
    }, error = function(e) return(NULL))
  })
  
res_df <- do.call(rbind, replicates)
  if (is.null(res_df) || nrow(res_df) == 0) return(NULL)

  data.frame(
    scenario_idx  = scenario_idx,
    K             = K,
    theta_1       = theta_1,
    theta_2       = theta_2,
    tau2_val      = tau2_val,
    rho_b         = rho_b,
    rho_w         = rho_w,
    delta_sim     = delta_sim,
    p1            = p1,
    select_type   = select_type,
    delta_est     = delta_est,
    
    N_Successful   = nrow(res_df), 
    
    Bias_Naive_Uni = mean(res_df$naive_uni, na.rm=TRUE) - true_theta,
    Bias_Naive_Biv = mean(res_df$naive_biv, na.rm=TRUE) - true_theta,
    
    Bias_Adj_Uni   = mean(res_df$uni, na.rm=TRUE) - true_theta,
    MSE_Adj_Uni    = mean((res_df$uni - true_theta)^2, na.rm=TRUE),
    Coverage_Uni   = mean(res_df$uni_ci_l <= true_theta & res_df$uni_ci_u >= true_theta, na.rm=TRUE),
    CI_Width_Uni   = mean(res_df$uni_ci_u - res_df$uni_ci_l, na.rm=TRUE),
    Avg_ESS_Uni    = mean(res_df$u_ess, na.rm=TRUE),

    Bias_Adj_Biv   = mean(res_df$biv, na.rm=TRUE) - true_theta,
    MSE_Adj_Biv    = mean((res_df$biv - true_theta)^2, na.rm=TRUE),
    Coverage_Biv   = mean(res_df$biv_ci_l <= true_theta & res_df$biv_ci_u >= true_theta, na.rm=TRUE),
    CI_Width_Biv   = mean(res_df$biv_ci_u - res_df$biv_ci_l, na.rm=TRUE),
    Avg_ESS_Biv    = mean(res_df$b_ess, na.rm=TRUE),
    
    Fail_Rate_Biv  = mean(res_df$b_f, na.rm=TRUE)
  )
}

final_results_list <- mclapply(
  X = 1:nrow(scenarios_grid), 
  FUN = function_to_create, 
  mc.cores = n_cores, 
  mc.preschedule = TRUE,
  mc.set.seed  = TRUE
)

final_simulation_metrics <- do.call(rbind, final_results_list)
saveRDS(final_simulation_metrics, "ORB_Full_Simulation.rds")