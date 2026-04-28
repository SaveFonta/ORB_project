#########################################
#   SE imputing check
#######################################

source("00.functions.R")
library(parallel)


RNGkind("L'Ecuyer-CMRG")
set.seed(1)

# Simulation parameters
n_sim <- 1000
true_theta <- 0.4
n_cores <- 25
M <- 200 # imputations




run_orb_se_comparison_scenario <- function(tau2, delta, orb_se, n_sim, true_theta, n_cores, M) {
  start <- Sys.time()
  
  results_list <- list()
  while (length(results_list) < n_sim) {
    needed <- n_sim - length(results_list)
    batch <- mclapply(seq_len(needed), function(i) {
      tryCatch({
        full_data <- generate_bivariate_ma(
          K = 25,
          theta = c(true_theta, true_theta),
          rho_w = 0.4,
          tau2 = c(tau2, tau2)
        )
        
        obs_data <- impose_orb(
          full_data,
          p1 = 0.4,
          delta_sim = delta,
          select_type = "zscore",
          orb.se = orb_se
        )
        
        if (all(!is.na(obs_data$O1_yi))) {
          return(NULL)
        }
        
        if (orb_se) {
          obs_data <- impute_missing_se(obs_data, target_theta_col = "O1_yi", target_se_col = "O1_sei")
        }
        
        mi_uni_new <- run_univariate_imputation(
          obs_data,
          theta_col = "O1_yi",
          se_col = "O1_sei",
          new_version = TRUE,
          method.re = "REML",
          m = M
        )
        est_uni_new <- adj_univariate(mi_uni_new, delta = delta, method.re = "REML")$Estimate
        naive.uni.reml <- as.numeric(mi_uni_new$res_naive$beta)
        
        mi_biv_new <- run_bivariate_imputation(
          obs_data,
          theta_cols = c("O1_yi", "O2_yi"),
          se_cols = c("O1_sei", "O2_sei"),
          new_version = TRUE,
          rho_w = 0.4,
          m = M
        )
        est_biv_new <- adj_bivariate(mi_biv_new, delta = delta)$Estimate[1]
        naive.biv <- as.numeric(mi_biv_new$res_naive$beta[1])
        
        data.frame(
          naive.uni.reml = naive.uni.reml,
          uni_new.reml = est_uni_new,
          naive.biv = naive.biv,
          biv_new = est_biv_new
        )
      }, error = function(e) NULL)
    }, mc.cores = n_cores, mc.preschedule = TRUE, mc.set.seed = TRUE)
    
    batch <- Filter(Negate(is.null), batch)
    if (length(batch) > 0) {
      results_list <- c(results_list, batch)
    }
  }
  
  end <- Sys.time()
  results_df <- do.call(rbind, results_list[seq_len(n_sim)])
  results_df$tau2 <- tau2
  results_df$delta <- delta
  results_df$orb_se <- orb_se
  results_df$time_min <- as.numeric(difftime(end, start, units = "mins"))
  
  results_df
}

orb_se_scenarios <- expand.grid(
  tau2 = c(0.06, 0.36),
  delta = c(0, 0.7, 1),
  orb_se = c(FALSE, TRUE)
)

orb_se_results <- vector("list", nrow(orb_se_scenarios))

for (s in 1:nrow(orb_se_scenarios)) {
  cat(sprintf(
    "\nRunning ORB-SE scenario %d/%d  (tau2 = %.2f, delta = %.1f, orb_se = %s)...\n",
    s,
    nrow(orb_se_scenarios),
    orb_se_scenarios$tau2[s],
    orb_se_scenarios$delta[s],
    orb_se_scenarios$orb_se[s]
  ))

  orb_se_results[[s]] <- run_orb_se_comparison_scenario(
    tau2 = orb_se_scenarios$tau2[s],
    delta = orb_se_scenarios$delta[s],
    orb_se = orb_se_scenarios$orb_se[s],
    n_sim = n_sim,
    true_theta = true_theta,
    n_cores = n_cores,
    M = M
  )
}

orb_se_results_df <- do.call(rbind, orb_se_results)













# ---------------
# summary
#-----------------
orb_se_summary_tbl <- do.call(rbind, lapply(1:nrow(orb_se_scenarios), function(s) {
  sub <- subset(
    orb_se_results_df,
    tau2 == orb_se_scenarios$tau2[s] &
      delta == orb_se_scenarios$delta[s] &
      orb_se == orb_se_scenarios$orb_se[s]
  )
  
  bias <- colMeans(sub[, c("naive.uni.reml", "uni_new.reml", "naive.biv", "biv_new")], na.rm = TRUE) - true_theta
  
  data.frame(
    tau2 = orb_se_scenarios$tau2[s],
    delta = orb_se_scenarios$delta[s],
    orb_se = orb_se_scenarios$orb_se[s],
    t(round(bias, 5)),
    time_min = round(sub$time_min[1], 2)
  )
}))

print(orb_se_summary_tbl, row.names = FALSE)
saveRDS(orb_se_summary_tbl, file = "data/unbiasedness.orb_se.comparison.rds")

