# ---------------------------------------------------------
# Comprehensive Unbiasedness Check 
# ---------------------------------------------------------



source("00.functions.R")
library(parallel)


RNGkind("L'Ecuyer-CMRG")
set.seed(1)

# Simulation parameters
n_sim <- 1000
true_theta <- 0.4
n_cores <- 25
M <- 200 # imputations

# Define the grid of scenarios
scenarios <- expand.grid(
  K     = c(25, 6),
  tau2  = c(0.06, 0.36),
  delta = c(0, 0.7, 1)
)

# ---------------------------------------------------------
#  Define the Scenario Function
# ---------------------------------------------------------
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
      obs_data <- impose_orb(full_data, p1 = 0.4, delta_sim = delta,  select_type = "zscore")
      

      # If no studies were suppressed, skip imputation entirely
      
      if (all(!is.na(obs_data$O1_yi))) {
        return(data.frame(
          K = K,
          full.reml=est_full_reml, full.pm=est_full_pm,
          naive.uni.reml=NA, naive.uni.pm=NA, naive.biv=NA,
          uni.reml=NA, uni.pm=NA, uni_new.reml=NA, uni_new.pm=NA,
          biv=NA, biv_new=NA
        ))
      }


      # --- Univariate (Old) - REML ---
      mi_uni_reml <- run_univariate_imputation(obs_data, theta_col = "O1_yi", se_col = "O1_sei", new_version = FALSE, method.re = "REML", m = M)
      est_uni_reml <- adj_univariate(mi_uni_reml, delta = delta, method.re = "REML")$Estimate
      naive.uni.reml <- as.numeric(mi_uni_reml$res_naive$beta) 

      # --- Univariate (Old) - PM ---
      mi_uni_pm <- run_univariate_imputation(obs_data, theta_col = "O1_yi", se_col = "O1_sei", new_version = FALSE, method.re = "PM", m = M)
      est_uni_pm <- adj_univariate(mi_uni_pm, delta = delta, method.re = "PM")$Estimate
      naive.uni.pm <- as.numeric(mi_uni_pm$res_naive$beta) 
      
      # --- Univariate (New) - REML ---
      mi_uni_new_reml <- run_univariate_imputation(obs_data, theta_col = "O1_yi", se_col = "O1_sei", new_version = TRUE, method.re = "REML", m = M)
      est_uni_new_reml <- adj_univariate(mi_uni_new_reml, delta = delta, method.re = "REML")$Estimate
      
      # --- Univariate (New) - PM ---
      mi_uni_new_pm <- run_univariate_imputation(obs_data, theta_col = "O1_yi", se_col = "O1_sei", new_version = TRUE, method.re = "PM", m = M)
      est_uni_new_pm <- adj_univariate(mi_uni_new_pm, delta = delta, method.re = "PM")$Estimate

      # --- Bivariate  ---
      mi_biv <- run_bivariate_imputation(obs_data, theta_cols = c("O1_yi", "O2_yi"), 
                                         se_cols = c("O1_sei", "O2_sei"), new_version = FALSE, rho_w = 0.4, m = M)
      est_biv <- adj_bivariate(mi_biv, delta = delta)$Estimate[1]
      naive.biv <- as.numeric(mi_biv$res_naive$beta[1]) 
      
      # --- Bivariate (New) ---
      mi_biv_new <- run_bivariate_imputation(obs_data, theta_cols = c("O1_yi", "O2_yi"), 
                                             se_cols = c("O1_sei", "O2_sei"), new_version = TRUE, rho_w = 0.4, m = M)
      est_biv_new <- adj_bivariate(mi_biv_new, delta = delta)$Estimate[1]
      
      # Return as a 1-row data frame
      return(data.frame(
        K              = K,
        full.reml      = est_full_reml, 
        full.pm        = est_full_pm,
        naive.uni.reml = naive.uni.reml, 
        naive.uni.pm   = naive.uni.pm,
        naive.biv      = naive.biv, 
        uni.reml       = est_uni_reml, 
        uni.pm         = est_uni_pm,
        uni_new.reml   = est_uni_new_reml, 
        uni_new.pm     = est_uni_new_pm,
        biv            = est_biv,
        biv_new        = est_biv_new
      ))
      
    }, error = function(e) {
      # Return NAs if convergence fails to avoid crashing the job
      return(data.frame(
        K = K,
        full.reml=NA, full.pm=NA, naive.uni.reml=NA, naive.uni.pm=NA, naive.biv=NA, 
        uni.reml=NA, uni.pm=NA, uni_new.reml=NA, uni_new.pm=NA, biv=NA, biv_new=NA
      ))
    })
  }, mc.cores = n_cores, mc.preschedule = TRUE, mc.set.seed = TRUE)
  
  end <- Sys.time()
  
  # Combine the replicates into a single data frame
  results_df        <- do.call(rbind, results_list)
  results_df$tau2   <- tau2
  results_df$delta  <- delta
  results_df$time_min <- as.numeric(difftime(end, start, units = "mins"))
  
  return(results_df)
}

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

# ---------------------------------------------------------
# Create the Summary Table
# ---------------------------------------------------------
summary_tbl <- do.call(rbind, lapply(1:nrow(scenarios), function(s) {
  # Subset to the current scenario
  sub  <- subset(all_results_df, K == scenarios$K[s] & tau2 == scenarios$tau2[s] & delta == scenarios$delta[s])
  
  # Calculate bias for all 11 estimators
  cols_to_avg <- c("full.reml", "full.pm", "naive.uni.reml", "naive.uni.pm", "naive.biv", 
                   "uni.reml", "uni.pm", "uni_new.reml", "uni_new.pm", 
                   "biv", "biv_new")
  bias <- colMeans(sub[, cols_to_avg], na.rm = TRUE) - true_theta
  
  # Build the row
  data.frame(K         = scenarios$K[s],
             tau2      = scenarios$tau2[s],
             delta     = scenarios$delta[s],
             t(round(bias, 5)),
             time_min  = round(sub$time_min[1], 2))
}))

# Print and save
print(summary_tbl, row.names = FALSE)
saveRDS(summary_tbl, file = "data/unbiasedness.comprehensive.rds")

















#################################
# ------------------------- 
# Results 
# -------------------------


res <- readRDS("data/unbiasedness.comprehensive.rds")

library(ggplot2)
library(tidyr)
library(dplyr)

method_cols <- c(
    "biv_new", "biv", "uni_new.pm", "uni_new.reml", "uni.pm", "uni.reml", "naive.biv", "naive.uni.pm", "naive.uni.reml",  "full.pm", "full.reml"
)

plot_df <- res %>%
  pivot_longer(cols = all_of(method_cols), names_to = "method", values_to = "bias") %>%
  mutate(
    scenario = paste0("K=", K, ", tau2=", tau2, ", delta=", delta),
    method = factor(method, levels = method_cols)
  )

p1 <- ggplot(plot_df, aes(x = method, y = bias, color = method, group = 1)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  geom_line(linewidth = 0.5, show.legend = FALSE) +
  geom_point(size = 2, show.legend = FALSE) +
  facet_wrap(~scenario, ncol = 3) +
  coord_flip() +
  labs(
    title = "Bias by Scenario",
    x = "Method",
    y = "Bias"
  ) +
  theme_bw(base_size = 11)

print(p1)


# put larger correlation p_W