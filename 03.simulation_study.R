source("00.functions.R")
set.seed(222)



# Generate bivariate data
cat("--- Generating Complete Data ---\n")
complete_data <- generate_bivariate_ma(K = 12, 
                                       theta = c(0.4, 0.4), 
                                       tau2 = c(0.06, 0.06), 
                                       rho_b = 0.4, 
                                       rho_w = 0.4, 
                                       n_arm = 50)
print(head(complete_data))


# Impose ORB
# 20% missingness on Outcome 1 with delta = 0.5
cat("\n--- Imposing ORB (Selection on z-score) ---\n")
orb_data <- impose_orb(data = complete_data, 
                       p1 = 0.2, 
                       delta_sim = 0.5, 
                       select_type = "zscore")
print(orb_data)





# --- Some cheeeck  ---

full_data <- generate_bivariate_ma(K = 25, theta = c(0.4, 0.4), tau2 = c(0.06, 0.06))

obs_data <- impose_orb(full_data, p1 = 0.4, delta_sim = 0, select_type = "zscore")


set.seed(1)
mi_uni <- run_univariate_imputation(obs_data, theta_col = "O1_yi", se_col = "O1_sei", new_version = TRUE)
mi_biv <- run_bivariate_imputation(obs_data, theta_cols = c("O1_yi", "O2_yi"), 
                                   se_cols = c("O1_sei", "O2_sei"), new_version = TRUE)

res_uni <- adj_univariate(mi_uni, delta = 0)
res_biv <- adj_bivariate(mi_biv, delta = 0)

res_uni$Estimate
res_biv$Estimate[1]




























































# ----------------------------------
# Univariate check for new version vs old version
# ----------------------------------

library(parallel)

RNGkind("L'Ecuyer-CMRG")
set.seed(1)

#RNGkind("Mersenne-Twister") then set.seed() to go back to normal


n_sim <- 1000
true_theta <- 0.4
scenarios <- expand.grid(
  tau2  = c(0.06, 0.36),
  delta = c(0, 0.7,1)
)

n_cores <- parallel::detectCores() - 1

run_scenario <- function(tau2, delta, n_sim, true_theta, n_cores) {

start <- Sys.time()

results_list <- mclapply(1:n_sim, function(i) {
  
  full_data <- generate_bivariate_ma(K = 25, theta = c(true_theta, true_theta), 
                                    rho_w = 0.4, tau2 = c(tau2, tau2))
  
  #   # this is the correct value we would get if there was no ORB
  res_full <- rma(yi = O1_yi, sei = O1_sei, data = full_data, method = "REML",
                  control=list(stepadj=0.1, rel.tol=1e-8, maxiter=1000))
  est_full <- as.numeric(res_full$beta)
  
  # PUT ORB
  obs_data <- impose_orb(full_data, p1 = 0.4, delta_sim = delta,  select_type = "zscore")
  
  # Univariate
  mi_uni <- run_univariate_imputation(obs_data, theta_col = "O1_yi", se_col = "O1_sei", new_version = FALSE)
  est_uni <- adj_univariate(mi_uni, delta = delta)$Estimate
  
  #naive univariate
  naive.uni <- as.numeric(mi_uni$res_naive$beta) 
  

  # New version
  mi_uni_new <- run_univariate_imputation(obs_data, theta_col = "O1_yi", se_col = "O1_sei", new_version = TRUE)
  est_uni_new <- adj_univariate(mi_uni_new, delta = delta)$Estimate
  

  
  # Return as a 1-row data frame
  return(data.frame(
    full = est_full, naive.uni = naive.uni, 
    uni = est_uni,
    uni_new = est_uni_new
  ))
  
}, mc.cores = n_cores, mc.preschedule = TRUE)
  end <- Sys.time()
  results_df        <- do.call(rbind, results_list)
  results_df$tau2   <- tau2
  results_df$delta  <- delta
  results_df$time_min <- as.numeric(difftime(end, start, units = "mins"))
  
  return(results_df)
}



# loop over scennarios
all_results <- vector("list", nrow(scenarios))

for (s in 1:nrow(scenarios)) {
  cat(sprintf("\nRunning scenario %d/%d  (tau2 = %.2f, delta = %.1f)...\n",
              s, nrow(scenarios), scenarios$tau2[s], scenarios$delta[s]))
  set.seed(s)
  all_results[[s]] <- run_scenario(tau2     = scenarios$tau2[s],
                                   delta    = scenarios$delta[s],
                                   n_sim    = n_sim,
                                   true_theta = true_theta,
                                   n_cores  = n_cores)
}

all_results <- do.call(rbind, all_results)


 #summ table
 summary_tbl <- do.call(rbind, lapply(1:nrow(scenarios), function(s) {
  sub  <- subset(all_results, tau2 == scenarios$tau2[s] & delta == scenarios$delta[s])
bias <- colMeans(sub[, c("full", "naive.uni", "uni", "uni_new")], na.rm = TRUE) - true_theta
  data.frame(tau2      = scenarios$tau2[s],
             delta     = scenarios$delta[s],
             t(round(bias, 5)),
             time_min  = round(sub$time_min[1], 2))
}))

print(summary_tbl, row.names = FALSE)
saveRDS(summary_tbl, file = "data/unbiasedness.uni.rds")


# RESULTS WITH p1 = 0.4: 

# tau2 delta     full    naive      uni  uni_new     time
# 0.06   0.0 -0.00346 -0.00254 -0.00264 -0.00259    56.49
# 0.36   0.0 -0.00019  0.00475  0.00479  0.00479    55.62
# 0.06   0.5 -0.00098  0.07204  0.02608 -0.01563    54.80
# 0.36   0.5  0.00048  0.25975  0.20300  0.04078    54.63
# 0.06   0.7 -0.00185  0.09278  0.03622 -0.00177    53.65
# 0.36   0.7 -0.00300  0.29470  0.22261  0.07577    57.43
# 0.06     1  0.00061  0.11086  0.04461  0.01424    60.06
# 0.36     1  0.00168  0.33691  0.25215  0.12960    57.79












































# ----------------------------------
# Unbiasedness check (both bivariate and univariate)
# ----------------------------------

RNGkind("L'Ecuyer-CMRG")
set.seed(1)

#RNGkind("Mersenne-Twister") then set.seed() to go back to normal


library(parallel)


n_sim <- 1000
true_theta <- 0.4
tau2 <- 0.36
delta = 0.5

n_cores <- parallel::detectCores() - 1

simulation_start <- Sys.time()

results_list <- mclapply(1:n_sim, function(i) {
  
  full_data <- generate_bivariate_ma(K = 25, theta = c(true_theta, true_theta), rho_w = 0.4, tau2 = c(tau2, tau2))
  
  #   # this is the correct value we would get if there was no ORB
  res_full <- rma(yi = O1_yi, sei = O1_sei, data = full_data, method = "REML")
  est_full <- as.numeric(res_full$beta)
  
  # PUT ORB
  obs_data <- impose_orb(full_data, p1 = 0.4, delta_sim = delta,  select_type = "zscore")
  
  # Univariate
  mi_uni <- run_univariate_imputation(obs_data, theta_col = "O1_yi", se_col = "O1_sei", new_version = FALSE)
  est_uni <- adj_univariate(mi_uni, delta = delta, theta_col = "O1_yi", se_col = "O1_sei")$Estimate
  
  #naive univariate
  naive.uni <- as.numeric(mi_uni$res_naive$beta) 


  # bivariate
  mi_biv <- run_bivariate_imputation(obs_data, theta_cols = c("O1_yi", "O2_yi"), 
                                     se_cols = c("O1_sei", "O2_sei"), new_version = FALSE, rho_w = 0.4)
  est_biv <- adj_bivariate(mi_biv, delta = delta)$Estimate[1]

  #naive biv
  naive.biv <- as.numeric(mi_biv$res_naive$beta[1]) 
  

  # New version
  mi_uni_new <- run_univariate_imputation(obs_data, theta_col = "O1_yi", se_col = "O1_sei", new_version = TRUE)
  est_uni_new <- adj_univariate(mi_uni_new, delta = delta, theta_col = "O1_yi", se_col = "O1_sei")$Estimate
  

  mi_biv_new <- run_bivariate_imputation(obs_data, theta_cols = c("O1_yi", "O2_yi"), 
                                         se_cols = c("O1_sei", "O2_sei"), new_version = TRUE, rho_w = 0.4)
  est_biv_new <- adj_bivariate(mi_biv_new, delta = delta)$Estimate[1]
  
  # Return as a 1-row data frame
  return(data.frame(
    full = est_full, naive.uni = naive.uni, 
    naive.biv = naive.biv, 
    uni = est_uni, biv = est_biv,
    uni_new = est_uni_new, biv_new = est_biv_new
  ))
  
}, mc.cores = n_cores, mc.preschedule = TRUE)

# Combine
results_df <- do.call(rbind, results_list)

simulation_end <- Sys.time()

bias_results <- colMeans(results_df, na.rm = TRUE) - true_theta
print(round(bias_results, 3))


# see the time
difftime(simulation_end, simulation_start, units = "mins"), "minutes\n")

saveRDS(bias_results, file = "data/unbiasedness.biv.rds")










































































# ---------------------------
# Let's see the time consumption 
# ---------------------------


library(profvis)


object <- profvis({
  true_theta <- 0.4
  tau2 <- 0.36
  delta <- 0.5
  
  full_data <- generate_bivariate_ma(K = 25, theta = c(true_theta, true_theta), rho_w = 0.4, tau2 = c(tau2, tau2))
  
  obs_data <- impose_orb(full_data, p1 = 0.4, delta_sim = delta,  select_type = "zscore")
  
  # Univariate
  mi_uni <- run_univariate_imputation(obs_data, theta_col = "O1_yi", se_col = "O1_sei", new_version = FALSE)
  est_uni <- adj_univariate(mi_uni, delta = delta)$Estimate
  
  # Bivariate 
  mi_biv <- run_bivariate_imputation(obs_data, theta_cols = c("O1_yi", "O2_yi"), 
                                     se_cols = c("O1_sei", "O2_sei"), new_version = FALSE, rho_w = 0.4)
  est_biv <- adj_bivariate(mi_biv, delta = delta)$Estimate[1]
})

#adj_univariate --> 4.3 secs
# adj_bivariate --> 9.0 secs

# so each one of the replicates take 15 secs (with the 1000 IS)

# but we want to do 1900 replications --> 15 seconds $\times$ 1900 replicates =  a lot




# RESULTS USING PM:

# PM cannot be used for rma.mv, try ML instead

object_PM <- profvis({
  true_theta <- 0.4
  tau2 <- 0.36
  delta <- 0.5
  
  full_data <- generate_bivariate_ma(K = 25, theta = c(true_theta, true_theta), rho_w = 0.4, tau2 = c(tau2, tau2))
  
  obs_data <- impose_orb(full_data, p1 = 0.4, delta_sim = delta,  select_type = "zscore")
  
  # Univariate
  mi_uni <- run_univariate_imputation(obs_data, theta_col = "O1_yi", se_col = "O1_sei", new_version = FALSE, method.re = "PM")
  est_uni <- adj_univariate(mi_uni, delta = delta,  method.re = "PM")$Estimate
  
  # Bivariate 
  mi_biv <- run_bivariate_imputation(obs_data, theta_cols = c("O1_yi", "O2_yi"), 
                                     se_cols = c("O1_sei", "O2_sei"), new_version = FALSE, rho_w = 0.4,  method.re = "ML")
  est_biv <- adj_bivariate(mi_biv, delta = delta,  method.re = "ML")$Estimate[1]
})


#adj_univariate --> 5.9 secs
# adj_bivariate --> 9.3 secs








library(microbenchmark)

set.seed(123)
N <- 2000

dat <- data.frame(
  study = factor(sample(1:200, N, replace = TRUE)),
  yi = rnorm(N),
  vi = runif(N, 0.05, 0.2)
)


bench <- microbenchmark(

  rma_REML = rma(yi, vi, data = dat, method = "REML"),

  rma_PM = rma(yi, vi, data = dat, method = "PM"),

  rma_mv_REML = rma.mv(
    yi, vi,
    random = ~ 1 | study,
    data = dat,
    method = "REML"
  ),

  rma_mv_ML = rma.mv(
    yi, vi,
    random = ~ 1 | study,
    data = dat,
    method = "ML"
  ),

  times = 50
)

print(bench)



























# ----------------------------------
# Final to send to the cluster:
# ----------------------------------


library(parallel)

scenarios_grid <- expand.grid(
  K = c(6, 12, 25),                      
  p1 = c(0.2, 0.4),  # paper not clear                  
  tau2_val = c(0, 0.02, 0.06, 0.36),     
  theta_1 = c(0, 0.4, 0.8), 
  theta_2 =  c(0, 0.4, 0.8),          
  rho_b = c(0, 0.4, 1),                  
  rho_w = c(0, 0.4),                     
  delta_sim = seq(0, 1, by = 0.2),       
  delta_est = seq(0, 1, by = 0.2),       
  select_type = c("zscore", "effect")
  )




function_to_create <- function(scenario){ 
  K <- scenario$K
  p1 <- scenario$p1
   (...) # maybe using with? 

  # generate data 
  full_data <- generate_bivariate_ma(K = K, 
                                    theta = c(theta_1, theta_2), 
                                    tau2 = c()
                                    rho_b = rho_b,
                                    rho_w = rho_w
                                    )
}






n_cores <- parallel::detectCores() - 1

final_results_list <- mclapply(
  X = 1:nrow(scenarios_grid), 
  FUN = function_to_create, 
  mc.cores = n_cores, 
  mc.preschedule = TRUE 
)

# iN THE END I want a 31,104-row summary
final_simulation_metrics <- do.call(rbind, final_results_list)