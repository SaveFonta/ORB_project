library(metafor)
library(parallel)

# Single simulation replicate - how long does it take?
source("00.functions.R")

set.seed(1)
full_data <- generate_bivariate_ma(K = 25, theta = c(0.4, 0.4), rho_w = 0.4, tau2 = c(0.06, 0.06))
obs_data  <- impose_orb(full_data, p1 = 0.4, delta_sim = 0.7, select_type = "zscore")

# Time each component separately
# t1 <- system.time({
#   mi_uni <- run_univariate_imputation(obs_data, theta_col = "O1_yi", se_col = "O1_sei", 
#                                       new_version = FALSE, method.re = "REML", m = 1000)
# })

# t2 <- system.time({
#   adj_univariate(mi_uni, delta = 0.7, method.re = "REML")
# })

# t3 <- system.time({
#   mi_biv <- run_bivariate_imputation(obs_data, theta_cols = c("O1_yi", "O2_yi"),
#                                      se_cols = c("O1_sei", "O2_sei"), 
#                                      new_version = FALSE, rho_w = 0.4, m = 1000)
# })

# t4 <- system.time({
#   adj_bivariate(mi_biv, delta = 0.7)
# } )

# cat("run_univariate_imputation: ", t1["elapsed"], "s\n")
# cat("adj_univariate:            ", t2["elapsed"], "s\n")
# cat("run_bivariate_imputation:  ", t3["elapsed"], "s\n")
# cat("adj_bivariate:             ", t4["elapsed"], "s\n")
# cat("Total per replicate:       ", sum(t1["elapsed"],t2["elapsed"],t3["elapsed"],t4["elapsed"]), "s\n")
# cat("Estimated total (1000 sim x 6 scenarios / n_cores): ",
   #  sum(t1["elapsed"],t2["elapsed"],t3["elapsed"],t4["elapsed"]) * 1000 * 6 / 20 / 60, "mins\n")


source("00.functions.R")

# --------------------------------------------------------
# wrapper: for errors to continue
# --------------------------------------------------------
safe_adj <- function(mi, m_use, delta) {
  if (is.null(mi)) return(NA_real_)
  tryCatch(
    adj_bivariate(within(mi, draws <- draws[1:m_use, ]), delta = delta)$Estimate[1],
    error = function(e) NA_real_
  )
}

# --------------------------------------------------------
# Core function: run M sensitivity for one scenario
# --------------------------------------------------------
run_M_sensitivity <- function(K, tau2, delta_sim, seeds = c(1, 2, 3)) {

  cat(sprintf("\nK = %d | tau2 = %.2f | delta = %.1f\n", K, tau2, delta_sim))

  rows <- lapply(seeds, function(seed) {
    set.seed(seed)

    full_data <- generate_bivariate_ma(K = K, theta = c(0.4, 0.4),
                                       rho_w = 0.4, tau2 = c(tau2, tau2))
    obs_data  <- impose_orb(full_data, p1 = 0.4, delta_sim = delta_sim,
                            select_type = "zscore")
    mi_biv    <- run_bivariate_imputation(obs_data,
                                          theta_cols = c("O1_yi", "O2_yi"),
                                          se_cols    = c("O1_sei", "O2_sei"),
                                          rho_w = 0.4, m = 1000)

    r50   <- safe_adj(mi_biv, m_use = 50,   delta = delta_sim)
    r200  <- safe_adj(mi_biv, m_use = 200,  delta = delta_sim)
    r1000 <- safe_adj(mi_biv, m_use = 1000, delta = delta_sim)

    cat(sprintf("  Seed %d  |  M=50: %.4f  |  M=200: %.4f  |  M=1000: %.4f\n",
                seed, r50, r200, r1000))

    data.frame(K = K, tau2 = tau2, delta = delta_sim,
               seed = seed, M50 = r50, M200 = r200, M1000 = r1000)
  })

  do.call(rbind, rows)
}

# --------------------------------------------------------
# Grid of scenarios
# --------------------------------------------------------
scenarios <- expand.grid(
  K          = c(6, 25),
  tau2       = c(0, 0.36),
  delta_sim  = c(0, 1),
  stringsAsFactors = FALSE
)

seeds <- c(1, 2, 3)

# --------------------------------------------------------
# Run all scenarios
# --------------------------------------------------------
all_results <- vector("list", nrow(scenarios))

for (s in seq_len(nrow(scenarios))) {
  all_results[[s]] <- run_M_sensitivity(
    K         = scenarios$K[s],
    tau2      = scenarios$tau2[s],
    delta_sim = scenarios$delta_sim[s],
    seeds     = seeds
  )
}

results_df <- do.call(rbind, all_results)
rownames(results_df) <- NULL

# --------------------------------------------------------
# Save and print summary
# --------------------------------------------------------
saveRDS(results_df, file = "data/M_sensitivity.rds")

cat("\n\n===== FULL RESULTS TABLE =====\n")
print(results_df)