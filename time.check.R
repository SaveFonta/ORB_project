library(metafor)
library(parallel)

# Single simulation replicate - how long does it take?
source("00.functions.R")

RNGkind("L'Ecuyer-CMRG")
set.seed(1)




# --------------------------------------------------------
# wrapper: for errors to continue
# --------------------------------------------------------
safe_adj <- function(mi, m_use, delta) {
  if (is.null(mi)) return(NA_real_)

  mi_subset$draws <- mi_subset$draws[1:m_use, , drop = FALSE]

  tryCatch({
    mi_subset <- mi
    mi_subset$draws <- mi_subset$draws[1:m_use, , drop = FALSE]
    
    # Extract the Estimate for Outcome 1
    res <- adj_bivariate(mi_subset, delta = delta, select_type = "zscore")
    return(res$Estimate[1])
  }, error = function(e) return(NA_real_))
}

# --------------------------------------------------------
# Core function: run M sensitivity for one scenario
# --------------------------------------------------------
run_M_sensitivity <- function(K, tau2, delta_sim, seeds = c(1, 2, 3)) {

  cat(sprintf("\nK = %d | tau2 = %.2f | delta = %.1f\n", K, tau2, delta_sim))

  rows <- lapply(seeds, function(seed) {
    set.seed(seed)

    # generate data
    full_data <- generate_bivariate_ma(K = K, theta = c(0.4, 0.4),
                                       rho_w = 0.4, tau2 = c(tau2, tau2))

    # impose orb
    obs_data  <- impose_orb(full_data, p1 = 0.4, delta_sim = delta_sim,
                            select_type = "zscore" , orb.se = TRUE)
 
    # impute se
    obs_data <- impute_missing_se(data = obs_data, 
                                  target_theta_col = "O1_yi", 
                                  target_se_col = "O1_sei", 
                                  n_col = "n_total")

    # bivariate
    mi_biv    <- run_bivariate_imputation(obs_data,
                                          theta_cols = c("O1_yi", "O2_yi"),
                                          se_cols    = c("O1_sei", "O2_sei"),
                                          rho_w = 0.4, m = 1000)

    # m
    r50   <- safe_adj(mi_biv, m_use = 50,   delta = delta_sim)
    r200  <- safe_adj(mi_biv, m_use = 200,  delta = delta_sim)
    r1000 <- safe_adj(mi_biv, m_use = 1000, delta = delta_sim)


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

library(dplyr)

obj <- readRDS("data/M_sensitivity.rds")

RES <- obj  |> 
  mutate (
    diff.50 = abs(M50 - M1000), 
    diff.200 = abs(M200 - M1000)
   )

round(RES,3)

mean(RES$diff.50, na.rm = TRUE)
mean(RES$diff.200, na.rm = TRUE)
library(ggplot2)
library(tidyr)

long <- obj  |> 
    pivot_longer(cols = c("M50", "M200", "M1000"),
                 names_to = "M_value",
                 values_to = "estimate") |>
    mutate(scenario = paste("K =", K, "| tau2 =", tau2, "| delta =", delta))

# Create faceted plot
ggplot(long, aes(x = M_value, y = estimate, color = as.factor(seed), group = seed)) +
  geom_point(size = 2) +
  geom_line() +
  facet_wrap(~scenario) +
  labs(title = "How much does M affect the results?",
       x = "Number of Imputations",
       y = "Estimate",
       color = "Seed") +
  theme_bw()
