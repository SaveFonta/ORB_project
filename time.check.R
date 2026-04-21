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




mi_biv <- run_bivariate_imputation(obs_data, theta_cols = c("O1_yi", "O2_yi"),
                                   se_cols = c("O1_sei", "O2_sei"), 
                                   new_version = FALSE, rho_w = 0.4, m = 1000)
# draw once with m=1000, then just use the first 50/200/1000 rows

t5 <- system.time({
  res_m50  <- adj_bivariate(mi_biv_50  <- within(mi_biv, draws <- draws[1:50,  ]), delta = 0.7)
})

t6 <- system.time({
  res_m200 <- adj_bivariate(mi_biv_200 <- within(mi_biv, draws <- draws[1:200, ]), delta = 0.7)
})

t7 <- system.time({
  res_m1000 <- adj_bivariate(mi_biv, delta = 0.7)
})

print(rbind(res_m50$Estimate, res_m200$Estimate, res_m1000$Estimate))
cat("adj_bivariate M=50:  ", t5["elapsed"], "s\n")
cat("adj_bivariate M=200: ", t6["elapsed"], "s\n")
cat("adj_bivariate M=1000:", t7["elapsed"], "s\n")

