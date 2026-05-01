source("00.functions.R")
df_topiramate <- readRDS("data/df_topiramate.RDS")




set.seed(1)
measures <- c("OR", "RR")
rhos <- c(0, 0.7) # that needs to be used for multivariate
m_imputations <- 1000
results_list_uni <- list()




# -------------------------------------------------------------------------
# Univariate Table Generation 
# -------------------------------------------------------------------------


#loop the two measures
for (meas in measures) {

  #loop the two outcomes
  for (out_num in 1:2) {
    
    cat(sprintf("Running Univariate Imputation: Measure = %s, Outcome = %s\n", meas, out_num))
    
    theta_name <- paste0("log_", meas, out_num)
    se_name <- paste0("se_", meas, out_num)
    
    # imputation
    mi <- run_univariate_imputation(
      data = df_topiramate, 
      theta_col = theta_name, 
      se_col = se_name, 
      m = m_imputations,
      new_version = FALSE
    )
    
    #  Extract Naive into a df matching adj_univariate output
    naive_df <- data.frame(
      Approach = "Naive estimate",
      Estimate = as.numeric(mi$res_naive$beta),
      SE       = mi$res_naive$se,
      CI_Lower = mi$res_naive$ci.lb,
      CI_Upper = mi$res_naive$ci.ub
    )
    
    #  Adj - Effect and Z-score
    a_eff <- adj_univariate(mi, delta = 0.5, select_type = "effect", track.ess = FALSE)
    a_z   <- adj_univariate(mi, delta = 0.5, select_type = "zscore", track.ess = FALSE)
    
    # Stack 
    tmp <- rbind(naive_df, a_eff, a_z)
    
    # add info columns (correlation is empyt)
    tmp$Correlation <- "" 
    tmp$Measure     <- paste0("log ", meas)
    tmp$Outcome     <- paste0("O", out_num) # e.g., "O1" or "O2"
    
    results_list_uni[[length(results_list_uni) + 1]] <- tmp
  }
}

final_uni_df <- do.call(rbind, results_list_uni)


# -------------------------------------------------------------------------
# Formatting 
# -------------------------------------------------------------------------

table_2_uni_repro <- final_uni_df  |> 
  dplyr::mutate(
    Outcome = ifelse(Outcome == "O1", "Outcome 1", "Outcome 2"),
    
        Est_SE = paste0(sprintf("%.2f", Estimate), " (", sprintf("%.2f", SE), ")"),
    CI_95  = paste0(sprintf("%.2f", CI_Lower), " to ", sprintf("%.2f", CI_Upper))
  )  |> 
  dplyr::select(Correlation, Measure, Approach, Outcome, Est_SE, CI_95)  |> 
  tidyr::pivot_wider(
    names_from = Outcome, 
    values_from = c(Est_SE, CI_95),
    names_glue = "{Outcome}_{.value}"
  )  |> 
  dplyr::select(Correlation, Measure, Approach, 
                `Outcome 1_Est_SE`, `Outcome 1_CI_95`, 
                `Outcome 2_Est_SE`, `Outcome 2_CI_95`)

# Change the Correlation column to say "Univariate" so it looks nice in the final table
table_2_uni_repro$Correlation <- "Univariate"






































# -------------------------------------------------------------------------
# Bivariate Table Generation (r = 0 and r = 0.7)
# -------------------------------------------------------------------------

results_list_biv <- list()


#loop over the rhos 
for (rho in rhos) {

  #loop over the two measures (here we dont loop for the outcomes since they are processed together)
  for (m in measures) {
    
    cat(sprintf("Running Bivariate Imputation: Measure = %s, rho_w = %s\n", m, rho))
    
    theta_cols <- paste0("log_", m, c(1, 2))
    se_cols <- paste0("se_", m, c(1, 2))
    
    # Run imputation
    mi_biv <- run_bivariate_imputation(
      data = df_topiramate, 
      theta_cols = theta_cols, 
      se_cols = se_cols, 
      rho_w = rho, 
      m = m_imputations,
      new_version = FALSE
    )
    
    # create naive df similar to the one produced by adj_bivariate
    naive_df <- data.frame(
          Outcome   = c("O1", "O2"),
          Approach  = "Naive estimate",
          Estimate  = as.numeric(mi_biv$res_naive$beta),
          SE        = sqrt(diag(mi_biv$res_naive$vb)),
          CI_Lower  = mi_biv$res_naive$ci.lb,
          CI_Upper  = mi_biv$res_naive$ci.ub
        )
    
    # Calculate both Adjustments (each return a 2-row dataframe for O1 and O2)
    a_eff <- adj_bivariate(mi_biv, delta = 0.5, select_type = "effect", track.ess = FALSE, track.failed.proportion = FALSE)
    a_z   <- adj_bivariate(mi_biv, delta = 0.5, select_type = "zscore", track.ess = FALSE, track.failed.proportion = FALSE)
    
    # stack them
    tmp <- rbind(naive_df, a_eff, a_z)

    #add explaining columns
    tmp$Correlation <- paste0("r = ", rho)
    tmp$Measure     <- paste0("log ", m)
    
    results_list_biv[[length(results_list_biv) + 1]] <- tmp
  }
}

final_biv_df <- do.call(rbind, results_list_biv)

# -------------------------------------------------------------------------
# Formatting
# -------------------------------------------------------------------------
library(dplyr)
library(tidyr)

table_2_biv_repro <- final_biv_df  |> 
  dplyr::mutate(
    Outcome = ifelse(Outcome == "O1", "Outcome 1", "Outcome 2"),
    
    Est_SE = paste0(sprintf("%.2f", Estimate), " (", sprintf("%.2f", SE), ")"),
    CI_95  = paste0(sprintf("%.2f", CI_Lower), " to ", sprintf("%.2f", CI_Upper))
  )  |> 
  dplyr::select(Correlation, Measure, Approach, Outcome, Est_SE, CI_95)  |> 
  tidyr::pivot_wider(
    names_from = Outcome, 
    values_from = c(Est_SE, CI_95),
    names_glue = "{Outcome}_{.value}"
  )  |> 
  dplyr::select(Correlation, Measure, Approach, 
                `Outcome 1_Est_SE`, `Outcome 1_CI_95`, 
                `Outcome 2_Est_SE`, `Outcome 2_CI_95`)

















# ----------------------------------
# FINAL TABLE 2
# ----------------------------------
FULL_TABLE_2 <- rbind(table_2_uni_repro, table_2_biv_repro)

saveRDS(FULL_TABLE_2, file = "data/Table_2.rds")
cat("Table 2 created!!")


# FULL_TABLE_2 <- readRDS("data/Table_2.rds")
# View(FULL_TABLE_2)

















# Next step is erplicate Section 3.2 in which you have rho that is estimated in various way. 
# Do this after the simulation study, when you are sure your method works, since now you dont really have a
# baseline

