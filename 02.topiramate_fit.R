source("00.functions.R")
df_topiramate <- readRDS("data/df_topiramate.RDS")



#quick test
biv_mi_results <- run_bivariate_imputation(
  data = df_topiramate, 
  theta_cols = c("log_RR1", "log_RR2"), 
  se_cols = c("se_RR1", "se_RR2"), 
  rho_w = 0.7, 
  m = 50, new_version = FALSE 
)


set.seed(1)
biv_adj_zscore <- adj_bivariate(
  mi_results = biv_mi_results, 
  delta = 0.5, 
  select_type = "zscore")



set.seed(1)















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
    a_eff <- adj_univariate(mi, delta = 0.5, select_type = "effect")
    a_z   <- adj_univariate(mi, delta = 0.5, select_type = "zscore")
    
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

table_2_uni_repro <- final_uni_df %>%
  dplyr::mutate(
    Outcome = ifelse(Outcome == "O1", "Outcome 1", "Outcome 2"),
    
        Est_SE = paste0(sprintf("%.2f", Estimate), " (", sprintf("%.2f", SE), ")"),
    CI_95  = paste0(sprintf("%.2f", CI_Lower), " to ", sprintf("%.2f", CI_Upper))
  ) %>%
  dplyr::select(Correlation, Measure, Approach, Outcome, Est_SE, CI_95) %>%
  tidyr::pivot_wider(
    names_from = Outcome, 
    values_from = c(Est_SE, CI_95),
    names_glue = "{Outcome}_{.value}"
  ) %>%
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
    a_eff <- adj_bivariate(mi_biv, delta = 0.5, select_type = "effect")
    a_z   <- adj_bivariate(mi_biv, delta = 0.5, select_type = "zscore")
    
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

table_2_biv_repro <- final_biv_df %>%
  dplyr::mutate(
    Outcome = ifelse(Outcome == "O1", "Outcome 1", "Outcome 2"),
    
    Est_SE = paste0(sprintf("%.2f", Estimate), " (", sprintf("%.2f", SE), ")"),
    CI_95  = paste0(sprintf("%.2f", CI_Lower), " to ", sprintf("%.2f", CI_Upper))
  ) %>%
  dplyr::select(Correlation, Measure, Approach, Outcome, Est_SE, CI_95) %>%
  tidyr::pivot_wider(
    names_from = Outcome, 
    values_from = c(Est_SE, CI_95),
    names_glue = "{Outcome}_{.value}"
  ) %>%
  dplyr::select(Correlation, Measure, Approach, 
                `Outcome 1_Est_SE`, `Outcome 1_CI_95`, 
                `Outcome 2_Est_SE`, `Outcome 2_CI_95`)

# View the final reproducible Bivariate block of Table 2
View(table_2_biv_repro)
















# ----------------------------------
# FINAL TABLE 2
# ----------------------------------
FULL_TABLE_2 <- rbind(table_2_uni_repro, table_2_biv_repro)
# View(FULL_TABLE_2)

saveRDS(FULL_TABLE_2, file = "data/Table_2.rds")




















# Next step is erplicate Section 3.2 in which you have rho that is estimated in various way. 
# Do this after the simulation study, when you are sure your method works, since now you dont really have a
# baseline



































# Now I need to create Figure 1 if possible. 
# on the x axis delta values that range from 0 to 1.6 with a 0.2 step 
# y axis the overall summary estimate

# one colum for log OR and one for log RR 
# the rows are selection on estimate and selection on Z score 
library(ggplot2)
library(dplyr)

# Define the grid
deltas <- seq(0, 1.6, by = 0.2)
measures <- c("OR", "RR")
select_types <- c("effect", "zscore")

# List to store the plotted data
plot_data_list <- list()

# Figure 1 specifically looks at Outcome 2 ("Seizure Freedom")
# because it has the highest missingness, showing the clearest ORB curve.
out_num <- 2 

for (meas in measures) {
  
  cat(sprintf("Running Initial Imputation for Measure = %s...\n", meas))
  
  theta_name <- paste0("log_", meas, out_num)
  se_name <- paste0("se_", meas, out_num)
  
  # Step 1: Run Imputation ONCE per measure.
  mi_results <- run_univariate_imputation(
    data = df_topiramate, 
    theta_col = theta_name, 
    se_col = se_name, 
    m = 1000,
    new_version = TRUE 
  )
  
  # Extract Naive baseline (Constant across all deltas)
  naive_est <- as.numeric(mi_results$res_naive$beta)
  naive_ci_lb <- mi_results$res_naive$ci.lb
  naive_ci_ub <- mi_results$res_naive$ci.ub
  
  # Step 2: Loop the Importance Sampling over Delta and Selection Type
  for (sel in select_types) {
    for (d in deltas) {
      
      # Calculate the adjusted estimate for this specific delta
      adj_res <- adj_univariate(
        mi_results = mi_results, 
        delta = d, 
        select_type = sel
      )
      
      # Store in our dataframe
      plot_data_list[[length(plot_data_list) + 1]] <- data.frame(
        Measure = meas,
        Selection = ifelse(sel == "effect", "Selection on log OR/log RR", "Selection on Z-score"),
        Delta = d,
        Adj_Est = adj_res$Estimate,
        Adj_CI_Lower = adj_res$CI_Lower,
        Adj_CI_Upper = adj_res$CI_Upper,
        Naive_Est = naive_est,
        Naive_CI_Lower = naive_ci_lb,
        Naive_CI_Upper = naive_ci_ub
      )
    }
  }
}

# Bind into a single dataframe for ggplot
df_fig1 <- do.call(rbind, plot_data_list)


















# Format the facet labels to exactly match the paper
df_fig1$Measure <- factor(df_fig1$Measure, levels = c("OR", "RR"), 
                          labels = c("Effect Estimate: Odds Ratio (OR)", 
                                     "Effect Estimate: Relative Risk (RR)"))

df_fig1$Selection <- factor(df_fig1$Selection, 
                            levels = c("Selection on log OR/log RR", "Selection on Z-score"))

# Generate the Plot
fig_1 <- ggplot(df_fig1, aes(x = Delta)) +
  
  # 1. Naive Estimate (Red dashed line and light red shaded CI)
  geom_ribbon(aes(ymin = Naive_CI_Lower, ymax = Naive_CI_Upper), fill = "red", alpha = 0.1) +
  geom_hline(aes(yintercept = Naive_Est), color = "red", linetype = "dashed", linewidth = 1) +
  
  # 2. Adjusted Estimate (Blue solid line with points and light blue shaded CI)
  geom_ribbon(aes(ymin = Adj_CI_Lower, ymax = Adj_CI_Upper), fill = "blue", alpha = 0.2) +
  geom_line(aes(y = Adj_Est), color = "blue", linewidth = 1) +
  geom_point(aes(y = Adj_Est), color = "blue", size = 2) +
  
  # 3. Facet Grid (Rows = Selection Type, Columns = Measure)
  facet_grid(Selection ~ Measure) +
  
  # 4. Aesthetics and Labels
  labs(
    x = expression(delta ~ "(Selection Weight)"),
    y = "Overall Summary Estimate"
  ) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "gray95", color = "black"),
    strip.text = element_text(size = 11, face = "bold"),
    axis.title = element_text(size = 12, face = "bold"),
    panel.grid.minor = element_blank()
  )

# View the plot
print(fig_1)









































library(ggplot2)
library(dplyr)

# -------------------------------------------------------------------------
# FIGURE 3
# -------------------------------------------------------------------------

#  grid
deltas <- seq(0, 1.6, by = 0.2)
measures <- c("OR", "RR")
select_types <- c("effect", "zscore")


out_num <- 2      # We only plot Seizure Freedom
rho_fig3 <- 0.5   # The paper uses r = 0.5  or 0.7?? 

plot_data_list_fig3 <- list()

for (meas in measures) {
  
  cat(sprintf("Running Imputations for Measure = %s...\n", meas))
  
  # Column names
  theta_name_uni <- paste0("log_", meas, out_num)
  se_name_uni <- paste0("se_", meas, out_num)
  
  theta_cols_biv <- paste0("log_", meas, c(1, 2))
  se_cols_biv <- paste0("se_", meas, c(1, 2))
  
  # run both imp
  mi_uni <- run_univariate_imputation(df_topiramate, theta_name_uni, se_name_uni, m = 1000, new_version = TRUE)
  mi_biv <- run_bivariate_imputation(df_topiramate, theta_cols_biv, se_cols_biv, rho_w = rho_fig3, m = 1000, new_version = TRUE)

  # Extract the naive   
  # Univariate beta is length 1
  naive_uni_est <- as.numeric(mi_uni$res_naive$beta) 
  
  # Bivariate beta is length 2 (Outcome 1, Outcome 2) --> We want the second one.
  naive_biv_est <- as.numeric(mi_biv$res_naive$beta)[2] 
  
  # loop the adj for each effect and delta
  for (sel in select_types) {
    for (d in deltas) {
      
      # Adjust Univariate
      adj_u <- adj_univariate(mi_uni, delta = d, select_type = sel)
      
      # Adjust Bivariate
      adj_b <- adj_bivariate(mi_biv, delta = d, select_type = sel)
      
      # Extract ONLY Outcome 2 from Bivariate results
      adj_b_O2 <- adj_b[adj_b$Outcome == "O2", ]
      
      # Format Selection Label
      sel_label <- ifelse(sel == "effect", "Selection on log OR/log RR", "Selection on Z-score")
      
      # Append Univariate Row
      plot_data_list_fig3[[length(plot_data_list_fig3) + 1]] <- data.frame(
        Measure = meas,
        Selection = sel_label,
        Delta = d,
        Model = "Univariate",
        Adj_Est = adj_u$Estimate,
        Adj_CI_Lower = adj_u$CI_Lower,
        Adj_CI_Upper = adj_u$CI_Upper,
        Naive_Est = naive_uni_est
      )
      
      # Append Bivariate Row
      plot_data_list_fig3[[length(plot_data_list_fig3) + 1]] <- data.frame(
        Measure = meas,
        Selection = sel_label,
        Delta = d,
        Model = "Bivariate",
        Adj_Est = adj_b_O2$Estimate,
        Adj_CI_Lower = adj_b_O2$CI_Lower,
        Adj_CI_Upper = adj_b_O2$CI_Upper,
        Naive_Est = naive_biv_est

      )
    }
  }
  
  print("One done")
}

df_fig3 <- do.call(rbind, plot_data_list_fig3)




df_fig3$Measure <- factor(df_fig3$Measure, levels = c("OR", "RR"), 
                          labels = c("Effect Estimate: Odds Ratio (OR)", 
                                     "Effect Estimate: Relative Risk (RR)"))

df_fig3$Selection <- factor(df_fig3$Selection, 
                            levels = c("Selection on log OR/log RR", "Selection on Z-score"))

# Generate the Plot
fig_3 <- ggplot(df_fig3, aes(x = Delta, group = Model, color = Model, fill = Model)) +
  
  # naive
  geom_hline(aes(yintercept = Naive_Est, color = Model), linetype = "dashed", linewidth = 1) +
  
  # put CI with geom_ribbon
  geom_ribbon(aes(ymin = Adj_CI_Lower, ymax = Adj_CI_Upper), alpha = 0.2, color = NA) +
  
  # adj est
  geom_line(aes(y = Adj_Est), linewidth = 1) +
  geom_point(aes(y = Adj_Est), size = 2) +
  
  facet_grid(Selection ~ Measure) +

  scale_color_manual(values = c("Univariate" = "mediumblue", "Bivariate" = "seagreen")) +
  scale_fill_manual(values = c("Univariate" = "mediumblue", "Bivariate" = "seagreen")) +
  
  labs(
    x = expression(delta ~ "(Selection Weight)"),
    y = "Overall Summary Measure",
    color = "Model Type",
    fill  = "Model Type"
  ) +
  theme_bw() 

print(fig_3)

