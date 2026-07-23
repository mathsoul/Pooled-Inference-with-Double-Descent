library(tidyverse)
library(Matrix)
library(MCMCpack)
library(doParallel)
library(doRNG)
library(data.table)

pkg_names <- c("tidyverse", "nlshrink", "MASS", "quadprog")
source("Code/func.R")

# Load helper functions and packages needed for the experiment.
n_cores = detectCores() - 1
cl <- makeCluster(n_cores)
registerDoParallel(cl)

# Retrieve all function names from the global environment in preparation for parallel computation 
all_objects <- ls(envir = .GlobalEnv)
all_functions <- all_objects[sapply(mget(all_objects, envir = .GlobalEnv), is.function)]

# Data Frequency ------------------------------------------------------------
data_freq = "Monthly"

# Load forecast, truth, and error files for the selected frequency.
f_data = fread(paste0("CleanedData/M4ScaledData/", data_freq, "_forecast.csv")) %>% dplyr::select(-V1)
true_data = fread(paste0("CleanedData/M4ScaledData/", data_freq, "_truth.csv")) %>% dplyr::select(-V1)
u_data = fread(paste0("CleanedData/M4ScaledData/", data_freq, "_err.csv")) %>% dplyr::select(-V1)

rank_vec = 1:8
methods = c("EW","Linear")
n_methods = length(methods)

group_size_vec = 4^(1:4)
ratio_mat = NULL

# Number of random shuffle simulations for robustness evaluation.
B = 100

for(i in 1:length(rank_vec)){
  rank_idx = rank_vec[i]
  time_series_names = getTimeSeriesNames(data_freq, rank_idx)
  if (length(time_series_names) < group_size_vec[length(group_size_vec)]) break
  
  n_prods = length(time_series_names)
  n_experts = 17
  n_periods = ncol(true_data) - 1
  
  # --- Keep 'id' for tracking ---
  true_comb_df = true_data %>% filter(id %in% time_series_names) %>% arrange(as.numeric(gsub("\\D", "", id)))
  f_comb_df = f_data %>% filter(id %in% time_series_names) %>% arrange(as.numeric(gsub("\\D", "", id)))
  u_comb_df = u_data %>% filter(id %in% time_series_names) %>% arrange(as.numeric(gsub("\\D", "", id)))
  
  # Save ordered IDs for aligning EW and PCR
  ordered_ids <- true_comb_df$id
  
  true_comb = as.data.frame(true_comb_df %>% dplyr::select(-id))
  f_comb = as.data.frame(f_comb_df %>% dplyr::select(-id, -group_id))
  u_comb = as.data.frame(u_comb_df %>% dplyr::select(-id, -group_id))
  
  # EW and S:linear -------------------------------------------------------------------
  L2_method_list = foreach(p = 1:n_prods, .options.RNG = 20240422,
                       .packages = pkg_names, .export = all_functions) %dorng% {
                         single_prod_idx = 1:n_experts + (p-1) * n_experts
                         true_sep = true_comb[p,]
                         f_sep = tibble(f_comb[single_prod_idx,]) 
                         u_sep = u_comb[single_prod_idx,]
                         u_method = suppressMessages(getSepUMethods(1:n_experts, methods, n = n_experts, f_sep, u_sep, true_sep))
                         apply(u_method, 1, function(x) mean(x^2))
                       }
  
  L2_method_mat = matrix(unlist(L2_method_list), nrow = n_prods, byrow = TRUE)
  colnames(L2_method_mat) = methods
  method_df <- data.frame(id = ordered_ids, L2_method_mat)
  
  # P:Ridgeless ----------------------------------------------------------
  Y_mat = t(apply(u_comb, 2, function(x) sumEachProd(x, n_prods, n_experts)/n_experts))
  Gamma0 = getGamma0(n_experts)
  Gamma = Matrix(diag(n_prods), sparse = TRUE) %x% Gamma0
  X_mat = as.matrix(t(u_comb) %*% Gamma)
  ridge0_test_err_std = suppressMessages(map_dfc(1:nrow(X_mat), ~ getRidge0TestErr(X_mat[-.x, , drop = FALSE], Y_mat[-.x, , drop = FALSE], X_mat[.x,,drop = FALSE],
                                                                                   Y_mat[.x,,drop = FALSE], stdize = TRUE)))
  L2_ridge0_std = apply(ridge0_test_err_std, 1, function(x) mean(x^2))
  P_Ridge_df <- data.frame(id = ordered_ids, P_Ridge = L2_ridge0_std)
  
  
  # S:Ridgeless, namely P:Ridgeless with M = 1 ------------------------------------------------------
  L2_R1_list <- foreach(p = 1:n_prods, .packages = c("tidyverse", "Matrix", "MASS", "quadprog"), .export = all_functions) %dorng% {
    single_prod_idx = 1:n_experts + (p-1) * n_experts
    u_sep = u_comb[single_prod_idx,]
    sub_n_prods = 1
    
    Y_mat <- t(apply(u_sep, 2, function(x) sumEachProd(x, sub_n_prods, n_experts) / n_experts))
    Y_mat <- matrix(Y_mat, nrow = n_periods, ncol = sub_n_prods, byrow = FALSE)

    Gamma0 <- getGamma0(n_experts)
    Gamma <- Matrix(diag(sub_n_prods), sparse = TRUE) %x% Gamma0
    X_mat <- as.matrix(t(u_sep) %*% Gamma)
    
    R1_test_err_std <- suppressMessages(
      map_dfc(1:nrow(X_mat), ~ getRidge0TestErr(
        X_mat[-.x, , drop = FALSE], Y_mat[-.x, , drop = FALSE], 
        X_mat[.x, , drop = FALSE], Y_mat[.x, , drop = FALSE], stdize = TRUE
      ))
    )
    
    L2_R1_std <- apply(R1_test_err_std, 1, function(x) mean(x^2))
    
    data.frame(id = ordered_ids[p], mse = L2_R1_std)
  }
  
  R1_df <- do.call(rbind, L2_R1_list) %>% arrange(as.numeric(gsub("\\D", "", id)))
  
  # Store the standardized error matrix from each shuffle simulation.
  sim_results_list <- list()
  
  # Repeat the experiment over many random shuffles to average out sampling variation.
  for(sim in 1:B) {
    # Unique seed for each iteration to guarantee reproducibility
    set.seed(2023 + sim) 
    shuffled_names = sample(time_series_names)
    # Group:Ridgeless ----------------------------------------------------------
    L2_ridge0_list_of_mats <- list()
    
    for (t in 1:length(group_size_vec)){
      group_size = group_size_vec[t] 
      
      n_groups = ceiling(length(shuffled_names) / group_size)
      groups = split(shuffled_names, rep(1:n_groups, each=group_size, length.out=length(shuffled_names)))
      
      export_functions <- c("getRidge0TestErr", "getGamma0", "sumEachProd", "replaceZeroWOne", "getTimeSeriesNames")
      
      group_results <- foreach(j = 1:n_groups, 
                               .packages = c("tidyverse", "Matrix", "MASS", "quadprog"),
                               .export = export_functions) %dorng% {
                                 
                                 sub_group <- groups[[j]]
                                 sub_n_prods <- length(sub_group)
                                 
                                 # Keep sub_ids to tracking
                                 sub_true_df <- true_data %>% filter(id %in% sub_group) %>% arrange(as.numeric(gsub("\\D", "", id)))
                                                                     sub_ids <- sub_true_df$id
                                                                     
                                                                    # true_comb_sub <- as.data.frame(sub_true_df %>% dplyr::select(-id))
                                                                    # f_comb_sub <- as.data.frame(f_data %>% filter(id %in% sub_group) %>% arrange(as.numeric(gsub("\\D", "", id))) %>% dplyr::select(-id, -group_id))
                                                                     u_comb_sub <- as.data.frame(u_data %>% filter(id %in% sub_group) %>% arrange(as.numeric(gsub("\\D", "", id))) %>% dplyr::select(-id, -group_id))
                                                                     
                                                                     Y_mat <- t(apply(u_comb_sub, 2, function(x) sumEachProd(x, sub_n_prods, n_experts) / n_experts))
                                                                     if(sub_n_prods == 1) {
                                                                       Y_mat <- matrix(Y_mat, nrow = n_periods, ncol = sub_n_prods, byrow = FALSE)
                                                                     }
                                                                     
                                                                     Gamma0 <- getGamma0(n_experts)
                                                                     Gamma <- Matrix(diag(sub_n_prods), sparse = TRUE) %x% Gamma0
                                                                     X_mat <- as.matrix(t(u_comb_sub) %*% Gamma)
                                                                     
                                                                     Gridge0_test_err_std <- suppressMessages(
                                                                       map_dfc(1:nrow(X_mat), ~ getRidge0TestErr(
                                                                         X_mat[-.x, , drop = FALSE], Y_mat[-.x, , drop = FALSE], 
                                                                         X_mat[.x, , drop = FALSE], Y_mat[.x, , drop = FALSE], stdize = TRUE
                                                                       ))
                                                                     )
                                                                     
                                                                     L2_Gridge0_std <- apply(Gridge0_test_err_std, 1, function(x) mean(x^2))
                                                                     
                                                                     data.frame(id = sub_ids, mse = L2_Gridge0_std)
                               }
      
      # Combine chunks and sort immediately by explicit ID
      col_name <-  paste0("R", group_size_vec[t])
      ridge_t_df <- do.call(rbind, group_results) %>% arrange(as.numeric(gsub("\\D", "", id)))
      
      if(t == 1) {
        L2_ridge0_df_combined <- ridge_t_df
        colnames(L2_ridge0_df_combined)[2] <- col_name
      } else {
        L2_ridge0_df_combined[[col_name]] <- ridge_t_df$mse
      }
    }
    
    # --- Merge Everything Safely By ID ---
    combined_sim_errors <- method_df %>% 
      inner_join(R1_df, by = "id") %>%
      inner_join(L2_ridge0_df_combined, by = "id") %>%
      inner_join(P_Ridge_df, by = "id") %>%
      arrange(as.numeric(gsub("\\D", "", id))) # Guaranteed fixed alignment across simulations
    
    sim_results_list[[sim]] <- as.matrix(combined_sim_errors %>% dplyr::select(-id))
  }
  
  # Average the results across all random shuffles.
  mean_error_mat <- reduce(sim_results_list, `+`) / B
  
  # Compute loss ratios using the aggregated average values.
  ratio_mat_temp = mean_error_mat / mean_error_mat[, ncol(mean_error_mat)]
  
  ratio_mat = rbind(ratio_mat, 
                    data.frame(ratio_mat_temp, dataset = paste0(data_freq, rank_idx)))
}

save(ratio_mat, file = paste0("Result/", data_freq, "double descent From", min(group_size_vec), "to", max(group_size_vec), ".Rdata"))
stopCluster(cl)


# Figure for double descent top8
library(ggplot2)
library(ggpubr)
library(patchwork)

load("Result/Monthlydouble descent From4to1024.Rdata")
ratio_all = data.frame(ratio_mat)

load("Result/Dailydouble descent From4to1024.Rdata")
ratio_all = data.frame(rbind(ratio_all,
                             ratio_mat))

ratio_select = ratio_all[,-9]
group_size_vec = 4^(1:5)
custom_names <- c("EW", "S:Linear", "M=1", paste0("M=", group_size_vec))
colnames(ratio_select)[1:8] = custom_names
theme_slides = theme(text=element_text(size=15),
                     legend.position = "top")

df_long <- ratio_select %>%
  as.data.frame() %>% # Ensure it's a data frame if it's currently a matrix
  pivot_longer(
    cols = custom_names,       # Reshape all 10 columns
    names_to = "Variable",     # This column will hold your original column names
    values_to = "Value"        # This column will hold the 21,790 numeric observations
  )

df_long$Variable = factor(df_long$Variable, 
                          levels = custom_names)

df_long = df_long %>% filter(dataset %in% c(paste0("Monthly",1:7), "Daily1"))

df_long$dataset = factor(df_long$dataset,
                         levels = c("Monthly1", "Monthly2", "Monthly3", "Monthly4", "Daily1", "Monthly5", "Monthly6", "Monthly7"),
                         labels = c("(a) Monthly1, #Pooled = 6782", "(b) Monthly2, #Pooled = 6258",
                                    "(c) Monthly3, #Pooled = 2643", "(d) Monthly4, #Pooled = 1517",
                                    "(e) Daily1, #Pooled = 1506", "(f) Monthly5, #Pooled = 1285",
                                    "(g) Monthly6, #Pooled = 1225", "(h) Monthly7, #Pooled = 1048"))

DS_top8 = ggplot(df_long, aes(x = Variable, y = Value)) +
  geom_boxplot(outlier.alpha = 0.05, notch = FALSE) +
  scale_y_log10(breaks = 10^(-3:3), labels = c("0.001", "0.010", "0.100", "1.000", "10.000", "100.000", "1000.000")) +
  annotation_logticks(sides = "l", short = unit(0.5, "mm"), mid = unit(0.5, "mm"), long = unit(1.5, "mm")) +
  coord_cartesian(ylim = c(0.001, 1000)) +
  geom_hline(yintercept = 1, color = "red") +
  labs(
    x = "Model Complexity",
    y = "Loss Ratio (Benchmarks/P:Ridgeless)"
  ) +
  theme_bw() + 
  theme(legend.position = "none",
        axis.text.x = element_text(size = 6, angle = 0, hjust = 0.5),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10)) + 
  theme_slides

p1 <- DS_top8 + facet_wrap(~ dataset, nrow = 4)
print(p1)

ggsave(
  filename = "Result/DoubleDescent_top8.pdf",
  plot = p1,
  width = 8.27,
  height = 11.69,
  units = "in",
  dpi = 300,
  device = "pdf"
)