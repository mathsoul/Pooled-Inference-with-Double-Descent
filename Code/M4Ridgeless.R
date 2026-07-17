library(tidyverse)
library(Matrix)
library(MCMCpack)
library(doParallel)
library(doRNG)
library(data.table)
pkg_names = c("tidyverse", "nlshrink", "MASS")

source("Code/func.R")


n_cores = detectCores() - 1
cl <- makeCluster(n_cores)
registerDoParallel(cl)

# Retrieve all function names from the global environment in preparation for parallel computation 
all_objects <- ls(envir = .GlobalEnv)
# Filter to include only functions
all_functions <- all_objects[sapply(mget(all_objects, envir = .GlobalEnv), is.function)]

data_freq_vec = c("Hourly", "Daily", "Monthly")

limits_df = data.frame(lower = rep(1,3), upper = c(1, 2, 68))

methods = c("EW", "Linear", "Cor", "Var")
n_methods = length(methods)

ratio_mat = NULL

for(k in 1:length(data_freq_vec)){
  data_freq = data_freq_vec[k]
  rank_vec = limits_df$lower[k]:limits_df$upper[k]
  
  f_data = fread(paste0("CleanedData/M4ScaledData/", data_freq, "_forecast.csv")) %>% dplyr::select(-V1)
  true_data = fread(paste0("CleanedData/M4ScaledData/", data_freq, "_truth.csv")) %>% dplyr::select(-V1)
  u_data = fread(paste0("CleanedData/M4ScaledData/", data_freq, "_err.csv")) %>% dplyr::select(-V1)
  
  if(data_freq == "Monthly"){
    pb <- txtProgressBar(min = min(rank_vec), max = max(rank_vec), style = 3, 
                         char = "=", width = 50)
  }
  
  # It takes about 40 mins using a MacBook Pro with M4 chip and 16GB RAM
  for(i in 1:length(rank_vec)){
    if(data_freq == "Monthly"){
      setTxtProgressBar(pb, i)
    }
    
    rank_idx = rank_vec[i]
    
    time_series_names = getTimeSeriesNames(data_freq, rank_idx)
    
    n_prods = length(time_series_names)
    n_experts = 17
    n_periods = ncol(true_data) - 1
    
    true_comb = as.data.frame(true_data %>% filter(id %in% time_series_names[1:n_prods]) %>% 
                                arrange(as.numeric(gsub("\\D", "", id))) %>% dplyr::select(-id))
    f_comb = as.data.frame(f_data %>% filter(id %in% time_series_names[1:n_prods]) %>% 
                             arrange(as.numeric(gsub("\\D", "", id))) %>% dplyr::select(-id, -group_id))
    u_comb = as.data.frame(u_data %>% filter(id %in% time_series_names[1:n_prods]) %>% 
                             arrange(as.numeric(gsub("\\D", "", id))) %>% dplyr::select(-id, -group_id))
    
    
    # P:Linear ----------------------------------------------------------------
    pool_idx = rep(1:n_experts, n_prods) + rep(0:(n_prods - 1) * n_experts, each = n_experts)
    
    # The following is too slow when n_prods is large
    # u_pool_Linear1 = getPooledU(pool_idx, "Linear", n_prods, n_experts, f_comb, u_comb, true_data)
    u_pool_Linear = suppressMessages(getPooledUFast(pool_idx, n_prods, n_experts, u_comb, true_comb))
    
    L2_PLinear = apply(u_pool_Linear, 1, function(x) mean(x^2))
    
    
    # P:Ridge0+ ---------------------------------------------------------------
    Y_mat = t(apply(u_comb, 2, function(x) sumEachProd(x, n_prods, n_experts)/n_experts))
    
    Gamma0 = getGamma0(n_experts)
    Gamma = Matrix(diag(n_prods), sparse = TRUE) %x% Gamma0
    
    X_mat = as.matrix(t(u_comb) %*% Gamma)
    
    ridge0_test_err_std = suppressMessages(map_dfc(1:nrow(X_mat), ~ getRidge0TestErr(X_mat[-.x,], Y_mat[-.x,], X_mat[.x,,drop = FALSE],
                                                                                     Y_mat[.x,,drop = FALSE], stdize = TRUE)))
    
    PCR80_test_err = suppressMessages(map_dfc(1:nrow(X_mat), ~ getPCRTestErr(X_mat[-.x,], Y_mat[-.x,], X_mat[.x,,drop = FALSE],
                                                                             Y_mat[.x,,drop = FALSE], threshold = 0.8)))
    
    PCR70_test_err = suppressMessages(map_dfc(1:nrow(X_mat), ~ getPCRTestErr(X_mat[-.x,], Y_mat[-.x,], X_mat[.x,,drop = FALSE],
                                                                             Y_mat[.x,,drop = FALSE], threshold = 0.7)))
    
    L2_ridge0_std = apply(ridge0_test_err_std, 1, function(x) mean(x^2))
    L2_PCR_80 = apply(PCR80_test_err, 1, function(x) mean(x^2))
    L2_PCR_70 = apply(PCR70_test_err, 1, function(x) mean(x^2))
    
    Gamma0_SeqDiff = getGamma0(n_experts, type = "sequential difference")
    Gamma_SeqDiff = Matrix(diag(n_prods), sparse = TRUE) %x% Gamma0_SeqDiff
    
    X_mat_SeqDiff = as.matrix(t(u_comb) %*% Gamma_SeqDiff)
    
    ridge0_test_err_SeqDiff = suppressMessages(map_dfc(1:nrow(X_mat_SeqDiff), ~ getRidge0TestErr(X_mat_SeqDiff[-.x,], Y_mat[-.x,], X_mat_SeqDiff[.x,,drop = FALSE],
                                                                                                 Y_mat[.x,,drop = FALSE], stdize = TRUE)))
    
    L2_ridge0_SeqDiff = apply(ridge0_test_err_SeqDiff, 1, function(x) mean(x^2))
    
    # Separate Methods --------------------------------------------------------
    
    L2_list = foreach(i = 1:n_prods,.options.RNG = 20240422,
                      .packages = pkg_names, .export = all_functions)%dorng%{
                        single_prod_idx = 1:n_experts + (i-1) * n_experts
                        
                        true_sep = true_comb[i,]
                        f_sep = tibble(f_comb[single_prod_idx,])
                        u_sep = u_comb[single_prod_idx,]
                        
                        u_methods = suppressMessages(getSepUMethods(1:n_experts, methods, n = n_experts, f_sep, u_sep, true_sep))
                        
                        apply(u_methods, 1, function(x) mean(x^2))
                      }
    
    L2_mat = matrix(unlist(L2_list), nrow = n_prods, byrow = TRUE)
    # L2_mat = matrix(NA, nrow = n_prods, ncol = n_methods)
    colnames(L2_mat) = methods
    
    ratio_mat_1start = cbind(L2_PLinear, L2_PCR_80, L2_PCR_70, L2_ridge0_SeqDiff, L2_mat)/L2_ridge0_std
    
    # print(summary(ratio_mat_1start))
    
    ratio_mat = rbind(ratio_mat, 
                      data.frame(ratio_mat_1start, dataset = paste0(data_freq, rank_idx)))
  }
}

save(ratio_mat, file = paste0("Result/M4.Rdata"))

summary(ratio_mat)

stopCluster(cl)

