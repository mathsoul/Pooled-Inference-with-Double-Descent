## ============================================================================
## Tests whether the log-variance of pairwise expert forecast-error contrasts
## decomposes additively into a variable effect and an expert effect, via a 
## two-way-fixed-effects (TWFE) regression:
##
##     log( var(x_{m,n}) ) = grand_mean + variable_effect[m] + expert_effect[n] + resid
##
## where m indexes variable (an M4 time series, or an HHS region for flu)
## and n indexes one of the expert
## ============================================================================

library(tidyverse)
library(Matrix)
library(data.table)

source("Code/func.R")

theme_slides = theme(text = element_text(size = 15), legend.position = "top")

## ----------------------------------------------------------------------------
## Flu dataset
## ----------------------------------------------------------------------------

csv_path = "CleanedData/FluForecasting/point_ests_adj-w20172018.csv"

df = fread(csv_path) %>%
  filter(target == "1 wk ahead", location != "US National",
         !model_name %in% c("ReichLab_kde", "UTAustin_edm",
                             "constant-weights", "equal-weights",
                             "target-and-region-based-weights",
                             "target-based-weights",
                             "target-type-based-weights"))

# df = fread(csv_path) %>%
#   filter(target == "2 wk ahead", location != "US National",
#          !model_name %in% c("ReichLab_kde", "UTAustin_edm",
#                             "constant-weights", "equal-weights",
#                             "target-and-region-based-weights",
#                             "target-based-weights",
#                             "target-type-based-weights"))

locations = sort(unique(df$location))
n_prods = length(locations)          # 10 HHS regions ("products")
seasons = sort(unique(df$Season))    # 8 flu seasons

results = list()
for(season in seasons){
  wide = prepareFluData(df, season)
  n_experts = length(unique(wide$model_name))

  err_mat = wide %>% dplyr::select(-location, -model_name) %>% as.matrix()
  var_mat = computeVarMat(err_mat, n_experts, n_prods)
  r2_vec = getTWFE_R2(var_mat)

  results[[season]] = data.frame(
    Season = season, n_experts = n_experts, n_prods = n_prods,
    n_periods = ncol(err_mat), R2 = r2_vec["R2"], adjR2 = r2_vec["adjR2"]
  )

  cat(sprintf("Flu 1 wk ahead season %s: n_experts=%d n_prods=%d weeks=%d R2=%.4f adjR2=%.4f\n",
              season, n_experts, n_prods, ncol(err_mat), r2_vec["R2"], r2_vec["adjR2"]))
  
  # cat(sprintf("Flu 2 wk ahead season %s: n_experts=%d n_prods=%d weeks=%d R2=%.4f adjR2=%.4f\n",
  #             season, n_experts, n_prods, ncol(err_mat), r2_vec["R2"], r2_vec["adjR2"]))
}

r2_df_flu = do.call(rbind, results)
rownames(r2_df_flu) = NULL
r2_df_flu$Season = factor(r2_df_flu$Season, levels = seasons)

p_bar = ggplot(r2_df_flu, aes(x = Season, y = adjR2)) +
  geom_col(fill = "grey70", color = "black") +
  geom_text(aes(label = round(adjR2, 3)), vjust = -0.5, size = 4) +
  theme_bw() + theme_slides +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_cartesian(ylim = c(0, 1)) +
  ylab("Adjusted R2") + xlab("Season")

print(p_bar)
ggsave("Result/FluGoodnessOfFit_bySeason_bar.pdf", plot = p_bar, width = 9, height = 6)

## ----------------------------------------------------------------------------
## M4 dataset
## ----------------------------------------------------------------------------

n_experts = 17
scaled_data_dir = "CleanedData/M4ScaledData"

results = list()

for(data_freq in c("Monthly", "Daily")){
  u_data = fread(file.path(scaled_data_dir, paste0(data_freq, "_err.csv"))) %>% dplyr::select(-V1)
  n_scenarios = if(data_freq == "Monthly") 68 else 2

  for(rank_idx in 1:n_scenarios){
    prep = prepareM4Scenario(data_freq, rank_idx, u_data)
    var_mat = computeVarMat(prep$err_mat, n_experts, prep$n_prods)
    r2_vec = getTWFE_R2(var_mat)

    scenario = paste0(data_freq, rank_idx)
    results[[scenario]] = data.frame(
      scenario = scenario, data_freq = data_freq, rank_idx = rank_idx,
      n_prods = prep$n_prods, n_periods = ncol(prep$err_mat),
      R2 = r2_vec["R2"], adjR2 = r2_vec["adjR2"]
    )

    cat(sprintf("%s: n_prods=%d R2=%.4f adjR2=%.4f\n",
                scenario, prep$n_prods, r2_vec["R2"], r2_vec["adjR2"]))
  }
}

r2_df_m4 = do.call(rbind, results)
rownames(r2_df_m4) = NULL

save(r2_df_m4, file = "Result/M4GoodnessOfFit_70scenarios.RData")
write.csv(r2_df_m4, "Result/M4GoodnessOfFit_70scenarios.csv", row.names = FALSE)

# Histogram of the ~70 adjusted R^2 values (linear x-axis). Saved as PDF.
p_hist = ggplot(r2_df_m4, aes(x = adjR2)) +
  geom_histogram(bins = 15, color = "black", fill = "white") +
  theme_bw() + theme_slides +
  xlab("Adjusted R2") + ylab("Count")

print(p_hist)
ggsave("Result/M4GoodnessOfFit_histogram.pdf", plot = p_hist, width = 8, height = 6)

# Scatterplot of adjusted R^2 values (y) vs M (x).
p_scatter = ggplot(r2_df_m4, aes(x = n_prods, y = adjR2)) +
  geom_point() + coord_cartesian(xlim = c(100,10000)) + 
  theme_bw() + theme_slides + scale_x_log10() + 
  ylab("Adjusted R2") + xlab("Number of Varibles") + annotation_logticks(sides = "b")

print(p_scatter)
