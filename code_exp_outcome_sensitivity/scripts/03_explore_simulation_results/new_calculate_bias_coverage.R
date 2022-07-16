
### SET UP ---------------------------------------------------------------------

# load packages
library(tidyverse)

# clear working environment
rm(list = ls())

# load simulation results
# set number of data frames
n <- 16
# set what results you want 
analysis_type <- "main"
#analysis_type <- "confounding"

if (analysis_type == "main") {
  for (df_num in 1:n) {
    assign(paste0("true_pop_effects_", df_num), read.csv(paste0("./data/simulated_pop_data/true_pop_or/true_pop_or_", df_num, ".csv"), row.names = c("mci", "hyp")))
    assign(paste0("list_effect_or_", df_num), readRDS(paste0("./data/simulation_results/main_results/analyses_effect_", df_num, ".RData")))
    assign(paste0("list_se_", df_num), readRDS(paste0("./data/simulation_results/main_results/analyses_se_", df_num, ".RData")))
  }  
} else if (analysis_type == "confounding") {
  for (df_num in 1:n) {
    assign(paste0("true_pop_effects_", df_num), read.csv(paste0("./data/simulated_pop_data/true_pop_or/true_pop_or_", df_num, ".csv"), row.names = c("mci", "hyp")))
    assign(paste0("list_effect_or_", df_num), readRDS(paste0("./data/simulation_results/confounding_sensitivity_no_age/analyses_effect_", df_num, ".RData")))
    assign(paste0("list_se_", df_num), readRDS(paste0("./data/simulation_results/confounding_sensitivity_no_age/analyses_se_", df_num, ".RData")))
  }    
}


### CALC AVG EFFECT --------------------------------------------------------------


check_if_in_interval <- function(lower, upper, true) {
  # create empty vector to store whether the repetition's 95% CI captured the true value
  in_interval <- rep(NA, length(lower))
  # assess whether each of the repetitions' 95% CIs capture the true value
  for (i in 1:length(lower)) {
    in_interval[i] <- true >= lower[i] & true <= upper[i]
  } 
  # take the mean over all the repetitions to get the coverage
  coverage <- mean(in_interval)
  # return coverage result
  return(coverage)
}

obtain_results <- function(outcome, list_effect_or, list_se, true_pop_effects){
  df_results <- data.frame(method = c("psm, unadjusted", "psm, adjusted", "cem, unadjusted", "cem, adjusted"))
  df_results$effect <- apply(list_effect_or[[outcome]], 2, mean)
  
  df_results$true <- rep(unlist(true_pop_effects[outcome,]), 2)
  df_results$bias <- df_results$effect - df_results$true
  
  check_if_in_interval <- function(lower, upper, true) {
    # create empty vector to store whether the repetition's 95% CI captured the true value
    in_interval <- rep(NA, length(lower))
    # assess whether each of the repetitions' 95% CIs capture the true value
    for (i in 1:length(lower)) {
      in_interval[i] <- true >= lower[i] & true <= upper[i]
    } 
    # take the mean over all the repetitions to get the coverage
    coverage <- mean(in_interval)
    # return coverage result
    return(coverage)
  }
  
  critical_z <- qnorm(0.975) 
  coverage_all <- rep(NA, 4)
  for (i in 1:4) {
    true <- rep(unlist(true_pop_effects[outcome,]), 2)[i]
    est_lower_log_or <- log(list_effect_or[[outcome]][,i]) - critical_z * list_se[[outcome]][,i] 
    est_upper_log_or <- log(list_effect_or[[outcome]][,i]) + critical_z * list_se[[outcome]][,i] 
    est_lower <- exp(est_lower_log_or)
    est_upper <- exp(est_upper_log_or)
    coverage_all[i] <- check_if_in_interval(est_lower, est_upper, true)
  }
  df_results$coverage <- coverage_all
  return(df_results)

}

df_scenarios <- expand.grid(c(5, 15, 25, 35), c(5, 15, 25, 35)) %>% 
  data.frame() %>% 
  rename(outcome_prev = Var1, exposure_prev = Var2) %>%
  mutate(scenario = as.numeric(row.names(.)))

list_results_mci <- vector("list", length = n)
list_results_hyp <- vector("list", length = n)

for (i in 1:n){
  list_effect_or <- get(paste0("list_effect_or_", i))
  list_se <- get(paste0("list_se_", i))
  true_pop_effects <- get(paste0("true_pop_effects_", i))
  list_results_mci[[i]] <- obtain_results("mci", list_effect_or, list_se, true_pop_effects) %>%
    mutate(scenario = rep(i, 4)) %>%
    left_join(df_scenarios, by = "scenario") 
  list_results_hyp[[i]] <- obtain_results("hyp", list_effect_or, list_se, true_pop_effects) %>%
    mutate(scenario = rep(i, 4)) %>%
    left_join(df_scenarios, by = "scenario")
}

df_mci <- do.call(rbind, list_results_mci) 
df_hyp <- do.call(rbind, list_results_hyp) 

### SAVE RESULTS ------------------------------------

saveRDS(df_mci, file = paste0("./data/simulation_results/bias_coverage_results/effects_mci_", analysis_type, ".RData"))
saveRDS(df_hyp, file = paste0("./data/simulation_results/bias_coverage_results/effects_hyp_", analysis_type, ".RData"))




