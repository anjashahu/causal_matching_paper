
### SET UP ---------------------------------------------------------------------

# load packages
library(tidyverse)

# clear working environment
rm(list = ls())

# specify which of the individual data frames you want to use
specify_df <- "_01"
specify_df

# load true population effects
true_pop_effects <- read.csv(paste0("./data/simulated_pop_data/true_pop_or", specify_df, ".csv"))
rownames(true_pop_effects) <- c("mci", "hyp")


# load simulation results

# main results
#load(paste0("./data/simulation_results/main_results/psm_analyses", specify_df, ".RData"))
#load(paste0("./data/simulation_results/main_results/cem_analyses", specify_df, ".RData"))
#load(paste0("./data/simulation_results/main_results/weighting_analyses", specify_df, ".RData"))

# no bmi results
#load(paste0("./data/simulation_results/confounding_sensitivity_results/remove_bmi_results/psm_analyses_no_bmi", specify_df,".RData"))
#load(paste0("./data/simulation_results/confounding_sensitivity_results/remove_bmi_results/cem_analyses_no_bmi", specify_df,".RData"))
#load(paste0("./data/simulation_results/confounding_sensitivity_results/remove_bmi_results/weighting_analyses_no_bmi", specify_df,".RData"))

# no age results
#load(paste0("./data/simulation_results/confounding_sensitivity_results/remove_age_results/psm_analyses_no_age", specify_df,".RData"))
#load(paste0("./data/simulation_results/confounding_sensitivity_results/remove_age_results/cem_analyses_no_age", specify_df,".RData"))
#load(paste0("./data/simulation_results/confounding_sensitivity_results/remove_age_results/weighting_analyses_no_age", specify_df,".RData"))

### CALC AVG EFFECT --------------------------------------------------------------

# list of outcomes
outcomes <- names(list_psm_effect_or)
# create list to hold results of both outcomes
list_avg_effect <- vector(mode = "list", length = length(outcomes))
names(list_avg_effect) <- outcomes

# loop through both outcomes to calculate average effect estimate for each outcome
for (outcome in outcomes) {
  
  # calculate average effect estimate of the outcome for each method and collect into one vector
  avg_effect_outcome <- c(
    # PSM
    apply(list_psm_effect_or[[outcome]], 2, mean),
    # weighting
    apply(list_weighting_effect_or[[outcome]], 2, mean),
    # CEM
    apply(list_cem_effect_or[[outcome]], 2, mean)
  )
  # make labels match for all the outcomes
  names(avg_effect_outcome) <- str_replace(names(avg_effect_outcome), paste0("_", outcome), "")

  # store into associated entry for the outcome in the list
  list_avg_effect[[outcome]] <- avg_effect_outcome
  
}

# check to see that the two outcomes calculated the same estimates in the same order
identical(names(list_avg_effect$mci), names(list_avg_effect$hyp))
# since the same estimates were calculated, collapse the list into a single data frame
avg_effect <- data.frame(
  estimate =  names(list_avg_effect$mci),
  matrix(unlist(list_avg_effect), ncol = length(list_avg_effect), byrow = FALSE)
)
# name the columns 
colnames(avg_effect) <- c("estimate", outcomes)
avg_effect

### CALC BIAS ------------------------------------------------------------------
# associate each estimate with a true pop effect 
marg_ate_estimates <- c(
  "psw_weighted_iptw_unadj_svy", 
  "psw_weighted_iptw_osw_unadj_svy", 
  "psw_cov_iptw_unadj_svy", 
  "psw_cov_iptw_osw_unadj_svy"
)
cond_ate_estimates <- c(
  "psw_weighted_iptw_adj_svy", 
  "psw_weighted_iptw_osw_adj_svy", 
  "psw_cov_iptw_adj_svy", 
  "psw_cov_iptw_osw_adj_svy"
)
marg_att_estimates <- avg_effect$estimate[!(avg_effect$estimate %in% c(marg_ate_estimates, cond_ate_estimates))] %>%
  str_extract(".*_unadj.*") %>%
  discard(is.na)
cond_att_estimates <- avg_effect$estimate[!(avg_effect$estimate %in% c(marg_ate_estimates, cond_ate_estimates))] %>%
  str_extract(".*_adj.*") %>%
  discard(is.na)

# add columns to data frame for true pop effect and bias
avg_effect <- avg_effect %>%
  mutate(
    effect = case_when(
      estimate %in% marg_ate_estimates ~ "marg_ate",
      estimate %in% cond_ate_estimates ~ "cond_ate",
      estimate %in% marg_att_estimates ~ "marg_att",
      estimate %in% cond_att_estimates ~ "cond_att",
    )
  ) %>% 
  rowwise() %>%
  mutate(
    true_effect_mci = true_pop_effects["mci", effect],
    true_effect_hyp = true_pop_effects["hyp", effect],
    bias_mci = mci - true_effect_mci,
    bias_hyp = hyp - true_effect_hyp
  )

# look at mci bias
avg_effect %>% 
  select(estimate, mci, bias_mci) 
# look at hypertension bias 
avg_effect %>%
  select(estimate, hyp, bias_hyp) 


### CALC COVERAGE --------------------------------------------------------------

# find critical z for calculating 95% CIs
critical_z <- qnorm(0.975) 
# combine all effect estimates and SEs into one list for both outcome
# effect estimates
list_comb_effect <- list(
  mci = data.frame(list_psm_effect_or[["mci"]], list_weighting_effect_or[["mci"]], list_cem_effect_or[["mci"]]),
  hyp = data.frame(list_psm_effect_or[["hyp"]], list_weighting_effect_or[["hyp"]], list_cem_effect_or[["hyp"]])
)
# SEs
list_comb_se <- list(
  mci = data.frame(list_psm_se[["mci"]], list_weighting_se[["mci"]], list_cem_se[["mci"]]),
  hyp = data.frame(list_psm_se[["hyp"]], list_weighting_se[["hyp"]], list_cem_se[["hyp"]])
)
# fix column names of each outcome data frame to match estimate names in avg_effect data frame
for (outcome in outcomes) {
  colnames(list_comb_effect[[outcome]]) <- str_replace(colnames(list_comb_effect[[outcome]]), paste0("_", outcome), "")
  colnames(list_comb_se[[outcome]]) <- str_replace(colnames(list_comb_se[[outcome]]), paste0("_", outcome), "")
}

# create function to check whether true value is captured in the interval
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

# create data frame to hold coverage results
coverage <- data.frame(estimate = avg_effect$estimate)
# loop through the outcomes and the estimates to calculate coverage for each estimate for each outcome
for (outcome in outcomes) {
  # create empty vector for the current outcome to store whether the different estimates' 95% CI captured the true value
  est_coverage <- rep(NA, length(coverage$estimate))
  # loop through the estimates to see what proportion of the 95% CIs for the repetitions of that estimate captured the true value
  for (est in coverage$estimate) {
    est_lower_log_or <- log(list_comb_effect[[outcome]][, est]) - critical_z * list_comb_se[[outcome]][, est]
    est_upper_log_or <- log(list_comb_effect[[outcome]][, est]) + critical_z * list_comb_se[[outcome]][, est]
    est_lower <- exp(est_lower_log_or)
    est_upper <- exp(est_upper_log_or)
    true_effect <- avg_effect %>% 
      filter(estimate == est) %>% 
      .[, paste0("true_effect_", outcome)]
    est_coverage[which(coverage$estimate == est)] <- check_if_in_interval(est_lower, est_upper, true_effect)
  }
  coverage[, paste0("cov_", outcome)] <- est_coverage
}

coverage

### COMBINE RESULTS -----------------------------------------------------------------------------

# join bias and coverage results together into one data frame
bias_coverage <- avg_effect %>% 
  left_join(coverage, by = "estimate") %>%
  rename(est_mci = mci, est_hyp = hyp)
# save
#save(bias_coverage, file = paste0("./data/simulation_results/bias_coverage_results/bias_cov_main", specify_df, ".RData"))
#save(bias_coverage, file = paste0("./data/simulation_results/bias_coverage_results/bias_cov_no_bmi", specify_df, ".RData"))
#save(bias_coverage, file = paste0("./data/simulation_results/bias_coverage_results/bias_cov_no_age", specify_df, ".RData"))


