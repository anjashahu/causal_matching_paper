
### SET UP ---------------------------------------------------------------------

# load packages
library(tidyverse)
library(survey)

# clear working environment
rm(list = ls())

# specify which of the individual data frames you want to use
specify_df <- "_01"
specify_df

# load data 
df_pop_bg <- read_csv("./data/simulated_pop_data/df_pop_bg.csv")
df_pop_hh <- read_csv("./data/simulated_pop_data/df_pop_hh.csv")
df_pop_indiv <- read_csv(paste0("./data/simulated_pop_data/df_pop_indiv", specify_df, ".csv")) # change this path to load the data frame of interest

# set numbers of repetitions
reps <- 1000
reps

# note: for this sensitivity analysis, we are pretending that we do not have BMI information
# to see the impact of confounding on our results

### SIMULATE SAMPLES --------------------------------------------------------
# run R script to obtain function to simulate samples
source("./scripts/02_conduct_simulation_analyses/01_simulate_samples_function.R")

# create empty list to store samples
list_samps <- vector(mode = "list", length = reps)
# create empty list to store survey objects of samples that use OSW
list_samps_svy <- vector(mode = "list", length = reps)

# simulate all the samples and their survey design object by looping through the number of repetitions and storing each sample into list
for (rep in 1:reps) {
  if((rep %% 100) == 0) print(rep)
  # simulate sample for the repetition and add to associated list
  list_samps[[rep]] <- sim_one_sample(df_pop_bg, df_pop_hh, df_pop_indiv, rep) 
  # simulate survey design object for repetition and add to associated list
  list_samps_svy[[rep]] <- svydesign(
    ids = ~ bg + hh, 
    weights = ~ final_weights, 
    strata = ~ strata, 
    data = list_samps[[rep]]
  )
}

# check to make sure that none of the strata end up with one BG only (as this will not allow survey regressions to be fit)
#problems <- rep(NA, reps)
#for (rep in 1:reps) {
#print(list_samps[[rep]] %>% count(strata, bg) %>% count(strata) %>% pull(n))
#check_rep <- list_samps[[rep]] %>% count(strata, bg) %>% count(strata) %>% pull(n)
#problems[rep] <- sum(check_rep == 0 | check_rep == 1)
#}
#sum(problems > 0) # so it looks like there won't be any issues

# remove unneeded objects
rm(rep, list = ls(pattern = "df_pop|cuts"), sim_one_sample)

### CONDUCT MATCHING -----------------------------------------------------------

# run R script to obtain function to conduct matching
source("./scripts/02_conduct_simulation_analyses/04_confounding_sensitivity_analyses/02_remove_bmi_analyses/02.01_conduct_psm_function_no_bmi.R")

# create empty lists to store matches using propensity scores 
# where propensity score is based on weighted logistic regression using OSW (original survey weights)
list_matches_psm_weighted <- vector(mode = "list", length = reps)
# where propensity score is based on logistic regression with OSW as covariate
list_matches_psm_cov <- vector(mode = "list", length = reps)
# loop through samples to obtain matches 
for (rep in 1:reps){
  if((rep %% 100) == 0) print(rep)
  psm_rep_results <- conduct_one_psm(list_samps[[rep]])
  list_matches_psm_weighted[[rep]] <- psm_rep_results[[1]]
  list_matches_psm_cov[[rep]] <- psm_rep_results[[2]]
}

# remove unneeded objects
rm(rep, psm_rep_results, get_inherited_weights, list_samps, conduct_one_psm)

### CREATE SURVEY DESIGN OBJECTS -------------------------------------------------

# for PSM 
# create empty list to store survey objects of PSM matches that use OSW
# where propensity score is calculated using weighted logistic regression 
list_matches_psm_weighted_osw_svy <- vector(mode = "list", length = reps)
# where propensity score is calculated using logistic regression with weight as covariate
list_matches_psm_cov_osw_svy <- vector(mode = "list", length = reps)
# create empty list to store survey objects of PSM matches that use ISW
# where propensity score is calculated using weighted logistic regression 
list_matches_psm_weighted_isw_svy <- vector(mode = "list", length = reps)
# where propensity score is calculated using logistic regression with weight as covariate
list_matches_psm_cov_isw_svy <- vector(mode = "list", length = reps)

# loop through the PSM matches, create survey design objects with different weights and add to associated list
for (rep in 1:reps) {
  if((rep %% 100) == 0) print(rep)
  
  # for PSM matches using OSW
  # where propensity score is calculated using weighted logistic regression 
  list_matches_psm_weighted_osw_svy[[rep]] <- svydesign(
    ids = ~ bg + hh, 
    weights = ~ final_weights, 
    strata = ~ strata, 
    data = list_matches_psm_weighted[[rep]]
  )
  # where propensity score is calculated using logistic regression with weight as covariate
  list_matches_psm_cov_osw_svy[[rep]] <- svydesign(
    ids = ~ bg + hh, 
    weights = ~ final_weights, 
    strata = ~ strata, 
    data = list_matches_psm_cov[[rep]]
  )
  
  # for PSM matches that use ISW
  # where propensity score is calculated using weighted logistic regression 
  list_matches_psm_weighted_isw_svy[[rep]] <- svydesign(
    ids = ~ bg + hh, 
    weights = ~ inherited_weights, 
    strata = ~ strata, 
    data = list_matches_psm_weighted[[rep]]
  )
  # where propensity score is calculated using logistic regression with weight as covariate
  list_matches_psm_cov_isw_svy[[rep]] <- svydesign(
    ids = ~ bg + hh, 
    weights = ~ inherited_weights, 
    strata = ~ strata, 
    data = list_matches_psm_cov[[rep]]
  )
  
}

# remove unneeded variables
rm(rep)

### CONDUCT ANALYSES MCI -------------------------------------------------------------

# list the type of analyses we will conduct for MCI for psm
psm_analyses_mci <- c(
  "psm_weighted_unadj_mci", "psm_weighted_unadj_svy_osw_mci", "psm_weighted_unadj_svy_isw_mci",
  "psm_weighted_adj_mci", "psm_weighted_adj_svy_osw_mci", "psm_weighted_adj_svy_isw_mci",
  "psm_cov_unadj_mci", "psm_cov_unadj_svy_osw_mci", "psm_cov_unadj_svy_isw_mci",
  "psm_cov_adj_mci", "psm_cov_adj_svy_osw_mci", "psm_cov_adj_svy_isw_mci"
)
# create matrix to hold effect estimates from the various analyses
psm_effect_mci <- matrix(NA, nrow = reps, ncol = length(psm_analyses_mci), dimnames = list(NULL, psm_analyses_mci))
# create matrix to hold SEs from the various analyses
psm_se_mci <- matrix(NA, nrow = reps, ncol = length(psm_analyses_mci), dimnames = list(NULL, psm_analyses_mci))

# loop through samples to conduct the various analyses on the matched data
for (rep in 1:reps) {
  if((rep %% 100) == 0) print(rep)
  
  # fit the various regressions
  
  # when propensity score is calculated using weighted logistic regression 
  # unadjusted regression without weights
  fit_psm_weighted_unadj_mci <- glm(
    mci ~ insomnia, 
    family = binomial(), 
    data = list_matches_psm_weighted[[rep]]
  )
  # survey unadjusted regression using OSW
  fit_psm_weighted_unadj_svy_osw_mci <- svyglm(
    mci ~ insomnia, 
    family = quasibinomial(), 
    design = list_matches_psm_weighted_osw_svy[[rep]]
  )
  # survey unadjusted regression using ISW
  fit_psm_weighted_unadj_svy_isw_mci <- svyglm(
    mci ~ insomnia, 
    family = quasibinomial(), 
    design = list_matches_psm_weighted_isw_svy[[rep]]
  )
  # adjusted regression without weights
  fit_psm_weighted_adj_mci <- glm(
    mci ~ insomnia + age, 
    family = binomial(), 
    data = list_matches_psm_weighted[[rep]]
  )
  # survey adjusted regression using OSW
  fit_psm_weighted_adj_svy_osw_mci <- svyglm(
    mci ~ insomnia + age, 
    family = quasibinomial(), 
    design = list_matches_psm_weighted_osw_svy[[rep]]
  )
  # survey adjusted regression using ISW
  fit_psm_weighted_adj_svy_isw_mci <- svyglm(
    mci ~ insomnia + age, 
    family = quasibinomial(), 
    design = list_matches_psm_weighted_isw_svy[[rep]]
  )
  
  # when propensity score is calculated using logistic regression with weight as covariate
  # unadjusted regression without weights
  fit_psm_cov_unadj_mci <- glm(
    mci ~ insomnia, 
    family = binomial(), 
    data = list_matches_psm_cov[[rep]]
  )
  # survey unadjusted regression using OSW
  fit_psm_cov_unadj_svy_osw_mci <- svyglm(
    mci ~ insomnia, 
    family = quasibinomial(), 
    design = list_matches_psm_cov_osw_svy[[rep]]
  )
  # survey unadjusted regression using ISW
  fit_psm_cov_unadj_svy_isw_mci <- svyglm(
    mci ~ insomnia, 
    family = quasibinomial(), 
    design = list_matches_psm_cov_isw_svy[[rep]]
  )
  # adjusted regression without weights
  fit_psm_cov_adj_mci <- glm(
    mci ~ insomnia + age, 
    family = binomial(), 
    data = list_matches_psm_cov[[rep]]
  )
  # survey adjusted regression using OSW
  fit_psm_cov_adj_svy_osw_mci <- svyglm(
    mci ~ insomnia + age, 
    family = quasibinomial(), 
    design = list_matches_psm_cov_osw_svy[[rep]]
  )
  # survey adjusted regression using ISW
  fit_psm_cov_adj_svy_isw_mci <- svyglm(
    mci ~ insomnia + age, 
    family = quasibinomial(), 
    design = list_matches_psm_cov_isw_svy[[rep]]
  )
  
  # add result to associated matrix
  # effect estimate results
  psm_effect_mci[rep, ] <- c(
    summary(fit_psm_weighted_unadj_mci)$coefficients["insomnia", "Estimate"],
    summary(fit_psm_weighted_unadj_svy_osw_mci)$coefficients["insomnia", "Estimate"],
    summary(fit_psm_weighted_unadj_svy_isw_mci)$coefficients["insomnia", "Estimate"],
    summary(fit_psm_weighted_adj_mci)$coefficients["insomnia", "Estimate"],
    summary(fit_psm_weighted_adj_svy_osw_mci)$coefficients["insomnia", "Estimate"],
    summary(fit_psm_weighted_adj_svy_isw_mci)$coefficients["insomnia", "Estimate"],
    
    summary(fit_psm_cov_unadj_mci)$coefficients["insomnia", "Estimate"],
    summary(fit_psm_cov_unadj_svy_osw_mci)$coefficients["insomnia", "Estimate"],
    summary(fit_psm_cov_unadj_svy_isw_mci)$coefficients["insomnia", "Estimate"],
    summary(fit_psm_cov_adj_mci)$coefficients["insomnia", "Estimate"],
    summary(fit_psm_cov_adj_svy_osw_mci)$coefficients["insomnia", "Estimate"],
    summary(fit_psm_cov_adj_svy_isw_mci)$coefficients["insomnia", "Estimate"]
  )
  # SE results
  psm_se_mci[rep, ] <- c(
    summary(fit_psm_weighted_unadj_mci)$coefficients["insomnia", "Std. Error"],
    summary(fit_psm_weighted_unadj_svy_osw_mci)$coefficients["insomnia", "Std. Error"],
    summary(fit_psm_weighted_unadj_svy_isw_mci)$coefficients["insomnia", "Std. Error"],
    summary(fit_psm_weighted_adj_mci)$coefficients["insomnia", "Std. Error"],
    summary(fit_psm_weighted_adj_svy_osw_mci)$coefficients["insomnia", "Std. Error"],
    summary(fit_psm_weighted_adj_svy_isw_mci)$coefficients["insomnia", "Std. Error"],
    
    summary(fit_psm_cov_unadj_mci)$coefficients["insomnia", "Std. Error"],
    summary(fit_psm_cov_unadj_svy_osw_mci)$coefficients["insomnia", "Std. Error"],
    summary(fit_psm_cov_unadj_svy_isw_mci)$coefficients["insomnia", "Std. Error"],
    summary(fit_psm_cov_adj_mci)$coefficients["insomnia", "Std. Error"],
    summary(fit_psm_cov_adj_svy_osw_mci)$coefficients["insomnia", "Std. Error"],
    summary(fit_psm_cov_adj_svy_isw_mci)$coefficients["insomnia", "Std. Error"]
  )
  
}

# look at results
head(psm_effect_mci)
head(psm_se_mci)

# remove unneeded objects
rm(rep, list = ls(pattern = "fit"), psm_analyses_mci)

### CONDUCT ANALYSES HYP -------------------------------------------------------------

# list the type of analyses we will conduct for hypertension for PSM
psm_analyses_hyp <- c(
  "psm_weighted_unadj_hyp", "psm_weighted_unadj_svy_osw_hyp", "psm_weighted_unadj_svy_isw_hyp",
  "psm_weighted_adj_hyp", "psm_weighted_adj_svy_osw_hyp", "psm_weighted_adj_svy_isw_hyp",
  "psm_cov_unadj_hyp", "psm_cov_unadj_svy_osw_hyp", "psm_cov_unadj_svy_isw_hyp",
  "psm_cov_adj_hyp", "psm_cov_adj_svy_osw_hyp", "psm_cov_adj_svy_isw_hyp"
)
# create matrix to hold effect estimates from the various analyses
psm_effect_hyp <- matrix(NA, nrow = reps, ncol = length(psm_analyses_hyp), dimnames = list(NULL, psm_analyses_hyp))
# create matrix to hold SEs from the various analyses
psm_se_hyp <- matrix(NA, nrow = reps, ncol = length(psm_analyses_hyp), dimnames = list(NULL, psm_analyses_hyp))

# loop through samples to conduct the various analyses on the matched data
for (rep in 1:reps) {
  if((rep %% 100) == 0) print(rep)
  
  # fit the various regressions
  
  # when propensity score is calculated using weighted logistic regression 
  # unadjusted regression without weights
  fit_psm_weighted_unadj_hyp <- glm(
    hyp_v2 ~ insomnia + offset(log(visit_years)), 
    family = poisson(), 
    data = list_matches_psm_weighted[[rep]] %>% filter(hyp_v1 == 0)
  )
  # survey unadjusted regression using OSW
  fit_psm_weighted_unadj_svy_osw_hyp <- svyglm(
    hyp_v2 ~ insomnia + offset(log(visit_years)), 
    family = poisson(), 
    design = list_matches_psm_weighted_osw_svy[[rep]], 
    subset = hyp_v1 == 0
  )
  # survey unadjusted regression using ISW
  fit_psm_weighted_unadj_svy_isw_hyp <- svyglm(
    hyp_v2 ~ insomnia + offset(log(visit_years)), 
    family = poisson(),  
    design = list_matches_psm_weighted_isw_svy[[rep]],
    subset = hyp_v1 == 0
  )
  # adjusted regression without weights
  fit_psm_weighted_adj_hyp <- glm(
    hyp_v2 ~ insomnia + age + offset(log(visit_years)), 
    family = poisson(), 
    data = list_matches_psm_weighted[[rep]] %>% filter(hyp_v1 == 0),
  )
  # survey adjusted regression using OSW
  fit_psm_weighted_adj_svy_osw_hyp <- svyglm(
    hyp_v2 ~ insomnia + age + offset(log(visit_years)), 
    family = poisson(), 
    design = list_matches_psm_weighted_osw_svy[[rep]],
    subset = hyp_v1 == 0
  )
  # survey adjusted regression using ISW
  fit_psm_weighted_adj_svy_isw_hyp <- svyglm(
    hyp_v2 ~ insomnia + age + offset(log(visit_years)), 
    family = poisson(), 
    design = list_matches_psm_weighted_isw_svy[[rep]],
    subset = hyp_v1 == 0
  )
  
  # when propensity score is calculated using logistic regression with weight as covariate
  # unadjusted regression without weights
  fit_psm_cov_unadj_hyp <- glm(
    hyp_v2 ~ insomnia + offset(log(visit_years)), 
    family = poisson(), 
    data = list_matches_psm_cov[[rep]] %>% filter(hyp_v1 == 0)
  )
  # survey unadjusted regression using OSW
  fit_psm_cov_unadj_svy_osw_hyp <- svyglm(
    hyp_v2 ~ insomnia + offset(log(visit_years)), 
    family = poisson(), 
    design = list_matches_psm_cov_osw_svy[[rep]], 
    subset = hyp_v1 == 0
  )
  # survey unadjusted regression using ISW
  fit_psm_cov_unadj_svy_isw_hyp <- svyglm(
    hyp_v2 ~ insomnia + offset(log(visit_years)), 
    family = poisson(),  
    design = list_matches_psm_cov_isw_svy[[rep]],
    subset = hyp_v1 == 0
  )
  # adjusted regression without weights
  fit_psm_cov_adj_hyp <- glm(
    hyp_v2 ~ insomnia + age + offset(log(visit_years)), 
    family = poisson(), 
    data = list_matches_psm_cov[[rep]] %>% filter(hyp_v1 == 0),
  )
  # survey adjusted regression using OSW
  fit_psm_cov_adj_svy_osw_hyp <- svyglm(
    hyp_v2 ~ insomnia + age + offset(log(visit_years)), 
    family = poisson(), 
    design = list_matches_psm_cov_osw_svy[[rep]],
    subset = hyp_v1 == 0
  )
  # survey adjusted regression using ISW
  fit_psm_cov_adj_svy_isw_hyp <- svyglm(
    hyp_v2 ~ insomnia + age + offset(log(visit_years)), 
    family = poisson(), 
    design = list_matches_psm_cov_isw_svy[[rep]],
    subset = hyp_v1 == 0
  )
  
  # add result to associated matrix
  # effect estimate results
  psm_effect_hyp[rep, ] <- c(
    summary(fit_psm_weighted_unadj_hyp)$coefficients["insomnia", "Estimate"],
    summary(fit_psm_weighted_unadj_svy_osw_hyp)$coefficients["insomnia", "Estimate"],
    summary(fit_psm_weighted_unadj_svy_isw_hyp)$coefficients["insomnia", "Estimate"],
    summary(fit_psm_weighted_adj_hyp)$coefficients["insomnia", "Estimate"],
    summary(fit_psm_weighted_adj_svy_osw_hyp)$coefficients["insomnia", "Estimate"],
    summary(fit_psm_weighted_adj_svy_isw_hyp)$coefficients["insomnia", "Estimate"],
    
    summary(fit_psm_cov_unadj_hyp)$coefficients["insomnia", "Estimate"],
    summary(fit_psm_cov_unadj_svy_osw_hyp)$coefficients["insomnia", "Estimate"],
    summary(fit_psm_cov_unadj_svy_isw_hyp)$coefficients["insomnia", "Estimate"],
    summary(fit_psm_cov_adj_hyp)$coefficients["insomnia", "Estimate"],
    summary(fit_psm_cov_adj_svy_osw_hyp)$coefficients["insomnia", "Estimate"],
    summary(fit_psm_cov_adj_svy_isw_hyp)$coefficients["insomnia", "Estimate"]
  )
  # SE results
  psm_se_hyp[rep, ] <- c(
    summary(fit_psm_weighted_unadj_hyp)$coefficients["insomnia", "Std. Error"],
    summary(fit_psm_weighted_unadj_svy_osw_hyp)$coefficients["insomnia", "Std. Error"],
    summary(fit_psm_weighted_unadj_svy_isw_hyp)$coefficients["insomnia", "Std. Error"],
    summary(fit_psm_weighted_adj_hyp)$coefficients["insomnia", "Std. Error"],
    summary(fit_psm_weighted_adj_svy_osw_hyp)$coefficients["insomnia", "Std. Error"],
    summary(fit_psm_weighted_adj_svy_isw_hyp)$coefficients["insomnia", "Std. Error"],
    
    summary(fit_psm_cov_unadj_hyp)$coefficients["insomnia", "Std. Error"],
    summary(fit_psm_cov_unadj_svy_osw_hyp)$coefficients["insomnia", "Std. Error"],
    summary(fit_psm_cov_unadj_svy_isw_hyp)$coefficients["insomnia", "Std. Error"],
    summary(fit_psm_cov_adj_hyp)$coefficients["insomnia", "Std. Error"],
    summary(fit_psm_cov_adj_svy_osw_hyp)$coefficients["insomnia", "Std. Error"],
    summary(fit_psm_cov_adj_svy_isw_hyp)$coefficients["insomnia", "Std. Error"]
  )
  
}

# look at results
head(psm_effect_hyp)
head(psm_se_hyp)

# remove unneeded objects
rm(rep, list = ls(pattern = "fit"), psm_analyses_hyp)

### CALCULATE BALANCE ----------------------------------------------------------------

# create function to calculate standardized mean difference (SMD) in a sample for PSM
calc_smd_psm <- function(var, rep) {
  # means by exposure group using different weighting schemes for weighting
  means <- tibble(
    
    # unweighted mean 
    # where propensity score is calculated using weighted logistic regression 
    mean_psm_weighted_unweighted = list_matches_psm_weighted[[rep]] %>% 
      select({{ var }}, insomnia, final_weights) %>%
      group_by(insomnia) %>%
      summarise(mean = mean({{ var }})) %>%
      ungroup() %>%
      pull(mean), 
    # where propensity score is calculated using logistic regression with weight as covariate
    mean_psm_cov_unweighted = list_matches_psm_cov[[rep]] %>% 
      select({{ var }}, insomnia, final_weights) %>%
      group_by(insomnia) %>%
      summarise(mean = mean({{ var }})) %>%
      ungroup() %>%
      pull(mean),
    
    
    # weighted mean
    
    # where propensity score is calculated using weighted logistic regression 
    # using OSW
    mean_psm_weighted_osw = svyby(~ {{ var }}, by = ~ insomnia, list_matches_psm_weighted_osw_svy[[rep]], svymean) %>% pull({{ var }}),
    # using ISW 
    mean_psm_weighted_isw = svyby(~ {{ var }}, by = ~ insomnia, list_matches_psm_weighted_isw_svy[[rep]], svymean) %>% pull({{ var }}),
    
    # where propensity score is calculated using logistic regression with weight as covariate
    # using OSW
    mean_psm_cov_osw = svyby(~ {{ var }}, by = ~ insomnia, list_matches_psm_cov_osw_svy[[rep]], svymean) %>% pull({{ var }}),
    # using ISW 
    mean_psm_cov_isw = svyby(~ {{ var }}, by = ~ insomnia, list_matches_psm_cov_isw_svy[[rep]], svymean) %>% pull({{ var }}),
    
  )
  
  # sd weighted by OSW in exposed group (before matching sd)
  sd_osw <- tibble(sd_osw = svyby(~ {{ var }}, by = ~ insomnia, list_samps_svy[[rep]], svyvar) %>% pull({{ var }}) %>% sqrt(.)) %>% pull(sd_osw) %>% .[2]
  
  # calculate smd
  # calculate (absolute) difference in means between exposure groups
  diff_mean <- abs(unlist(means[2, ] - means[1, ]))
  # calculate smd by dividing difference in means by before matching sd
  smd <- diff_mean / sd_osw
  
  # return result
  return(smd)
}

# create function to create data frame with smd calculations for all covariates in a sample and add column with global smd
calc_smd_all_global_psm <- function(rep) {
  # calculate smd for all three covariates for each method of weights to calculate mean
  smd_all_global <- data.frame(
    age_smd = calc_smd_psm(age, rep), 
    bmi_smd = calc_smd_psm(bmi, rep), 
    visit_years_smd = calc_smd_psm(visit_years, rep)
  ) %>% 
    # calculate the mean smd of the three variables
    rowwise(.) %>%
    mutate(global_smd = (age_smd + bmi_smd + visit_years_smd) / 3) %>%
    # translate to data.frame object to allow for row names
    data.frame(row.names = names(calc_smd_psm(age, rep)))
  
  # return result
  return(smd_all_global)
}

# create empty list to store smd calculations for each sample
list_balance_psm <- vector(mode = "list", length = reps)

# loop through the samples to obtain the smd calculations for each sample
for (rep in 1:reps) {
  if((rep %% 100) == 0) print(rep)
  list_balance_psm[[rep]] <- calc_smd_all_global_psm(rep)
}

# find the average balance metrics across the samples
balance_avg_psm <- Reduce(`+`, list_balance_psm) / reps
balance_avg_psm

# remove unneeded objects
rm(rep, calc_smd_psm, calc_smd_all_global_psm)

### SAVE RESULTS ---------------------------------------------------------------------

# create list to hold effect estimate and SE results for the two outcomes
# where effect is given on the OR scale
list_psm_effect_or <- list(mci = exp(psm_effect_mci), hyp = exp(psm_effect_hyp))
list_psm_se <- list(mci = psm_se_mci, hyp = psm_se_hyp)

# save effect estimate and SE results for the two outcomes
save(list_psm_effect_or, list_psm_se, file = paste0("./data/simulation_results/confounding_sensitivity_results/remove_bmi_results/psm_analyses_no_bmi", specify_df, ".RData")) # adjust path to specify which simulated data frame was used

# save balance metrics
save(balance_avg_psm, file = paste0("./data/simulation_results/confounding_sensitivity_results/remove_bmi_results/balance_avg_psm_no_bmi", specify_df, ".RData")) # adjust path to specify which simulated data frame was used



