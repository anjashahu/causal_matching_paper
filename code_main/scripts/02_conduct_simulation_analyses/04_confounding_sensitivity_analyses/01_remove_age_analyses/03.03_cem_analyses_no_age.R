
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

# note: for this sensitivity analysis, we are pretending that we do not have age information
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

### CONDUCT CEM --------------------------------------------------------------

# run R script to obtain function to conduct matching
source("./scripts/02_conduct_simulation_analyses/04_confounding_sensitivity_analyses/01_remove_age_analyses/02.02_conduct_cem_function_no_age.R")

# create empty lists to store matches using coarsened exact matching (CEM)
# when CEM is based on binning of coarsened covariates
list_matches_cem_cov <- vector(mode = "list", length = reps)
# when CEM is based on binning of coarsened covariates + coarsened OSW
list_matches_cem_weighted <- vector(mode = "list", length = reps)
# loop through samples to conduct cem for each sample
for (rep in 1:reps){
  if((rep %% 100) == 0) print(rep)
  cem_rep_results <- conduct_cem(list_samps[[rep]])
  list_matches_cem_cov[[rep]] <- cem_rep_results[[1]]
  list_matches_cem_weighted[[rep]] <- cem_rep_results[[2]]
}

# remove unneeded objects
rm(rep, conduct_cem, cem_rep_results)

### CREATE SURVEY DESIGN OBJECTS -------------------------------------------------

# for CEM

# create empty list to store survey objects of CEM matches that use CEMW
# where bins are defined by coarsened covariates
list_matches_cem_cov_cemw_svy <- vector(mode = "list", length = reps)
# where bins are defined by coarsened covariates + coarsened OSW
list_matches_cem_weighted_cemw_svy <- vector(mode = "list", length = reps)

# create empty list to store survey objects of CEM matches that use CEMW x OSW
# where bins are defined by coarsened covariates
list_matches_cem_cov_cemw_osw_svy <- vector(mode = "list", length = reps)
# where bins are defined by coarsened covariates + coarsened OSW
list_matches_cem_weighted_cemw_osw_svy <- vector(mode = "list", length = reps)

# loop through the CEM matches, create survey design objects with different weights and add to associated list
for (rep in 1:reps) {
  if((rep %% 100) == 0) print(rep)
  
  # for CEM matches using CEMW
  # where bins are defined by coarsened covariates
  list_matches_cem_cov_cemw_svy[[rep]] <- svydesign(
    ids = ~ bg + hh, 
    weights = ~ cem_weights, 
    strata = ~ strata, 
    data = list_matches_cem_cov[[rep]]
  )
  # where bins are defined by coarsened covariates + coarsened OSW
  list_matches_cem_weighted_cemw_svy[[rep]] <- svydesign(
    ids = ~ bg + hh, 
    weights = ~ cem_weights, 
    strata = ~ strata, 
    data = list_matches_cem_weighted[[rep]]
  )
  
  # for CEM matches using CEMW x OSW
  list_matches_cem_cov_cemw_osw_svy[[rep]] <- svydesign(
    ids = ~ bg + hh, 
    weights = ~ cem_osw_weights, 
    strata = ~ strata, 
    data = list_matches_cem_cov[[rep]]
  )
  list_matches_cem_weighted_cemw_osw_svy[[rep]] <- svydesign(
    ids = ~ bg + hh, 
    weights = ~ cem_osw_weights, 
    strata = ~ strata, 
    data = list_matches_cem_weighted[[rep]]
  )
}

# remove unneeded variables
rm(rep)

### CONDUCT ANALYSES MCI ------------------------------------------------------

# list the type of analyses we will conduct for MCI for CEM
cem_analyses_mci <- c(
  "cem_cov_cemw_unadj_svy_mci", "cem_cov_cemw_osw_unadj_svy_mci", 
  "cem_cov_cemw_adj_svy_mci", "cem_cov_cemw_osw_adj_svy_mci",
  "cem_weighted_cemw_unadj_svy_mci", "cem_weighted_cemw_osw_unadj_svy_mci", 
  "cem_weighted_cemw_adj_svy_mci", "cem_weighted_cemw_osw_adj_svy_mci"
)
# create matrix to hold effect estimates from the various analyses
cem_effect_mci <- matrix(NA, nrow = reps, ncol = length(cem_analyses_mci), dimnames = list(NULL, cem_analyses_mci))
# create matrix to hold SEs from the various analyses
cem_se_mci <- matrix(NA, nrow = reps, ncol = length(cem_analyses_mci), dimnames = list(NULL, cem_analyses_mci))

# loop through samples to conduct the various analyses on the full sample
for (rep in 1:reps) {
  if((rep %% 100) == 0) print(rep)
  
  # fit the various regressions
  
  # when bins are defined by coarsened covariates
  # unadjusted survey regression with CEMW
  fit_cem_cov_cemw_unadj_svy_mci <- svyglm(
    mci ~ insomnia, 
    family = quasibinomial(), 
    design = list_matches_cem_cov_cemw_svy[[rep]]
  )
  # unadjusted survey regression with CEMW x OSW
  fit_cem_cov_cemw_osw_unadj_svy_mci <- svyglm(
    mci ~ insomnia, 
    family = quasibinomial(), 
    design = list_matches_cem_cov_cemw_osw_svy[[rep]]
  )
  # adjusted survey regression with CEMW
  fit_cem_cov_cemw_adj_svy_mci <- svyglm(
    mci ~ insomnia + bmi, 
    family = quasibinomial(), 
    design = list_matches_cem_cov_cemw_svy[[rep]]
  )
  # adjusted survey regression with CEMW x OSW
  fit_cem_cov_cemw_osw_adj_svy_mci <- svyglm(
    mci ~ insomnia + bmi, 
    family = quasibinomial(), 
    design = list_matches_cem_cov_cemw_osw_svy[[rep]]
  )
  
  # when bins are defined by coarsened covariates + coarsened OSW
  # unadjusted survey regression with CEMW
  fit_cem_weighted_cemw_unadj_svy_mci <- svyglm(
    mci ~ insomnia, 
    family = quasibinomial(), 
    design = list_matches_cem_weighted_cemw_svy[[rep]]
  )
  # unadjusted survey regression with CEMW x OSW
  fit_cem_weighted_cemw_osw_unadj_svy_mci <- svyglm(
    mci ~ insomnia, 
    family = quasibinomial(), 
    design = list_matches_cem_weighted_cemw_osw_svy[[rep]]
  )
  # adjusted survey regression with CEMW
  fit_cem_weighted_cemw_adj_svy_mci <- svyglm(
    mci ~ insomnia + bmi, 
    family = quasibinomial(), 
    design = list_matches_cem_weighted_cemw_svy[[rep]]
  )
  # adjusted survey regression with CEMW x OSW
  fit_cem_weighted_cemw_osw_adj_svy_mci <- svyglm(
    mci ~ insomnia + bmi, 
    family = quasibinomial(), 
    design = list_matches_cem_weighted_cemw_osw_svy[[rep]]
  )
  
  # add result to associated matrix
  # effect estimate results
  cem_effect_mci[rep, ] <- c(
    summary(fit_cem_cov_cemw_unadj_svy_mci)$coefficients["insomnia", "Estimate"],
    summary(fit_cem_cov_cemw_osw_unadj_svy_mci)$coefficients["insomnia", "Estimate"],
    summary(fit_cem_cov_cemw_adj_svy_mci)$coefficients["insomnia", "Estimate"],
    summary(fit_cem_cov_cemw_osw_adj_svy_mci)$coefficients["insomnia", "Estimate"],
    
    summary(fit_cem_weighted_cemw_unadj_svy_mci)$coefficients["insomnia", "Estimate"],
    summary(fit_cem_weighted_cemw_osw_unadj_svy_mci)$coefficients["insomnia", "Estimate"],
    summary(fit_cem_weighted_cemw_adj_svy_mci)$coefficients["insomnia", "Estimate"],
    summary(fit_cem_weighted_cemw_osw_adj_svy_mci)$coefficients["insomnia", "Estimate"]
  )
  
  cem_se_mci[rep, ] <- c(
    summary(fit_cem_cov_cemw_unadj_svy_mci)$coefficients["insomnia", "Std. Error"],
    summary(fit_cem_cov_cemw_osw_unadj_svy_mci)$coefficients["insomnia", "Std. Error"],
    summary(fit_cem_cov_cemw_adj_svy_mci)$coefficients["insomnia", "Std. Error"],
    summary(fit_cem_cov_cemw_osw_adj_svy_mci)$coefficients["insomnia", "Std. Error"],
    
    summary(fit_cem_weighted_cemw_unadj_svy_mci)$coefficients["insomnia", "Std. Error"],
    summary(fit_cem_weighted_cemw_osw_unadj_svy_mci)$coefficients["insomnia", "Std. Error"],
    summary(fit_cem_weighted_cemw_adj_svy_mci)$coefficients["insomnia", "Std. Error"],
    summary(fit_cem_weighted_cemw_osw_adj_svy_mci)$coefficients["insomnia", "Std. Error"]
  )
}

# look at results
head(cem_effect_mci)
head(cem_se_mci)

# remove unneeded objects
rm(rep, list = ls(pattern = "fit"), cem_analyses_mci)

### CONDUCT ANALYSES HYP ------------------------------------------------------

# list the type of analyses we will conduct for hypertension for CEM
cem_analyses_hyp <- c(
  "cem_cov_cemw_unadj_svy_hyp", "cem_cov_cemw_osw_unadj_svy_hyp", 
  "cem_cov_cemw_adj_svy_hyp", "cem_cov_cemw_osw_adj_svy_hyp",
  "cem_weighted_cemw_unadj_svy_hyp", "cem_weighted_cemw_osw_unadj_svy_hyp", 
  "cem_weighted_cemw_adj_svy_hyp", "cem_weighted_cemw_osw_adj_svy_hyp"
)
# create matrix to hold effect estimates from the various analyses
cem_effect_hyp <- matrix(NA, nrow = reps, ncol = length(cem_analyses_hyp), dimnames = list(NULL, cem_analyses_hyp))
# create matrix to hold SEs from the various analyses
cem_se_hyp <- matrix(NA, nrow = reps, ncol = length(cem_analyses_hyp), dimnames = list(NULL, cem_analyses_hyp))

# loop through samples to conduct the various analyses on the full sample
for (rep in 1:reps) {
  if((rep %% 100) == 0) print(rep)
  
  # fit the various regressions
  
  # when bins are defined by coarsened covariates
  # unadjusted survey regression with CEMW
  fit_cem_cov_cemw_unadj_svy_hyp <- svyglm(
    hyp_v2 ~ insomnia + offset(log(visit_years)), 
    family = poisson(), 
    design = list_matches_cem_cov_cemw_svy[[rep]],
    subset = hyp_v1 == 0
  )
  
  # unadjusted survey regression with CEMW x OSW
  fit_cem_cov_cemw_osw_unadj_svy_hyp <- svyglm(
    hyp_v2 ~ insomnia + offset(log(visit_years)), 
    family = poisson(),  
    design = list_matches_cem_cov_cemw_osw_svy[[rep]],
    subset = hyp_v1 == 0
  )
  # adjusted survey regression with CEMW
  fit_cem_cov_cemw_adj_svy_hyp <- svyglm(
    hyp_v2 ~ insomnia + bmi + offset(log(visit_years)), 
    family = poisson(), 
    design = list_matches_cem_cov_cemw_svy[[rep]],
    subset = hyp_v1 == 0
  )
  # adjusted survey regression with CEMW x OSW
  fit_cem_cov_cemw_osw_adj_svy_hyp <- svyglm(
    hyp_v2 ~ insomnia + bmi + offset(log(visit_years)), 
    family = poisson(),  
    design = list_matches_cem_cov_cemw_osw_svy[[rep]],
    subset = hyp_v1 == 0
  )
  
  # when bins are defined by coarsened covariates + coarsened OSW
  # unadjusted survey regression with CEMW
  fit_cem_weighted_cemw_unadj_svy_hyp <- svyglm(
    hyp_v2 ~ insomnia + offset(log(visit_years)), 
    family = poisson(),  
    design = list_matches_cem_weighted_cemw_svy[[rep]],
    subset = hyp_v1 == 0
  )
  # unadjusted survey regression with CEMW x OSW
  fit_cem_weighted_cemw_osw_unadj_svy_hyp <- svyglm(
    hyp_v2 ~ insomnia + offset(log(visit_years)), 
    family = poisson(),  
    design = list_matches_cem_weighted_cemw_osw_svy[[rep]],
    subset = hyp_v1 == 0
  )
  # adjusted survey regression with CEMW
  fit_cem_weighted_cemw_adj_svy_hyp <- svyglm(
    hyp_v2 ~ insomnia + bmi + offset(log(visit_years)), 
    family = poisson(), 
    design = list_matches_cem_weighted_cemw_svy[[rep]],
    subset = hyp_v1 == 0
  )
  # adjusted survey regression with CEMW x OSW
  fit_cem_weighted_cemw_osw_adj_svy_hyp <- svyglm(
    hyp_v2 ~ insomnia + bmi + offset(log(visit_years)), 
    family = poisson(), 
    design = list_matches_cem_weighted_cemw_osw_svy[[rep]],
    subset = hyp_v1 == 0
  )
  
  
  # add result to associated matrix
  # effect estimate results
  cem_effect_hyp[rep, ] <- c(
    summary(fit_cem_cov_cemw_unadj_svy_hyp)$coefficients["insomnia", "Estimate"],
    summary(fit_cem_cov_cemw_osw_unadj_svy_hyp)$coefficients["insomnia", "Estimate"],
    summary(fit_cem_cov_cemw_adj_svy_hyp)$coefficients["insomnia", "Estimate"],
    summary(fit_cem_cov_cemw_osw_adj_svy_hyp)$coefficients["insomnia", "Estimate"],
    
    summary(fit_cem_weighted_cemw_unadj_svy_hyp)$coefficients["insomnia", "Estimate"],
    summary(fit_cem_weighted_cemw_osw_unadj_svy_hyp)$coefficients["insomnia", "Estimate"],
    summary(fit_cem_weighted_cemw_adj_svy_hyp)$coefficients["insomnia", "Estimate"],
    summary(fit_cem_weighted_cemw_osw_adj_svy_hyp)$coefficients["insomnia", "Estimate"]
  )
  
  cem_se_hyp[rep, ] <- c(
    summary(fit_cem_cov_cemw_unadj_svy_hyp)$coefficients["insomnia", "Std. Error"],
    summary(fit_cem_cov_cemw_osw_unadj_svy_hyp)$coefficients["insomnia", "Std. Error"],
    summary(fit_cem_cov_cemw_adj_svy_hyp)$coefficients["insomnia", "Std. Error"],
    summary(fit_cem_cov_cemw_osw_adj_svy_hyp)$coefficients["insomnia", "Std. Error"],
    
    summary(fit_cem_weighted_cemw_unadj_svy_hyp)$coefficients["insomnia", "Std. Error"],
    summary(fit_cem_weighted_cemw_osw_unadj_svy_hyp)$coefficients["insomnia", "Std. Error"],
    summary(fit_cem_weighted_cemw_adj_svy_hyp)$coefficients["insomnia", "Std. Error"],
    summary(fit_cem_weighted_cemw_osw_adj_svy_hyp)$coefficients["insomnia", "Std. Error"]
  )
}

# look at results
head(cem_effect_hyp)
head(cem_se_hyp)

# remove unneeded objects
rm(rep, list = ls(pattern = "fit"), cem_analyses_hyp)

### CALCULATE BALANCE ----------------------------------------------------------------

# create function to calculate standardized mean difference (SMD) in a sample for CEM
calc_smd_cem <- function(var, rep) {
  # means by exposure group using different weighting schemes for CEM
  means <- tibble(
    # mean weighted by CEMW
    # where bins are defined by coarsened covariates
    mean_cem_cov_cemw = svyby(~ {{ var }}, by = ~ insomnia, list_matches_cem_cov_cemw_svy[[rep]], svymean) %>% pull({{ var }}),
    # where bins are defined by coarsened covariates + coarsened OSW
    mean_cem_weighted_cemw = svyby(~ {{ var }}, by = ~ insomnia, list_matches_cem_weighted_cemw_svy[[rep]], svymean) %>% pull({{ var }}),
    
    # mean weighted by CEMW x OSW
    # where bins are defined by coarsened covariates
    mean_cem_cov_cemw_osw = svyby(~ {{ var }}, by = ~ insomnia, list_matches_cem_cov_cemw_osw_svy[[rep]], svymean) %>% pull({{ var }}),
    # where bins are defined by coarsened covariates + coarsened OSW
    mean_cem_weighted_cemw_osw = svyby(~ {{ var }}, by = ~ insomnia, list_matches_cem_weighted_cemw_osw_svy[[rep]], svymean) %>% pull({{ var }}),
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
calc_smd_all_global_cem <- function(rep) {
  # calculate smd for all three covariates for each method of using weights to calculate the means
  smd_all_global <- data.frame(
    age_smd = calc_smd_cem(age, rep), 
    bmi_smd = calc_smd_cem(bmi, rep), 
    visit_years_smd = calc_smd_cem(visit_years, rep)
  ) %>% 
    # calculate the mean smd of the three variables 
    rowwise(.) %>%
    mutate(global_smd = (age_smd + bmi_smd + visit_years_smd) / 3) %>%
    # translate to data.frame object to allow for row names
    data.frame(row.names = names(calc_smd_cem(age, rep)))
  
  # return result
  return(smd_all_global)
}

# create empty list to store smd calculations for each sample
list_balance_cem <- vector(mode = "list", length = reps)

# loop through the samples to obtain the smd calculations for each sample
for (rep in 1:reps) {
  if((rep %% 100) == 0) print(rep)
  list_balance_cem[[rep]] <- calc_smd_all_global_cem(rep)
}

# find the average balance metrics across the samples for weighting
balance_avg_cem <- Reduce(`+`, list_balance_cem) / reps
balance_avg_cem

# remove unneeded objects
rm(rep, calc_smd_cem, calc_smd_all_global_cem)

### SAVE RESULTS ---------------------------------------------------------------------

# create list to hold effect estimate and SE results for the two outcomes
# where effect is given on the OR scale
list_cem_effect_or <- list(mci = exp(cem_effect_mci), hyp = exp(cem_effect_hyp))
list_cem_se <- list(mci = cem_se_mci, hyp = cem_se_hyp)

# save effect estimate and SE results for the two outcomes
save(list_cem_effect_or, list_cem_se, file = paste0("./data/simulation_results/confounding_sensitivity_results/remove_age_results/cem_analyses_no_age", specify_df, ".RData")) # adjust path to specify which simulated data frame was used

# save balance metrics
save(balance_avg_cem, file = paste0("./data/simulation_results/confounding_sensitivity_results/remove_age_results/balance_avg_cem_no_age", specify_df, ".RData")) # adjust path to specify which simulated data frame was used



