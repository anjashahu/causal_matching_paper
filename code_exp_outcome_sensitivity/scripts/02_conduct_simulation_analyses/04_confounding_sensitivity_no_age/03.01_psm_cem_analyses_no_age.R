
### SET UP ---------------------------------------------------------------------

# load packages
library(tidyverse)
library(survey)

# clear working environment
rm(list = ls())

# specify which of the individual data frames you want to use
df_num <- 1
df_num

# load data 
df_pop_bg <- read_csv("./data/simulated_pop_data/df_pop_bg.csv")
df_pop_hh <- read_csv("./data/simulated_pop_data/df_pop_hh.csv")
df_pop_indiv <- read_csv(paste0("./data/simulated_pop_data/df_pop_indiv/df_pop_indiv_", df_num, ".csv")) 

# set numbers of repetitions
reps <- 1000
reps

# note: we have removed age from this analysis

### SIMULATE SAMPLES --------------------------------------------------------
# run R script to obtain function to simulate samples
source("./scripts/02_conduct_simulation_analyses/01_simulate_samples_function.R")

# create empty list to store samples
list_samps <- vector(mode = "list", length = reps)

# simulate all the samples and their survey design object by looping through the number of repetitions and storing each sample into list
for (rep in 1:reps) {
  if((rep %% 100) == 0) print(rep)
  # simulate sample for the repetition and add to associated list
  list_samps[[rep]] <- sim_one_sample(df_pop_bg, df_pop_hh, df_pop_indiv, rep) 
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

### CONDUCT PSM -----------------------------------------------------------

# run R script to obtain function to conduct matching
source("./scripts/02_conduct_simulation_analyses/04_confounding_sensitivity_no_age/02.01_conduct_psm_function_no_age.R")

# create empty lists to store matches using propensity scores 
# where propensity score is based on logistic regression with OSW as covariate
list_matches_psm_cov <- vector(mode = "list", length = reps)
# loop through samples to obtain matches 
for (rep in 1:reps){
  if((rep %% 100) == 0) print(rep)
  list_matches_psm_cov[[rep]] <- conduct_one_psm(list_samps[[rep]])
}

# remove unneeded objects
rm(rep, get_inherited_weights, conduct_one_psm)

### CONDUCT CEM --------------------------------------------------------------

# run R script to obtain function to conduct matching
source("./scripts/02_conduct_simulation_analyses/04_confounding_sensitivity_no_age/02.02_conduct_cem_function_no_age.R")

# create empty lists to store matches using coarsened exact matching (CEM)
# when CEM is based on binning of coarsened covariates + coarsened OSW
list_matches_cem_weighted <- vector(mode = "list", length = reps)
# loop through samples to conduct cem for each sample
for (rep in 1:reps){
  if((rep %% 100) == 0) print(rep)
  list_matches_cem_weighted[[rep]] <- conduct_cem(list_samps[[rep]])
}

# remove unneeded objects
rm(rep, conduct_cem, list_samps)

### CREATE SURVEY DESIGN OBJECTS -------------------------------------------------

# for PSM 
# create empty list to store survey objects of PSM matches that use OSW
# where propensity score is calculated using logistic regression with weight as covariate
list_matches_psm_cov_osw_svy <- vector(mode = "list", length = reps)

# for CEM
# create empty list to store survey objects of CEM matches that use CEMW x OSW
# where bins are defined by coarsened covariates + coarsened OSW
list_matches_cem_weighted_cemw_osw_svy <- vector(mode = "list", length = reps)

# loop through the matches, create survey design objects with different weights and add to associated list
for (rep in 1:reps) {
  if((rep %% 100) == 0) print(rep)
  
  # for PSM matches using OSW
  # where propensity score is calculated using logistic regression with weight as covariate
  list_matches_psm_cov_osw_svy[[rep]] <- svydesign(
    ids = ~ bg + hh, 
    weights = ~ final_weights, 
    strata = ~ strata, 
    data = list_matches_psm_cov[[rep]]
  )
  
  # for CEM matches using CEMW x OSW
  # where matching is based on coarsened covariates + coarsened weights
  list_matches_cem_weighted_cemw_osw_svy[[rep]] <- svydesign(
    ids = ~ bg + hh, 
    weights = ~ cem_osw_weights, 
    strata = ~ strata, 
    data = list_matches_cem_weighted[[rep]]
  )
  
}

# remove unneeded variables
rm(rep)

### CONDUCT ANALYSES MCI -------------------------------------------------------------

# list the type of analyses we will conduct for MCI
analyses_mci <- c("psm_cov_unadj_svy_osw_mci", "psm_cov_adj_svy_osw_mci",
                  "cem_weighted_cemw_osw_unadj_svy_mci", "cem_weighted_cemw_osw_adj_svy_mci")
# create matrix to hold effect estimates from the various analyses
effect_mci <- matrix(NA, nrow = reps, ncol = length(analyses_mci), dimnames = list(NULL, analyses_mci))
# create matrix to hold SEs from the various analyses
se_mci <- matrix(NA, nrow = reps, ncol = length(analyses_mci), dimnames = list(NULL, analyses_mci))

# loop through samples to conduct the various analyses on the matched data
for (rep in 1:reps) {
  if((rep %% 100) == 0) print(rep)
  
  # fit the various regressions 
  
  # using PSM data
  # when propensity score is calculated using logistic regression with weight as covariate
  # survey unadjusted regression using OSW
  fit_psm_cov_unadj_svy_osw_mci <- svyglm(
    mci ~ insomnia, 
    family = quasibinomial(), 
    design = list_matches_psm_cov_osw_svy[[rep]]
  )
  # survey adjusted regression using OSW
  fit_psm_cov_adj_svy_osw_mci <- svyglm(
    mci ~ insomnia + bmi, 
    family = quasibinomial(), 
    design = list_matches_psm_cov_osw_svy[[rep]]
  )
  
  # using the CEM data
  # when bins are defined by coarsened covariates + coarsened OSW
  # unadjusted survey regression with CEMW x OSW
  fit_cem_weighted_cemw_osw_unadj_svy_mci <- svyglm(
    mci ~ insomnia, 
    family = quasibinomial(), 
    design = list_matches_cem_weighted_cemw_osw_svy[[rep]]
  )
  # adjusted survey regression with CEMW x OSW
  fit_cem_weighted_cemw_osw_adj_svy_mci <- svyglm(
    mci ~ insomnia + bmi, 
    family = quasibinomial(), 
    design = list_matches_cem_weighted_cemw_osw_svy[[rep]]
  )
  
  # add results to associated matrix
  # effect estimate results
  effect_mci[rep, ] <- c(
    summary(fit_psm_cov_unadj_svy_osw_mci)$coefficients["insomnia", "Estimate"],
    summary(fit_psm_cov_adj_svy_osw_mci)$coefficients["insomnia", "Estimate"],
    summary(fit_cem_weighted_cemw_osw_unadj_svy_mci)$coefficients["insomnia", "Estimate"],
    summary(fit_cem_weighted_cemw_osw_adj_svy_mci)$coefficients["insomnia", "Estimate"]
  )
  # SE results
  se_mci[rep, ] <- c(
    summary(fit_psm_cov_unadj_svy_osw_mci)$coefficients["insomnia", "Std. Error"],
    summary(fit_psm_cov_adj_svy_osw_mci)$coefficients["insomnia", "Std. Error"],
    summary(fit_cem_weighted_cemw_osw_unadj_svy_mci)$coefficients["insomnia", "Std. Error"],
    summary(fit_cem_weighted_cemw_osw_adj_svy_mci)$coefficients["insomnia", "Std. Error"]
  )

}

# look at results
head(effect_mci)
head(se_mci)

# remove unneeded objects
rm(rep, list = ls(pattern = "fit"), analyses_mci)

### CONDUCT ANALYSES HYP -------------------------------------------------------------

# list the type of analyses we will conduct for hypertension 
analyses_hyp <- c("psm_cov_unadj_svy_osw_hyp", "psm_cov_adj_svy_osw_hyp",
                      "cem_weighted_cemw_osw_unadj_svy_hyp", "cem_weighted_cemw_osw_adj_svy_hyp")
# create matrix to hold effect estimates from the various analyses
effect_hyp <- matrix(NA, nrow = reps, ncol = length(analyses_hyp), dimnames = list(NULL, analyses_hyp))
# create matrix to hold SEs from the various analyses
se_hyp <- matrix(NA, nrow = reps, ncol = length(analyses_hyp), dimnames = list(NULL, analyses_hyp))

# loop through samples to conduct the various analyses on the matched data
for (rep in 1:reps) {
  if((rep %% 100) == 0) print(rep)
  
  # fit the various regressions 
  
  # using PSM data
  # when propensity score is calculated using logistic regression with weight as covariate
  # survey unadjusted regression using OSW
  fit_psm_cov_unadj_svy_osw_hyp <- svyglm(
    hyp_v2 ~ insomnia + offset(log(visit_years)), 
    family = poisson(), 
    design = list_matches_psm_cov_osw_svy[[rep]], 
    subset = hyp_v1 == 0
  )
  # survey adjusted regression using OSW
  fit_psm_cov_adj_svy_osw_hyp <- svyglm(
    hyp_v2 ~ insomnia + bmi + offset(log(visit_years)), 
    family = poisson(), 
    design = list_matches_psm_cov_osw_svy[[rep]],
    subset = hyp_v1 == 0
  )
  
  # using CEM data
  # when bins are defined by coarsened covariates + coarsened OSW
  # unadjusted survey regression with CEMW x OSW
  fit_cem_weighted_cemw_osw_unadj_svy_hyp <- svyglm(
    hyp_v2 ~ insomnia + offset(log(visit_years)), 
    family = poisson(),  
    design = list_matches_cem_weighted_cemw_osw_svy[[rep]],
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
  effect_hyp[rep, ] <- c(
    summary(fit_psm_cov_unadj_svy_osw_hyp)$coefficients["insomnia", "Estimate"],
    summary(fit_psm_cov_adj_svy_osw_hyp)$coefficients["insomnia", "Estimate"],
    summary(fit_cem_weighted_cemw_osw_unadj_svy_hyp)$coefficients["insomnia", "Estimate"],
    summary(fit_cem_weighted_cemw_osw_adj_svy_hyp)$coefficients["insomnia", "Estimate"]
  )
  # SE results
  se_hyp[rep, ] <- c(
    summary(fit_psm_cov_unadj_svy_osw_hyp)$coefficients["insomnia", "Std. Error"],
    summary(fit_psm_cov_adj_svy_osw_hyp)$coefficients["insomnia", "Std. Error"],
    summary(fit_cem_weighted_cemw_osw_unadj_svy_hyp)$coefficients["insomnia", "Std. Error"],
    summary(fit_cem_weighted_cemw_osw_adj_svy_hyp)$coefficients["insomnia", "Std. Error"]
  )
  
}

# look at results
head(effect_hyp)
head(se_hyp)

# remove unneeded objects
rm(rep, list = ls(pattern = "fit"), analyses_hyp)

### SAVE RESULTS ---------------------------------------------------------------------

# create list to hold effect estimate and SE results for the two outcomes
# where effect is given on the OR scale
list_effect_or <- list(mci = exp(effect_mci), hyp = exp(effect_hyp))
list_se <- list(mci = se_mci, hyp = se_hyp)

# save effect estimate and SE results for the two outcomes
saveRDS(list_effect_or, file =  paste0("./data/simulation_results/confounding_sensitivity_no_age/analyses_effect_", df_num, ".RData"))
saveRDS(list_se, file =  paste0("./data/simulation_results/confounding_sensitivity_no_age/analyses_se_", df_num, ".RData"))


