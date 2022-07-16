
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
reps <- 3
reps

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

### CALCULATE PSW -----------------------------------------------------------

# run R script to obtain function to calculate propensity scores and weights
source("./scripts/02_conduct_simulation_analyses/02.03_conduct_weighting_function.R")

# loop through samples to obtain propensity scores and weights for each data frame
for (rep in 1:reps){
  if((rep %% 100) == 0) print(rep)
  list_samps[[rep]] <- conduct_one_weighting(list_samps[[rep]])
}

# remove unneeded objects
rm(rep, conduct_one_weighting)


### CREATE SURVEY DESIGN OBJECT -----------------------------------------------

# for weighting

# create empty list to store survey objects of samples that use PSW

# where propensity score is calculated using weighted logistic regression 
# for IPTW
list_samps_psw_weighted_iptw_svy <- vector(mode = "list", length = reps)
# for weighting by the odds
list_samps_psw_weighted_odds_svy <- vector(mode = "list", length = reps)

# where propensity score is calculated using logistic regression with weight as covariate
# for IPTW
list_samps_psw_cov_iptw_svy <- vector(mode = "list", length = reps)
# for weighting by the odds
list_samps_psw_cov_odds_svy <- vector(mode = "list", length = reps)


# create empty list to store survey objects of samples that use PSW x OSW

# where propensity score is calculated using weighted logistic regression 
# for IPTW
list_samps_psw_weighted_osw_iptw_svy <- vector(mode = "list", length = reps)
# for weighting by the odds
list_samps_psw_weighted_osw_odds_svy <- vector(mode = "list", length = reps)

# where propensity score is calculated using logistic regression with weight as covariate
# for IPTW
list_samps_psw_cov_osw_iptw_svy <- vector(mode = "list", length = reps)
# for weighting by the odds
list_samps_psw_cov_osw_odds_svy <- vector(mode = "list", length = reps)


# loop through the samples, create survey design objects with different weights and add to associated list
for (rep in 1:reps) {
  if((rep %% 100) == 0) print(rep)
  
  # for weighting using PSW
  
  # where propensity score is calculated using weighted logistic regression 
  # for IPTW
  list_samps_psw_weighted_iptw_svy[[rep]] <- svydesign(
    ids = ~ bg + hh, 
    weights = ~ psw_weighted_iptw_weights, 
    strata = ~ strata, 
    data = list_samps[[rep]]
  )  
  # for weighting by the odds
  list_samps_psw_weighted_odds_svy[[rep]] <- svydesign(
    ids = ~ bg + hh, 
    weights = ~ psw_weighted_odds_weights, 
    strata = ~ strata, 
    data = list_samps[[rep]]
  )
  
  # where propensity score is calculated using logistic regression with weight as covariate
  # for IPTW
  list_samps_psw_cov_iptw_svy[[rep]] <- svydesign(
    ids = ~ bg + hh, 
    weights = ~ psw_cov_iptw_weights, 
    strata = ~ strata, 
    data = list_samps[[rep]]
  ) 
  # for weighting by the odds
  list_samps_psw_cov_odds_svy[[rep]] <- svydesign(
    ids = ~ bg + hh, 
    weights = ~ psw_cov_odds_weights, 
    strata = ~ strata, 
    data = list_samps[[rep]]
  ) 
  
  
  # for weighting using PSW x OSW
  
  # where propensity score is calculated using weighted logistic regression 
  # for IPTW
  list_samps_psw_weighted_osw_iptw_svy[[rep]] <- svydesign(
    ids = ~ bg + hh, 
    weights = ~ psw_weighted_osw_iptw_weights, 
    strata = ~ strata, 
    data = list_samps[[rep]]
  ) 
  # for weighting by the odds
  list_samps_psw_weighted_osw_odds_svy[[rep]] <- svydesign(
    ids = ~ bg + hh, 
    weights = ~ psw_weighted_osw_odds_weights, 
    strata = ~ strata, 
    data = list_samps[[rep]]
  ) 
  
  # where propensity score is calculated using logistic regression with weight as covariate
  # for IPTW
  list_samps_psw_cov_osw_iptw_svy[[rep]] <- svydesign(
    ids = ~ bg + hh, 
    weights = ~ psw_cov_osw_iptw_weights, 
    strata = ~ strata, 
    data = list_samps[[rep]]
  ) 
  # for weighting by the odds
  list_samps_psw_cov_osw_odds_svy[[rep]] <- svydesign(
    ids = ~ bg + hh, 
    weights = ~ psw_cov_osw_odds_weights, 
    strata = ~ strata, 
    data = list_samps[[rep]]
  ) 
  
}

# remove unneeded variables
rm(rep)

### CONDUCT ANALYSES MCI ------------------------------------------------------

# list the type of analyses we will conduct for MCI for weighting
weighting_analyses_mci <- c(
  "psw_weighted_iptw_unadj_svy_mci", "psw_weighted_iptw_osw_unadj_svy_mci", 
  "psw_weighted_odds_unadj_svy_mci", "psw_weighted_odds_osw_unadj_svy_mci",
  "psw_weighted_iptw_adj_svy_mci", "psw_weighted_iptw_osw_adj_svy_mci", 
  "psw_weighted_odds_adj_svy_mci", "psw_weighted_odds_osw_adj_svy_mci",
  "psw_cov_iptw_unadj_svy_mci", "psw_cov_iptw_osw_unadj_svy_mci", 
  "psw_cov_odds_unadj_svy_mci", "psw_cov_odds_osw_unadj_svy_mci",
  "psw_cov_iptw_adj_svy_mci", "psw_cov_iptw_osw_adj_svy_mci", 
  "psw_cov_odds_adj_svy_mci", "psw_cov_odds_osw_adj_svy_mci"
)
# create matrix to hold effect estimates from the various analyses
weighting_effect_mci <- matrix(NA, nrow = reps, ncol = length(weighting_analyses_mci), dimnames = list(NULL, weighting_analyses_mci))
# create matrix to hold SEs from the various analyses
weighting_se_mci <- matrix(NA, nrow = reps, ncol = length(weighting_analyses_mci), dimnames = list(NULL, weighting_analyses_mci))

# loop through samples to conduct the various analyses on the full sample
for (rep in 1:reps) {
  if((rep %% 100) == 0) print(rep)
  
  # fit the various regressions
  
  # when propensity score is calculated using weighted logistic regression 
  # unadjusted survey regression with PSW using IPTW
  fit_psw_weighted_iptw_unadj_svy_mci <- svyglm(
    mci ~ insomnia, 
    family = quasibinomial(), 
    design = list_samps_psw_weighted_iptw_svy[[rep]]
  )
  # unadjusted survey regression with PSW x OSW using IPTW
  fit_psw_weighted_iptw_osw_unadj_svy_mci <- svyglm(
    mci ~ insomnia, 
    family = quasibinomial(), 
    design = list_samps_psw_weighted_osw_iptw_svy[[rep]]
  )
  # unadjusted survey regression with PSW using weighting by the odds
  fit_psw_weighted_odds_unadj_svy_mci <-  svyglm(
    mci ~ insomnia, 
    family = quasibinomial(), 
    design = list_samps_psw_weighted_odds_svy[[rep]]
  )
  # unadjusted survey regression with PSW x OSW using weighting by the odds
  fit_psw_weighted_odds_osw_unadj_svy_mci <-  svyglm(
    mci ~ insomnia, 
    family = quasibinomial(), 
    design = list_samps_psw_weighted_osw_odds_svy[[rep]]
  )
  # adjusted survey regression with PSW using IPTW
  fit_psw_weighted_iptw_adj_svy_mci <- svyglm(
    mci ~ insomnia + bmi + age, 
    family = quasibinomial(), 
    design = list_samps_psw_weighted_iptw_svy[[rep]]
  )
  
  # adjusted survey regression with PSW x OSW using IPTW
  fit_psw_weighted_iptw_osw_adj_svy_mci <- svyglm(
    mci ~ insomnia + bmi + age, 
    family = quasibinomial(), 
    design = list_samps_psw_weighted_osw_iptw_svy[[rep]]
  )
  # adjusted survey regression with PSW using weighting by the odds
  fit_psw_weighted_odds_adj_svy_mci <-  svyglm(
    mci ~ insomnia + bmi + age, 
    family = quasibinomial(), 
    design = list_samps_psw_weighted_odds_svy[[rep]]
  )
  # adjusted survey regression with PSW x OSW using weighting by the odds
  fit_psw_weighted_odds_osw_adj_svy_mci <-  svyglm(
    mci ~ insomnia + bmi + age, 
    family = quasibinomial(), 
    design = list_samps_psw_weighted_osw_odds_svy[[rep]]
  )

  # when propensity score is calculated using logistic regression with weight as covariate
  # unadjusted survey regression with PSW using IPTW
  fit_psw_cov_iptw_unadj_svy_mci <- svyglm(
    mci ~ insomnia, 
    family = quasibinomial(), 
    design = list_samps_psw_cov_iptw_svy[[rep]]
  )
  # unadjusted survey regression with PSW x OSW using IPTW
  fit_psw_cov_iptw_osw_unadj_svy_mci <- svyglm(
    mci ~ insomnia, 
    family = quasibinomial(), 
    design = list_samps_psw_cov_osw_iptw_svy[[rep]]
  )
  # unadjusted survey regression with PSW using weighting by the odds
  fit_psw_cov_odds_unadj_svy_mci <-  svyglm(
    mci ~ insomnia, 
    family = quasibinomial(), 
    design = list_samps_psw_cov_odds_svy[[rep]]
  )
  # unadjusted survey regression with PSW x OSW using weighting by the odds
  fit_psw_cov_odds_osw_unadj_svy_mci <-  svyglm(
    mci ~ insomnia, 
    family = quasibinomial(), 
    design = list_samps_psw_cov_osw_odds_svy[[rep]]
  )
  # adjusted survey regression with PSW using IPTW
  fit_psw_cov_iptw_adj_svy_mci <- svyglm(
    mci ~ insomnia + bmi + age, 
    family = quasibinomial(), 
    design = list_samps_psw_cov_iptw_svy[[rep]]
  )
  # adjusted survey regression with PSW x OSW using IPTW
  fit_psw_cov_iptw_osw_adj_svy_mci <- svyglm(
    mci ~ insomnia + bmi + age, 
    family = quasibinomial(), 
    design = list_samps_psw_cov_osw_iptw_svy[[rep]]
  )
  # adjusted survey regression with PSW using weighting by the odds
  fit_psw_cov_odds_adj_svy_mci <-  svyglm(
    mci ~ insomnia + bmi + age, 
    family = quasibinomial(), 
    design = list_samps_psw_cov_odds_svy[[rep]]
  )
  # adjusted survey regression with PSW x OSW using weighting by the odds
  fit_psw_cov_odds_osw_adj_svy_mci <-  svyglm(
    mci ~ insomnia + bmi + age, 
    family = quasibinomial(), 
    design = list_samps_psw_cov_osw_odds_svy[[rep]]
  )
  
  # add result to associated matrix
  # effect estimate results
  weighting_effect_mci[rep, ] <- c(
    summary(fit_psw_weighted_iptw_unadj_svy_mci)$coefficients["insomnia", "Estimate"],
    summary(fit_psw_weighted_iptw_osw_unadj_svy_mci)$coefficients["insomnia", "Estimate"],
    summary(fit_psw_weighted_odds_unadj_svy_mci)$coefficients["insomnia", "Estimate"],
    summary(fit_psw_weighted_odds_osw_unadj_svy_mci)$coefficients["insomnia", "Estimate"],
    summary(fit_psw_weighted_iptw_adj_svy_mci)$coefficients["insomnia", "Estimate"],
    summary(fit_psw_weighted_iptw_osw_adj_svy_mci)$coefficients["insomnia", "Estimate"],
    summary(fit_psw_weighted_odds_adj_svy_mci)$coefficients["insomnia", "Estimate"],
    summary(fit_psw_weighted_odds_osw_adj_svy_mci)$coefficients["insomnia", "Estimate"],
    
    summary(fit_psw_cov_iptw_unadj_svy_mci)$coefficients["insomnia", "Estimate"],
    summary(fit_psw_cov_iptw_osw_unadj_svy_mci)$coefficients["insomnia", "Estimate"],
    summary(fit_psw_cov_odds_unadj_svy_mci)$coefficients["insomnia", "Estimate"],
    summary(fit_psw_cov_odds_osw_unadj_svy_mci)$coefficients["insomnia", "Estimate"],
    summary(fit_psw_cov_iptw_adj_svy_mci)$coefficients["insomnia", "Estimate"],
    summary(fit_psw_cov_iptw_osw_adj_svy_mci)$coefficients["insomnia", "Estimate"],
    summary(fit_psw_cov_odds_adj_svy_mci)$coefficients["insomnia", "Estimate"],
    summary(fit_psw_cov_odds_osw_adj_svy_mci)$coefficients["insomnia", "Estimate"]
  )
  # SE results
  weighting_se_mci[rep, ] <- c(
    summary(fit_psw_weighted_iptw_unadj_svy_mci)$coefficients["insomnia", "Std. Error"],
    summary(fit_psw_weighted_iptw_osw_unadj_svy_mci)$coefficients["insomnia", "Std. Error"],
    summary(fit_psw_weighted_odds_unadj_svy_mci)$coefficients["insomnia", "Std. Error"],
    summary(fit_psw_weighted_odds_osw_unadj_svy_mci)$coefficients["insomnia", "Std. Error"],
    summary(fit_psw_weighted_iptw_adj_svy_mci)$coefficients["insomnia", "Std. Error"],
    summary(fit_psw_weighted_iptw_osw_adj_svy_mci)$coefficients["insomnia", "Std. Error"],
    summary(fit_psw_weighted_odds_adj_svy_mci)$coefficients["insomnia", "Std. Error"],
    summary(fit_psw_weighted_odds_osw_adj_svy_mci)$coefficients["insomnia", "Std. Error"],
    
    summary(fit_psw_cov_iptw_unadj_svy_mci)$coefficients["insomnia", "Std. Error"],
    summary(fit_psw_cov_iptw_osw_unadj_svy_mci)$coefficients["insomnia", "Std. Error"],
    summary(fit_psw_cov_odds_unadj_svy_mci)$coefficients["insomnia", "Std. Error"],
    summary(fit_psw_cov_odds_osw_unadj_svy_mci)$coefficients["insomnia", "Std. Error"],
    summary(fit_psw_cov_iptw_adj_svy_mci)$coefficients["insomnia", "Std. Error"],
    summary(fit_psw_cov_iptw_osw_adj_svy_mci)$coefficients["insomnia", "Std. Error"],
    summary(fit_psw_cov_odds_adj_svy_mci)$coefficients["insomnia", "Std. Error"],
    summary(fit_psw_cov_odds_osw_adj_svy_mci)$coefficients["insomnia", "Std. Error"]
  )
  
}

# look at results
head(weighting_effect_mci)
head(weighting_se_mci)

# remove unneeded objects
rm(rep, list = ls(pattern = "fit"), weighting_analyses_mci)

### CONDUCT ANALYSES HYP ------------------------------------------------------

# list the type of analyses we will conduct for hypertension for weighting
weighting_analyses_hyp <- c(
  "psw_weighted_iptw_unadj_svy_hyp", "psw_weighted_iptw_osw_unadj_svy_hyp", 
  "psw_weighted_odds_unadj_svy_hyp", "psw_weighted_odds_osw_unadj_svy_hyp",
  "psw_weighted_iptw_adj_svy_hyp", "psw_weighted_iptw_osw_adj_svy_hyp", 
  "psw_weighted_odds_adj_svy_hyp", "psw_weighted_odds_osw_adj_svy_hyp",
  "psw_cov_iptw_unadj_svy_hyp", "psw_cov_iptw_osw_unadj_svy_hyp", 
  "psw_cov_odds_unadj_svy_hyp", "psw_cov_odds_osw_unadj_svy_hyp",
  "psw_cov_iptw_adj_svy_hyp", "psw_cov_iptw_osw_adj_svy_hyp", 
  "psw_cov_odds_adj_svy_hyp", "psw_cov_odds_osw_adj_svy_hyp"
)
# create matrix to hold effect estimates from the various analyses
weighting_effect_hyp <- matrix(NA, nrow = reps, ncol = length(weighting_analyses_hyp), dimnames = list(NULL, weighting_analyses_hyp))
# create matrix to hold SEs from the various analyses
weighting_se_hyp <- matrix(NA, nrow = reps, ncol = length(weighting_analyses_hyp), dimnames = list(NULL, weighting_analyses_hyp))

# loop through samples to conduct the various analyses on the full sample
for (rep in 1:reps) {
  if((rep %% 100) == 0) print(rep)
  
  # fit the various regressions
  
  # when propensity score is calculated using weighted logistic regression 
  # unadjusted survey regression with PSW using IPTW
  fit_psw_weighted_iptw_unadj_svy_hyp <- svyglm(
    hyp_v2 ~ insomnia + offset(log(visit_years)), 
    family = poisson(), 
    design = list_samps_psw_weighted_iptw_svy[[rep]], 
    subset = hyp_v1 == 0
  )
  # unadjusted survey regression with PSW x OSW using IPTW
  fit_psw_weighted_iptw_osw_unadj_svy_hyp <- svyglm(
    hyp_v2 ~ insomnia + offset(log(visit_years)), 
    family = poisson(), 
    design = list_samps_psw_weighted_osw_iptw_svy[[rep]], 
    subset = hyp_v1 == 0
  )
  # unadjusted survey regression with PSW using weighting by the odds
  fit_psw_weighted_odds_unadj_svy_hyp <-  svyglm(
    hyp_v2 ~ insomnia + offset(log(visit_years)), 
    family = poisson(), 
    design = list_samps_psw_weighted_odds_svy[[rep]], 
    subset = hyp_v1 == 0
  )
  # unadjusted survey regression with PSW x OSW using weighting by the odds
  fit_psw_weighted_odds_osw_unadj_svy_hyp <-  svyglm(
    hyp_v2 ~ insomnia + offset(log(visit_years)), 
    family = poisson(), 
    design = list_samps_psw_weighted_osw_odds_svy[[rep]], 
    subset = hyp_v1 == 0
  )
  # adjusted survey regression with PSW using IPTW
  fit_psw_weighted_iptw_adj_svy_hyp <- svyglm(
    hyp_v2 ~ insomnia + bmi + age + offset(log(visit_years)), 
    family = poisson(), 
    design = list_samps_psw_weighted_iptw_svy[[rep]], 
    subset = hyp_v1 == 0
  )
  
  # adjusted survey regression with PSW x OSW using IPTW
  fit_psw_weighted_iptw_osw_adj_svy_hyp <- svyglm(
    hyp_v2 ~ insomnia + bmi + age + offset(log(visit_years)), 
    family = poisson(), 
    design = list_samps_psw_weighted_osw_iptw_svy[[rep]], 
    subset = hyp_v1 == 0
  )
  # adjusted survey regression with PSW using weighting by the odds
  fit_psw_weighted_odds_adj_svy_hyp <-  svyglm(
    hyp_v2 ~ insomnia + bmi + age + offset(log(visit_years)), 
    family = poisson(), 
    design = list_samps_psw_weighted_odds_svy[[rep]], 
    subset = hyp_v1 == 0
  )
  # adjusted survey regression with PSW x OSW using weighting by the odds
  fit_psw_weighted_odds_osw_adj_svy_hyp <-  svyglm(
    hyp_v2 ~ insomnia + bmi + age + offset(log(visit_years)), 
    family = poisson(), 
    design = list_samps_psw_weighted_osw_odds_svy[[rep]], 
    subset = hyp_v1 == 0
  )
  
  # when propensity score is calculated using logistic regression with weight as covariate
  # unadjusted survey regression with PSW using IPTW
  fit_psw_cov_iptw_unadj_svy_hyp <- svyglm(
    hyp_v2 ~ insomnia + offset(log(visit_years)), 
    family = poisson(), 
    design = list_samps_psw_cov_iptw_svy[[rep]], 
    subset = hyp_v1 == 0
  )
  # unadjusted survey regression with PSW x OSW using IPTW
  fit_psw_cov_iptw_osw_unadj_svy_hyp <- svyglm(
    hyp_v2 ~ insomnia + offset(log(visit_years)), 
    family = poisson(), 
    design = list_samps_psw_cov_osw_iptw_svy[[rep]], 
    subset = hyp_v1 == 0
  )
  # unadjusted survey regression with PSW using weighting by the odds
  fit_psw_cov_odds_unadj_svy_hyp <-  svyglm(
    hyp_v2 ~ insomnia + offset(log(visit_years)), 
    family = poisson(),  
    design = list_samps_psw_cov_odds_svy[[rep]], 
    subset = hyp_v1 == 0
  )
  # unadjusted survey regression with PSW x OSW using weighting by the odds
  fit_psw_cov_odds_osw_unadj_svy_hyp <-  svyglm(
    hyp_v2 ~ insomnia + offset(log(visit_years)), 
    family = poisson(), 
    design = list_samps_psw_cov_osw_odds_svy[[rep]], 
    subset = hyp_v1 == 0
  )
  # adjusted survey regression with PSW using IPTW
  fit_psw_cov_iptw_adj_svy_hyp <- svyglm(
    hyp_v2 ~ insomnia + bmi + age + offset(log(visit_years)), 
    family = poisson(), 
    design = list_samps_psw_cov_iptw_svy[[rep]], 
    subset = hyp_v1 == 0
  )
  # adjusted survey regression with PSW x OSW using IPTW
  fit_psw_cov_iptw_osw_adj_svy_hyp <- svyglm(
    hyp_v2 ~ insomnia + bmi + age + offset(log(visit_years)), 
    family = poisson(), 
    design = list_samps_psw_cov_osw_iptw_svy[[rep]], 
    subset = hyp_v1 == 0
  )
  # adjusted survey regression with PSW using weighting by the odds
  fit_psw_cov_odds_adj_svy_hyp <-  svyglm(
    hyp_v2 ~ insomnia + bmi + age + offset(log(visit_years)), 
    family = poisson(), 
    design = list_samps_psw_cov_odds_svy[[rep]], 
    subset = hyp_v1 == 0
  )
  # adjusted survey regression with PSW x OSW using weighting by the odds
  fit_psw_cov_odds_osw_adj_svy_hyp <-  svyglm(
    hyp_v2 ~ insomnia + bmi + age + offset(log(visit_years)), 
    family = poisson(), 
    design = list_samps_psw_cov_osw_odds_svy[[rep]], 
    subset = hyp_v1 == 0
  )
  
  # add result to associated matrix
  # effect estimate results
  weighting_effect_hyp[rep, ] <- c(
    summary(fit_psw_weighted_iptw_unadj_svy_hyp)$coefficients["insomnia", "Estimate"],
    summary(fit_psw_weighted_iptw_osw_unadj_svy_hyp)$coefficients["insomnia", "Estimate"],
    summary(fit_psw_weighted_odds_unadj_svy_hyp)$coefficients["insomnia", "Estimate"],
    summary(fit_psw_weighted_odds_osw_unadj_svy_hyp)$coefficients["insomnia", "Estimate"],
    summary(fit_psw_weighted_iptw_adj_svy_hyp)$coefficients["insomnia", "Estimate"],
    summary(fit_psw_weighted_iptw_osw_adj_svy_hyp)$coefficients["insomnia", "Estimate"],
    summary(fit_psw_weighted_odds_adj_svy_hyp)$coefficients["insomnia", "Estimate"],
    summary(fit_psw_weighted_odds_osw_adj_svy_hyp)$coefficients["insomnia", "Estimate"],
    
    summary(fit_psw_cov_iptw_unadj_svy_hyp)$coefficients["insomnia", "Estimate"],
    summary(fit_psw_cov_iptw_osw_unadj_svy_hyp)$coefficients["insomnia", "Estimate"],
    summary(fit_psw_cov_odds_unadj_svy_hyp)$coefficients["insomnia", "Estimate"],
    summary(fit_psw_cov_odds_osw_unadj_svy_hyp)$coefficients["insomnia", "Estimate"],
    summary(fit_psw_cov_iptw_adj_svy_hyp)$coefficients["insomnia", "Estimate"],
    summary(fit_psw_cov_iptw_osw_adj_svy_hyp)$coefficients["insomnia", "Estimate"],
    summary(fit_psw_cov_odds_adj_svy_hyp)$coefficients["insomnia", "Estimate"],
    summary(fit_psw_cov_odds_osw_adj_svy_hyp)$coefficients["insomnia", "Estimate"]
  )
  # SE results
  weighting_se_hyp[rep, ] <- c(
    summary(fit_psw_weighted_iptw_unadj_svy_hyp)$coefficients["insomnia", "Std. Error"],
    summary(fit_psw_weighted_iptw_osw_unadj_svy_hyp)$coefficients["insomnia", "Std. Error"],
    summary(fit_psw_weighted_odds_unadj_svy_hyp)$coefficients["insomnia", "Std. Error"],
    summary(fit_psw_weighted_odds_osw_unadj_svy_hyp)$coefficients["insomnia", "Std. Error"],
    summary(fit_psw_weighted_iptw_adj_svy_hyp)$coefficients["insomnia", "Std. Error"],
    summary(fit_psw_weighted_iptw_osw_adj_svy_hyp)$coefficients["insomnia", "Std. Error"],
    summary(fit_psw_weighted_odds_adj_svy_hyp)$coefficients["insomnia", "Std. Error"],
    summary(fit_psw_weighted_odds_osw_adj_svy_hyp)$coefficients["insomnia", "Std. Error"],
    
    summary(fit_psw_cov_iptw_unadj_svy_hyp)$coefficients["insomnia", "Std. Error"],
    summary(fit_psw_cov_iptw_osw_unadj_svy_hyp)$coefficients["insomnia", "Std. Error"],
    summary(fit_psw_cov_odds_unadj_svy_hyp)$coefficients["insomnia", "Std. Error"],
    summary(fit_psw_cov_odds_osw_unadj_svy_hyp)$coefficients["insomnia", "Std. Error"],
    summary(fit_psw_cov_iptw_adj_svy_hyp)$coefficients["insomnia", "Std. Error"],
    summary(fit_psw_cov_iptw_osw_adj_svy_hyp)$coefficients["insomnia", "Std. Error"],
    summary(fit_psw_cov_odds_adj_svy_hyp)$coefficients["insomnia", "Std. Error"],
    summary(fit_psw_cov_odds_osw_adj_svy_hyp)$coefficients["insomnia", "Std. Error"]
  )
  
}

# look at results
head(weighting_effect_hyp)
head(weighting_se_hyp)

# remove unneeded objects
rm(rep, list = ls(pattern = "fit"), weighting_analyses_hyp)

### CALCULATE BALANCE ----------------------------------------------------------------

# create function to calculate standardized mean difference (SMD) in a sample for weighting
calc_smd_weighting <- function(var, rep) {
  # means by exposure group using different weighting schemes for weighting
  means <- tibble(
    # mean weighted by PSW
    
    # where propensity score is calculated using weighted logistic regression 
    # for IPTW
    mean_psw_weighted_iptw = svyby(~ {{ var }}, by = ~ insomnia, list_samps_psw_weighted_iptw_svy[[rep]], svymean) %>% pull({{ var }}),
    # for weighting by the odds
    mean_psw_weighted_odds = svyby(~ {{ var }}, by = ~ insomnia, list_samps_psw_weighted_odds_svy[[rep]], svymean) %>% pull({{ var }}),
    
    # where propensity score is calculated using logistic regression with weight as covariate
    # for IPTW
    mean_psw_cov_iptw = svyby(~ {{ var }}, by = ~ insomnia, list_samps_psw_cov_iptw_svy[[rep]], svymean) %>% pull({{ var }}),
    # for weighting by the odds
    mean_psw_cov_odds = svyby(~ {{ var }}, by = ~ insomnia, list_samps_psw_cov_odds_svy[[rep]], svymean) %>% pull({{ var }}),
    
    # mean weighted by PSW x OSW
    
    # where propensity score is calculated using weighted logistic regression 
    # for IPTW
    mean_psw_weighted_osw_iptw = svyby(~ {{ var }}, by = ~ insomnia, list_samps_psw_weighted_osw_iptw_svy[[rep]], svymean) %>% pull({{ var }}),
    # for weighting by the odds
    mean_psw_weighted_osw_odds = svyby(~ {{ var }}, by = ~ insomnia, list_samps_psw_weighted_osw_odds_svy[[rep]], svymean) %>% pull({{ var }}),
    
    # where propensity score is calculated using logistic regression with weight as covariate
    # for IPTW
    mean_psw_cov_osw_iptw = svyby(~ {{ var }}, by = ~ insomnia, list_samps_psw_cov_osw_iptw_svy[[rep]], svymean) %>% pull({{ var }}),
    
    # for weighting by the odds
    mean_psw_cov_osw_odds = svyby(~ {{ var }}, by = ~ insomnia, list_samps_psw_cov_osw_odds_svy[[rep]], svymean) %>% pull({{ var }})
    
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
calc_smd_all_global_weighting <- function(rep) {
  # calculate smd for all three covariates for each method of using weights to calculate the means
  smd_all_global <- data.frame(
    age_smd = calc_smd_weighting(age, rep), 
    bmi_smd = calc_smd_weighting(bmi, rep), 
    visit_years_smd = calc_smd_weighting(visit_years, rep)
  ) %>% 
    # calculate the mean smd of the three variables 
    rowwise(.) %>%
    mutate(global_smd = (age_smd + bmi_smd + visit_years_smd) / 3) %>%
    # translate to data.frame object to allow for row names
    data.frame(row.names = names(calc_smd_weighting(age, rep)))
  
  # return result
  return(smd_all_global)
}

# create empty list to store smd calculations for each sample
list_balance_weighting <- vector(mode = "list", length = reps)

# loop through the samples to obtain the smd calculations for each sample
for (rep in 1:reps) {
  if((rep %% 100) == 0) print(rep)
  list_balance_weighting[[rep]] <- calc_smd_all_global_weighting(rep)
}

# find the average balance metrics across the samples for weighting
balance_avg_weighting <- Reduce(`+`, list_balance_weighting) / reps
balance_avg_weighting

# remove unneeded objects
rm(rep, calc_smd_weighting, calc_smd_all_global_weighting)

### SAVE RESULTS ---------------------------------------------------------------------

# create list to hold effect estimate and SE results for the two outcomes
# where effect is given on the OR scale
list_weighting_effect_or <- list(mci = exp(weighting_effect_mci), hyp = exp(weighting_effect_hyp))
list_weighting_se <- list(mci = weighting_se_mci, hyp = weighting_se_hyp)

# save effect estimate and SE results for the two outcomes
save(list_weighting_effect_or, list_weighting_se, file = paste0("./data/simulation_results/main_results/weighting_analyses", specify_df, ".RData")) # adjust path to specify which simulated data frame was used

# save balance metrics
save(balance_avg_weighting, file = paste0("./data/simulation_results/main_results/balance_avg_weighting", specify_df, ".RData")) # adjust path to specify which simulated data frame was used


