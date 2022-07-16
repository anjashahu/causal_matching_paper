
### SET UP ---------------------------------------------------------------------

# load packages
library(tidyverse)

# clear working environment
rm(list = ls())

# run R script to obtain the initialized population data frames
source("./scripts/01_simulate_pop/02_generate_pop_samp_design.R")

### GENERATE NON-EXPOSURE PREDICTORS -------------------------------------------

# set seed 
set.seed(100)

# generate remaining covariates (BMI and years between visits) and add to individual data frame
df_pop_indiv <- df_pop_indiv %>%
  mutate(
    # generate BMI from N(29, 81) distribution truncated to range 15 to 63
    bmi = rnorm(nrow(.), mean = 29, sd = 9),
    bmi = case_when(
      bmi < 15 ~ 15,
      bmi > 63 ~ 63,
      TRUE ~ bmi
    ),
    # generate years between visits from N(6, 0.25) distribution truncated to range 3 to 9
    visit_years = rnorm(nrow(.), mean = 6, sd = 0.5),
    visit_years = case_when(
      visit_years < 3 ~ 3, 
      visit_years > 9 ~ 9,
      TRUE ~ visit_years
    )
  )

### CALCULATE SAMPLING WEIGHTS --------------------------------------------------

# calculate weights
df_pop_indiv <- df_pop_indiv %>%
  mutate(
    # calculate sampling probability as the product of sampling probability of each stage
    samp_prob = bg_sampling_prob * hh_sampling_prob,
    # calculate base weights as the inverse of the sampling probability
    base_weights = 1 / samp_prob,
    # calculate final weights so that sum of the weights equals the population size
    final_weights = base_weights / mean(base_weights)
  )

### CATEGORIZE CONT VARIABLES ---------------------------------------------------

# generate categorical versions of continuous covariates 
# set cut points to categorize continuous variables 
age_cuts <- unname(quantile(df_pop_indiv$age))
bmi_cuts <- c(min(df_pop_indiv$bmi), 18.5, 25, 30, max(df_pop_indiv$bmi))
final_weights_cuts <- unname(quantile(df_pop_indiv$final_weights))
visit_years_cuts <- unname(quantile(df_pop_indiv$visit_years))
# add categorical variables to individual population data frame
df_pop_indiv <- df_pop_indiv %>% mutate(
  age_cat = as.numeric(cut(
    age, 
    age_cuts, 
    include.lowest = TRUE,
    right = FALSE
  )), 
  bmi_cat = as.numeric(cut(
    bmi, 
    bmi_cuts, 
    include.lowest = TRUE, 
    right = FALSE
  )), 
  final_weights_cat = as.numeric(cut(
    final_weights, final_weights_cuts, 
    include.lowest = TRUE, 
    right = FALSE
  )),
  visit_years_cat =  as.numeric(cut(
    visit_years, visit_years_cuts, 
    include.lowest = TRUE, 
    right = FALSE
  ))
)

### GENERATE CLUSTERING EFFECTS --------------------------------------------------

# generate clustering effects at BG and HH levels for outcome generation
df_pop_indiv <- df_pop_indiv %>% 
  mutate(
    bg_ce = rnorm(nrow(df_pop_bg), mean = 0, sd = 0.5)[bg],
    hh_ce = rnorm(nrow(df_pop_hh), mean = 0, sd = 1)[hh]
  )

### SAVE DATA --------------------------------------------------------------------

# save BG data frame for population
write_csv(df_pop_bg, "./data/simulated_pop_data/df_pop_bg.csv")
# save HH data frame for population
write_csv(df_pop_hh, "./data/simulated_pop_data/df_pop_hh.csv")
# save individual data frame for population
write_csv(df_pop_indiv, "./data/simulated_pop_data/df_pop_indiv_no_exposure_outcome.csv")

