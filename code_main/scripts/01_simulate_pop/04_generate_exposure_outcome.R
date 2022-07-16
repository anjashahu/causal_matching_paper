
### SET UP ---------------------------------------------------------------------

# load packages
library(tidyverse)

# clear working environment
rm(list = ls())

# create expit function to calculate probabilities 
expit <- function(x) {exp(x)/(1 + exp(x))}

# load the data frames that have non-exposure predictors 
df_pop_indiv <- read_csv("./data/simulated_pop_data/df_pop_indiv_no_exposure_outcome.csv")

# load the coef values (on OR scale) to use when generating exposure and outcome
load("./data/simulated_pop_data/hchs_or.RData")

### CREATE FUNCTIONS -----------------------------------------------------------

# create function to generate exposure
generate_exposure <- function(sim_or_vector, df_pop_indiv) { 
  # obtain coefficient values for generating exposure
  sim_coef_vector <- log(sim_or_vector)
  alpha0 <- sim_coef_vector[1]
  alpha1 <- sim_coef_vector[2]
  alpha2 <- sim_coef_vector[3]
  # generate insomnia and add to individual data frame
  df_pop_indiv <- df_pop_indiv %>%
    mutate(
      # generate probability of having insomnia that is dependent on age and bmi
      insomnia_p = expit(alpha0 + alpha1*bmi + alpha2*age),
      # generate binary (0/1) insomnia status
      insomnia = rbinom(nrow(.), size = 1, prob = insomnia_p)
    )
  return(df_pop_indiv)
}

# create function to generate counterfactual/potential outcomes and observed outcomes
generate_outcomes <- function(sim_or_df, df_pop_indiv) { # sim_coef is df holding coef values (as ORs)
  
  # calculate probability of outcome and translate to binary 0/1 variable using bernoulli(p)
  # need to generate both counterfactual outcomes and the observed outcomes for each of the two outcomes
  # for MCI simply generate for visit 2 only 
  # for hypertension generate for both visit 1 and visit 2
  
  # obtain coefficient values for generating outcomes
  sim_coef <- log(sim_or_df)
  # for mci
  beta0 <- sim_coef["intercept", "mci"]
  beta1 <- sim_coef["insomnia", "mci"]
  beta2 <- sim_coef["bmi", "mci"]
  beta3 <- sim_coef["age", "mci"]
  # for hyp
  # visit 1
  gamma0 <- sim_coef["intercept", "hyp_v1"]
  gamma1 <-  sim_coef["insomnia", "hyp_v1"]
  gamma2 <- sim_coef["bmi", "hyp_v1"]
  gamma3 <- sim_coef["age", "hyp_v1"]
  # visit 2
  phi0 <- sim_coef["intercept", "hyp_v2"]
  phi1 <- sim_coef["insomnia", "hyp_v2"]
  phi2 <- sim_coef["bmi", "hyp_v2"]
  phi3 <- sim_coef["age", "hyp_v2"]
  phi4 <- sim_coef["visit_years", "hyp_v2"]
  
  # now generate outcomes
  # MCI
  df_pop_indiv <- df_pop_indiv %>%
    mutate(
      # generate probability of mci under no insomnia
      mci_p0 = expit(beta0 + beta1*0 + beta2*bmi + beta3*age + bg_ce + hh_ce),
      # generate probability of mci under insomnia
      mci_p1 = expit(beta0 + beta1*1 + beta2*bmi + beta3*age + bg_ce + hh_ce),
      # translate to binary (0/1) mci status under no insomnia
      mci_0 = rbinom(n = nrow(.), size = 1, prob = mci_p0),
      # translate to binary (0/1) mci status under insomnia
      mci_1 = rbinom(n = nrow(.), size = 1, prob = mci_p1),
      # identify the outcome that actually occurred under the observed exposure
      mci = ifelse(insomnia == 0, mci_0, mci_1)
    )
  
  # hypertension
  df_pop_indiv <- df_pop_indiv %>%
    mutate(
      # for visit 1
      # generate probability of hypertension under no insomnia
      hyp_p0_v1 = expit(gamma0 + gamma1*0 + gamma2*bmi + gamma3*age + bg_ce + hh_ce),
      # generate probability of hypertension under insomnia
      hyp_p1_v1 = expit(gamma0 + gamma1*1 + gamma2*bmi + gamma3*age + bg_ce + hh_ce),
      # translate to binary (0/1) hypertension status under no insomnia
      hyp_0_v1 = rbinom(n = nrow(.), size = 1, prob = hyp_p0_v1),
      # translate to binary (0/1) hypertension status under insomnia
      hyp_1_v1 = rbinom(n = nrow(.), size = 1, prob = hyp_p1_v1),
      # identify the outcome that actually occurred under the observed exposure
      hyp_v1 = ifelse(insomnia == 0, hyp_0_v1, hyp_1_v1),
      
      # for visit 2
      # generate probability of hypertension under no insomnia
      hyp_p0_v2 = expit(phi0 + phi1*0 + phi2*bmi + phi3*age + phi4*visit_years + bg_ce + hh_ce),
      # generate probability of hypertension under insomnia
      hyp_p1_v2 = expit(phi0 + phi1*1 + phi2*bmi + phi3*age + phi4*visit_years + bg_ce + hh_ce),
      # translate to binary (0/1) hypertension status under no insomnia
      hyp_0_v2 = rbinom(n = nrow(.), size = 1, prob = hyp_p0_v2),
      # translate to binary (0/1) hypertension status under insomnia
      hyp_1_v2 = rbinom(n = nrow(.), size = 1, prob = hyp_p1_v2),
      # identify the outcome that actually occurred under the observed exposure
      hyp_v2 = ifelse(insomnia == 0, hyp_0_v2, hyp_1_v2)
    )
  
  return(df_pop_indiv)
}

### GENERATE DATA FRAMES -----------------------------------------------------------

# set insomnia coef values 
# based on hchs data 
sim_or_vector_hchs <- c(hchs_or["intercept", "insomnia"], hchs_or["bmi", "insomnia"], hchs_or["age", "insomnia"])
# based on other specification
sim_or_vector_other <- c(0.001, 1.07, 1.06)

# set outcome coef values
# based on hchs data
sim_or_df_hchs <- hchs_or[-4]
# based on other specification
sim_or_df_other <- sim_or_df_hchs
sim_or_df_other["insomnia", "hyp_v2"] <- 4
sim_or_df_other["bmi", "hyp_v2"] <- 1.11
sim_or_df_other["age", "hyp_v2"] <- 1.12
sim_or_df_other["insomnia", "mci"] <- 1.1
sim_or_df_other["age", "mci"] <- 1.03

# data frame 1: hchs-based exposure, hchs-based outcome
set.seed(55)

df_hchs_hchs <- generate_exposure(sim_or_vector_hchs, df_pop_indiv)
round(mean(df_hchs_hchs$insomnia_p)*100, 1)

df_hchs_hchs <- generate_outcomes(sim_or_df_hchs, df_hchs_hchs)
# round(mean(df_hchs_hchs$mci_0)*100, digits = 1) # probability of having mci under no insomnia
# round(mean(df_hchs_hchs$mci_1)*100, digits = 1) #  probability of having mci under insomnia
round(mean(df_hchs_hchs$mci)*100, digits = 1) # probability of having mci under observed exposure status
# round(mean(df_hchs_hchs$hyp_0_v1)*100, digits = 1) # probability of having hypertension under no insomnia for v1
# round(mean(df_hchs_hchs$hyp_1_v1)*100, digits = 1) # probability of having hypertension under insomnia for v1
round(mean(df_hchs_hchs$hyp_v1)*100, digits = 1) # probability of having hypertension under observed exposure status for v1
# round(mean(df_hchs_hchs$hyp_0_v2)*100, digits = 1) # probability of having hypertension under no insomnia for v2
# round(mean(df_hchs_hchs$hyp_1_v2)*100, digits = 1) # probability of having hypertension under insomnia for v2
round(mean(df_hchs_hchs$hyp_v2)*100, digits = 1) # probability of having hypertension under observed exposure status for v2

# data frame 2: hchs-based exposure, other outcome
set.seed(55)
df_hchs_other <- generate_exposure(sim_or_vector_hchs, df_pop_indiv)
round(mean(df_hchs_other$insomnia_p)*100, 1)

df_hchs_other <- generate_outcomes(sim_or_df_other, df_hchs_other)
# round(mean(df_hchs_other$mci_0)*100, digits = 1) # probability of having mci under no insomnia
# round(mean(df_hchs_other$mci_1)*100, digits = 1) #  probability of having mci under insomnia
round(mean(df_hchs_other$mci)*100, digits = 1) # probability of having mci under observed exposure status
# round(mean(df_hchs_other$hyp_0_v1)*100, digits = 1) # probability of having hypertension under no insomnia for v1
# round(mean(df_hchs_other$hyp_1_v1)*100, digits = 1) # probability of having hypertension under insomnia for v1
round(mean(df_hchs_other$hyp_v1)*100, digits = 1) # probability of having hypertension under observed exposure status for v1
# round(mean(df_hchs_other$hyp_0_v2)*100, digits = 1) # probability of having hypertension under no insomnia for v2
# round(mean(df_hchs_other$hyp_1_v2)*100, digits = 1) # probability of having hypertension under insomnia for v2
round(mean(df_hchs_other$hyp_v2)*100, digits = 1) # probability of having hypertension under observed exposure status for v2


# data frame 3: other exposure, hchs-based outcome
set.seed(55)
df_other_hchs <- generate_exposure(sim_or_vector_other, df_pop_indiv)
round(mean(df_other_hchs$insomnia_p)*100, 1)

df_other_hchs <- generate_outcomes(sim_or_df_hchs, df_other_hchs)
# round(mean(df_other_hchs$mci_0)*100, digits = 1) # probability of having mci under no insomnia
# round(mean(df_other_hchs$mci_1)*100, digits = 1) #  probability of having mci under insomnia
round(mean(df_other_hchs$mci)*100, digits = 1) # probability of having mci under observed exposure status
# round(mean(df_other_hchs$hyp_0_v1)*100, digits = 1) # probability of having hypertension under no insomnia for v1
# round(mean(df_other_hchs$hyp_1_v1)*100, digits = 1) # probability of having hypertension under insomnia for v1
round(mean(df_other_hchs$hyp_v1)*100, digits = 1) # probability of having hypertension under observed exposure status for v1
# round(mean(df_other_hchs$hyp_0_v2)*100, digits = 1) # probability of having hypertension under no insomnia for v2
# round(mean(df_other_hchs$hyp_1_v2)*100, digits = 1) # probability of having hypertension under insomnia for v2
round(mean(df_other_hchs$hyp_v2)*100, digits = 1) # probability of having hypertension under observed exposure status for v2


# data frame 4: other exposure, other outcome
set.seed(55)
df_other_other <- generate_exposure(sim_or_vector_other, df_pop_indiv)
round(mean(df_other_other$insomnia_p)*100, 1)

df_other_other <- generate_outcomes(sim_or_df_other, df_other_other)
# round(mean(df_other_other$mci_0)*100, digits = 1) # probability of having mci under no insomnia
# round(mean(df_other_other$mci_1)*100, digits = 1) #  probability of having mci under insomnia
round(mean(df_other_other$mci)*100, digits = 1) # probability of having mci under observed exposure status
# round(mean(df_other_other$hyp_0_v1)*100, digits = 1) # probability of having hypertension under no insomnia for v1
# round(mean(df_other_other$hyp_1_v1)*100, digits = 1) # probability of having hypertension under insomnia for v1
round(mean(df_other_other$hyp_v1)*100, digits = 1) # probability of having hypertension under observed exposure status for v1
# round(mean(df_other_other$hyp_0_v2)*100, digits = 1) # probability of having hypertension under no insomnia for v2
# round(mean(df_other_other$hyp_1_v2)*100, digits = 1) # probability of having hypertension under insomnia for v2
round(mean(df_other_other$hyp_v2)*100, digits = 1) # probability of having hypertension under observed exposure status for v2


### ALTER CONFOUNDERS STRENGTH -----------------------------------------------

sim_or_vector_confounder_new <- c(hchs_or["intercept", "insomnia"], 1.04, 1.01)
sim_or_df_confounder_new <- sim_or_df_hchs
sim_or_df_confounder_new["bmi", "hyp_v2"] <- 1.182
sim_or_df_confounder_new["age", "hyp_v2"] <- 1.03
sim_or_df_confounder_new["bmi", "hyp_v1"] <- 1.180
sim_or_df_confounder_new["age", "hyp_v1"] <- 1.02
sim_or_df_confounder_new["bmi", "mci"] <- 1.07
sim_or_df_confounder_new["age", "mci"] <- 1.02

set.seed(55)
df_confounder_new <- generate_exposure(sim_or_vector_confounder_new, df_pop_indiv)
round(mean(df_confounder_new$insomnia_p)*100, 1)

df_confounder_new <- generate_outcomes(sim_or_df_confounder_new, df_confounder_new)
# round(mean(df_confounder_new$mci_0)*100, digits = 1) # probability of having mci under no insomnia
# round(mean(df_confounder_new$mci_1)*100, digits = 1) #  probability of having mci under insomnia
round(mean(df_confounder_new$mci)*100, digits = 1) # probability of having mci under observed exposure status
# round(mean(df_confounder_new$hyp_0_v1)*100, digits = 1) # probability of having hypertension under no insomnia for v1
# round(mean(df_confounder_new$hyp_1_v1)*100, digits = 1) # probability of having hypertension under insomnia for v1
round(mean(df_confounder_new$hyp_v1)*100, digits = 1) # probability of having hypertension under observed exposure status for v1
# round(mean(df_confounder_new$hyp_0_v2)*100, digits = 1) # probability of having hypertension under no insomnia for v2
# round(mean(df_confounder_new$hyp_1_v2)*100, digits = 1) # probability of having hypertension under insomnia for v2
round(mean(df_confounder_new$hyp_v2)*100, digits = 1) # probability of having hypertension under observed exposure status for v2

write_csv(df_confounder_new, "./data/simulated_pop_data/df_pop_indiv_05.csv")

### SAVE ---------------------------------------------------------------------

write_csv(df_hchs_hchs, "./data/simulated_pop_data/df_pop_indiv_01.csv")
write_csv(df_hchs_other, "./data/simulated_pop_data/df_pop_indiv_02.csv")
write_csv(df_other_hchs, "./data/simulated_pop_data/df_pop_indiv_03.csv")
write_csv(df_other_other, "./data/simulated_pop_data/df_pop_indiv_04.csv")


