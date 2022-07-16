
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

# we will vary intercept only to achieve prevalences of interest
# other coefficient values will come from the HCHS data

# EXPOSURE 0.05 -----
# data frame 1: exposure 0.05, outcome (v2) 0.05
#
# set insomnia coef values 
sim_or_vector <- c(hchs_or["intercept", "insomnia"], hchs_or["bmi", "insomnia"], hchs_or["age", "insomnia"])
sim_or_vector[1] <- 0.0113
# set outcome coef values
sim_or_df <- hchs_or[-4]
sim_or_df[1, 3] <- 0.0016
sim_or_df[1, 2] <- 0.000031
sim_or_df[1, 1] <- 0.000025
#
# set seed
set.seed(55)
# generate exposure
df_pop_indiv_1 <- generate_exposure(sim_or_vector, df_pop_indiv)
round(mean(df_pop_indiv_1$insomnia)*100, 1) # probability of having insomnia
# generate outcome
df_pop_indiv_1 <- generate_outcomes(sim_or_df, df_pop_indiv_1)
round(mean(df_pop_indiv_1$mci)*100, digits = 1) # probability of having mci under observed exposure status
round(mean(df_pop_indiv_1$hyp_v1)*100, digits = 1) # probability of having hypertension under observed exposure status for v1
round(mean(df_pop_indiv_1$hyp_v2)*100, digits = 1) # probability of having hypertension under observed exposure status for v2

# data frame 2: exposure 0.05, outcome (v2) 0.15
#
# set insomnia coef values 
sim_or_vector <- c(hchs_or["intercept", "insomnia"], hchs_or["bmi", "insomnia"], hchs_or["age", "insomnia"])
sim_or_vector[1] <- 0.0113
# set outcome coef values
sim_or_df <- hchs_or[-4]
sim_or_df[1, 3] <- 0.0067
sim_or_df[1, 2] <- 0.000222
sim_or_df[1, 1] <- 0.000122
#
# set seed
set.seed(55)
# generate exposure
df_pop_indiv_2 <- generate_exposure(sim_or_vector, df_pop_indiv)
round(mean(df_pop_indiv_2$insomnia)*100, 1) # probability of having insomnia
# generate outcome
df_pop_indiv_2 <- generate_outcomes(sim_or_df, df_pop_indiv_2)
round(mean(df_pop_indiv_2$mci)*100, digits = 1) # probability of having mci under observed exposure status
round(mean(df_pop_indiv_2$hyp_v1)*100, digits = 1) # probability of having hypertension under observed exposure status for v1
round(mean(df_pop_indiv_2$hyp_v2)*100, digits = 1) # probability of having hypertension under observed exposure status for v2

# data frame 3: exposure 0.05, outcome (v2) 0.25
#
# set insomnia coef values 
sim_or_vector <- c(hchs_or["intercept", "insomnia"], hchs_or["bmi", "insomnia"], hchs_or["age", "insomnia"])
sim_or_vector[1] <- 0.0113
# set outcome coef values
sim_or_df <- hchs_or[-4]
sim_or_df[1, 3] <- 0.015
sim_or_df[1, 2] <- 0.000587
sim_or_df[1, 1] <- 0.000302
#
# set seed
set.seed(55)
# generate exposure
df_pop_indiv_3 <- generate_exposure(sim_or_vector, df_pop_indiv)
round(mean(df_pop_indiv_3$insomnia)*100, 1) # probability of having insomnia
# generate outcome
df_pop_indiv_3 <- generate_outcomes(sim_or_df, df_pop_indiv_3)
round(mean(df_pop_indiv_3$mci)*100, digits = 1) # probability of having mci under observed exposure status
round(mean(df_pop_indiv_3$hyp_v1)*100, digits = 1) # probability of having hypertension under observed exposure status for v1
round(mean(df_pop_indiv_3$hyp_v2)*100, digits = 1) # probability of having hypertension under observed exposure status for v2

# data frame 4: exposure 0.05, outcome (v2) 0.35
#
# set insomnia coef values 
sim_or_vector <- c(hchs_or["intercept", "insomnia"], hchs_or["bmi", "insomnia"], hchs_or["age", "insomnia"])
sim_or_vector[1] <- 0.0113
# set outcome coef values
sim_or_df <- hchs_or[-4]
sim_or_df[1, 3] <- 0.0283
sim_or_df[1, 2] <- 0.001215
sim_or_df[1, 1] <- 0.000623
#
# set seed
set.seed(55)
# generate exposure
df_pop_indiv_4 <- generate_exposure(sim_or_vector, df_pop_indiv)
round(mean(df_pop_indiv_4$insomnia)*100, 1) # probability of having insomnia
# generate outcome
df_pop_indiv_4 <- generate_outcomes(sim_or_df, df_pop_indiv_4)
round(mean(df_pop_indiv_4$mci)*100, digits = 1) # probability of having mci under observed exposure status
round(mean(df_pop_indiv_4$hyp_v1)*100, digits = 1) # probability of having hypertension under observed exposure status for v1
round(mean(df_pop_indiv_4$hyp_v2)*100, digits = 1) # probability of having hypertension under observed exposure status for v2

# EXPOSURE 0.15 -----
# data frame 5: exposure 0.15, outcome (v2) 0.05
#
# set insomnia coef values 
sim_or_vector <- c(hchs_or["intercept", "insomnia"], hchs_or["bmi", "insomnia"], hchs_or["age", "insomnia"])
sim_or_vector[1] <- 0.0383
# set outcome coef values
sim_or_df <- hchs_or[-4]
sim_or_df[1, 3] <- 0.00153
sim_or_df[1, 2] <- 0.000031
sim_or_df[1, 1] <- 0.000024
#
# set seed
set.seed(55)
# generate exposure
df_pop_indiv_5 <- generate_exposure(sim_or_vector, df_pop_indiv)
round(mean(df_pop_indiv_5$insomnia)*100, 1) # probability of having insomnia
# generate outcome
df_pop_indiv_5 <- generate_outcomes(sim_or_df, df_pop_indiv_5)
round(mean(df_pop_indiv_5$mci)*100, digits = 1) # probability of having mci under observed exposure status
round(mean(df_pop_indiv_5$hyp_v1)*100, digits = 1) # probability of having hypertension under observed exposure status for v1
round(mean(df_pop_indiv_5$hyp_v2)*100, digits = 1) # probability of having hypertension under observed exposure status for v2

# data frame 6: exposure 0.15, outcome (v2) 0.15
#
# set insomnia coef values 
sim_or_vector <- c(hchs_or["intercept", "insomnia"], hchs_or["bmi", "insomnia"], hchs_or["age", "insomnia"])
sim_or_vector[1] <- 0.0383
# set outcome coef values
sim_or_df <- hchs_or[-4]
sim_or_df[1, 3] <- 0.00634
sim_or_df[1, 2] <- 0.000222
sim_or_df[1, 1] <- 0.000119
#
# set seed
set.seed(55)
# generate exposure
df_pop_indiv_6 <- generate_exposure(sim_or_vector, df_pop_indiv)
round(mean(df_pop_indiv_6$insomnia)*100, 1) # probability of having insomnia
# generate outcome
df_pop_indiv_6 <- generate_outcomes(sim_or_df, df_pop_indiv_6)
round(mean(df_pop_indiv_6$mci)*100, digits = 1) # probability of having mci under observed exposure status
round(mean(df_pop_indiv_6$hyp_v1)*100, digits = 1) # probability of having hypertension under observed exposure status for v1
round(mean(df_pop_indiv_6$hyp_v2)*100, digits = 1) # probability of having hypertension under observed exposure status for v2

# data frame 7: exposure 0.15, outcome (v2) 0.25
#
# set insomnia coef values 
sim_or_vector <- c(hchs_or["intercept", "insomnia"], hchs_or["bmi", "insomnia"], hchs_or["age", "insomnia"])
sim_or_vector[1] <- 0.0383
# set outcome coef values
sim_or_df <- hchs_or[-4]
sim_or_df[1, 3] <- 0.0143
sim_or_df[1, 2] <- 0.000582
sim_or_df[1, 1] <- 0.000295
#
# set seed
set.seed(55)
# generate exposure
df_pop_indiv_7 <- generate_exposure(sim_or_vector, df_pop_indiv)
round(mean(df_pop_indiv_7$insomnia)*100, 1) # probability of having insomnia
# generate outcome
df_pop_indiv_7 <- generate_outcomes(sim_or_df, df_pop_indiv_7)
round(mean(df_pop_indiv_7$mci)*100, digits = 1) # probability of having mci under observed exposure status
round(mean(df_pop_indiv_7$hyp_v1)*100, digits = 1) # probability of having hypertension under observed exposure status for v1
round(mean(df_pop_indiv_7$hyp_v2)*100, digits = 1) # probability of having hypertension under observed exposure status for v2

# data frame 8: exposure 0.15, outcome (v2) 0.35
#
# set insomnia coef values 
sim_or_vector <- c(hchs_or["intercept", "insomnia"], hchs_or["bmi", "insomnia"], hchs_or["age", "insomnia"])
sim_or_vector[1] <- 0.0383
# set outcome coef values
sim_or_df <- hchs_or[-4]
sim_or_df[1, 3] <- 0.027
sim_or_df[1, 2] <- 0.00121
sim_or_df[1, 1] <- 0.000609
#
# set seed
set.seed(55)
# generate exposure
df_pop_indiv_8 <- generate_exposure(sim_or_vector, df_pop_indiv)
round(mean(df_pop_indiv_8$insomnia)*100, 1) # probability of having insomnia
# generate outcome
df_pop_indiv_8 <- generate_outcomes(sim_or_df, df_pop_indiv_8)
round(mean(df_pop_indiv_8$mci)*100, digits = 1) # probability of having mci under observed exposure status
round(mean(df_pop_indiv_8$hyp_v1)*100, digits = 1) # probability of having hypertension under observed exposure status for v1
round(mean(df_pop_indiv_8$hyp_v2)*100, digits = 1) # probability of having hypertension under observed exposure status for v2


# EXPOSURE 0.25 -----
# data frame 9: exposure 0.25, outcome (v2) 0.05
#
# set insomnia coef values 
sim_or_vector <- c(hchs_or["intercept", "insomnia"], hchs_or["bmi", "insomnia"], hchs_or["age", "insomnia"])
sim_or_vector[1] <- 0.073
# set outcome coef values
sim_or_df <- hchs_or[-4]
sim_or_df[1, 3] <- 0.00145
sim_or_df[1, 2] <- 0.000031
sim_or_df[1, 1] <- 0.0000238
#
# set seed
set.seed(55)
# generate exposure
df_pop_indiv_9 <- generate_exposure(sim_or_vector, df_pop_indiv)
round(mean(df_pop_indiv_9$insomnia)*100, 1) # probability of having insomnia
# generate outcome
df_pop_indiv_9 <- generate_outcomes(sim_or_df, df_pop_indiv_9)
round(mean(df_pop_indiv_9$mci)*100, digits = 1) # probability of having mci under observed exposure status
round(mean(df_pop_indiv_9$hyp_v1)*100, digits = 1) # probability of having hypertension under observed exposure status for v1
round(mean(df_pop_indiv_9$hyp_v2)*100, digits = 1) # probability of having hypertension under observed exposure status for v2

# data frame 10: exposure 0.25, outcome (v2) 0.15
#
# set insomnia coef values 
sim_or_vector <- c(hchs_or["intercept", "insomnia"], hchs_or["bmi", "insomnia"], hchs_or["age", "insomnia"])
sim_or_vector[1] <- 0.073
# set outcome coef values
sim_or_df <- hchs_or[-4]
sim_or_df[1, 3] <- 0.00605
sim_or_df[1, 2] <- 0.00022
sim_or_df[1, 1] <- 0.000116
#
# set seed
set.seed(55)
# generate exposure
df_pop_indiv_10 <- generate_exposure(sim_or_vector, df_pop_indiv)
round(mean(df_pop_indiv_10$insomnia)*100, 1) # probability of having insomnia
# generate outcome
df_pop_indiv_10 <- generate_outcomes(sim_or_df, df_pop_indiv_10)
round(mean(df_pop_indiv_10$mci)*100, digits = 1) # probability of having mci under observed exposure status
round(mean(df_pop_indiv_10$hyp_v1)*100, digits = 1) # probability of having hypertension under observed exposure status for v1
round(mean(df_pop_indiv_10$hyp_v2)*100, digits = 1) # probability of having hypertension under observed exposure status for v2

# data frame 11: exposure 0.25, outcome (v2) 0.25
#
# set insomnia coef values 
sim_or_vector <- c(hchs_or["intercept", "insomnia"], hchs_or["bmi", "insomnia"], hchs_or["age", "insomnia"])
sim_or_vector[1] <- 0.073
# set outcome coef values
sim_or_df <- hchs_or[-4]
sim_or_df[1, 3] <- 0.0136
sim_or_df[1, 2] <- 0.00058
sim_or_df[1, 1] <- 0.000288
#
# set seed
set.seed(55)
# generate exposure
df_pop_indiv_11 <- generate_exposure(sim_or_vector, df_pop_indiv)
round(mean(df_pop_indiv_11$insomnia)*100, 1) # probability of having insomnia
# generate outcome
df_pop_indiv_11 <- generate_outcomes(sim_or_df, df_pop_indiv_11)
round(mean(df_pop_indiv_11$mci)*100, digits = 1) # probability of having mci under observed exposure status
round(mean(df_pop_indiv_11$hyp_v1)*100, digits = 1) # probability of having hypertension under observed exposure status for v1
round(mean(df_pop_indiv_11$hyp_v2)*100, digits = 1) # probability of having hypertension under observed exposure status for v2

# data frame 12: exposure 0.25, outcome (v2) 0.35
#
# set insomnia coef values 
sim_or_vector <- c(hchs_or["intercept", "insomnia"], hchs_or["bmi", "insomnia"], hchs_or["age", "insomnia"])
sim_or_vector[1] <- 0.073
# set outcome coef values
sim_or_df <- hchs_or[-4]
sim_or_df[1, 3] <- 0.0257
sim_or_df[1, 2] <- 0.0012
sim_or_df[1, 1] <- 0.000595
#
# set seed
set.seed(55)
# generate exposure
df_pop_indiv_12 <- generate_exposure(sim_or_vector, df_pop_indiv)
round(mean(df_pop_indiv_12$insomnia)*100, 1) # probability of having insomnia
# generate outcome
df_pop_indiv_12 <- generate_outcomes(sim_or_df, df_pop_indiv_12)
round(mean(df_pop_indiv_12$mci)*100, digits = 1) # probability of having mci under observed exposure status
round(mean(df_pop_indiv_12$hyp_v1)*100, digits = 1) # probability of having hypertension under observed exposure status for v1
round(mean(df_pop_indiv_12$hyp_v2)*100, digits = 1) # probability of having hypertension under observed exposure status for v2


# EXPOSURE 0.35 -----
# data frame 13: exposure 0.35, outcome (v2) 0.05
#
# set insomnia coef values 
sim_or_vector <- c(hchs_or["intercept", "insomnia"], hchs_or["bmi", "insomnia"], hchs_or["age", "insomnia"])
sim_or_vector[1] <- 0.12
# set outcome coef values
sim_or_df <- hchs_or[-4]
sim_or_df[1, 3] <- 0.00138
sim_or_df[1, 2] <- 0.000031
sim_or_df[1, 1] <- 0.000023
#
# set seed
set.seed(55)
# generate exposure
df_pop_indiv_13 <- generate_exposure(sim_or_vector, df_pop_indiv)
round(mean(df_pop_indiv_13$insomnia)*100, 1) # probability of having insomnia
# generate outcome
df_pop_indiv_13 <- generate_outcomes(sim_or_df, df_pop_indiv_13)
round(mean(df_pop_indiv_13$mci)*100, digits = 1) # probability of having mci under observed exposure status
round(mean(df_pop_indiv_13$hyp_v1)*100, digits = 1) # probability of having hypertension under observed exposure status for v1
round(mean(df_pop_indiv_13$hyp_v2)*100, digits = 1) # probability of having hypertension under observed exposure status for v2

# data frame 14: exposure 0.35, outcome (v2) 0.15
#
# set insomnia coef values 
sim_or_vector <- c(hchs_or["intercept", "insomnia"], hchs_or["bmi", "insomnia"], hchs_or["age", "insomnia"])
sim_or_vector[1] <- 0.12
# set outcome coef values
sim_or_df <- hchs_or[-4]
sim_or_df[1, 3] <- 0.0058
sim_or_df[1, 2] <- 0.000218
sim_or_df[1, 1] <- 0.000113
#
# set seed
set.seed(55)
# generate exposure
df_pop_indiv_14 <- generate_exposure(sim_or_vector, df_pop_indiv)
round(mean(df_pop_indiv_14$insomnia)*100, 1) # probability of having insomnia
# generate outcome
df_pop_indiv_14 <- generate_outcomes(sim_or_df, df_pop_indiv_14)
round(mean(df_pop_indiv_14$mci)*100, digits = 1) # probability of having mci under observed exposure status
round(mean(df_pop_indiv_14$hyp_v1)*100, digits = 1) # probability of having hypertension under observed exposure status for v1
round(mean(df_pop_indiv_14$hyp_v2)*100, digits = 1) # probability of having hypertension under observed exposure status for v2

# data frame 15: exposure 0.35, outcome (v2) 0.25
#
# set insomnia coef values 
sim_or_vector <- c(hchs_or["intercept", "insomnia"], hchs_or["bmi", "insomnia"], hchs_or["age", "insomnia"])
sim_or_vector[1] <- 0.12
# set outcome coef values
sim_or_df <- hchs_or[-4]
sim_or_df[1, 3] <- 0.013
sim_or_df[1, 2] <- 0.000575
sim_or_df[1, 1] <- 0.000282
#
# set seed
set.seed(55)
# generate exposure
df_pop_indiv_15 <- generate_exposure(sim_or_vector, df_pop_indiv)
round(mean(df_pop_indiv_15$insomnia)*100, 1) # probability of having insomnia
# generate outcome
df_pop_indiv_15 <- generate_outcomes(sim_or_df, df_pop_indiv_15)
round(mean(df_pop_indiv_15$mci)*100, digits = 1) # probability of having mci under observed exposure status
round(mean(df_pop_indiv_15$hyp_v1)*100, digits = 1) # probability of having hypertension under observed exposure status for v1
round(mean(df_pop_indiv_15$hyp_v2)*100, digits = 1) # probability of having hypertension under observed exposure status for v2

# data frame 16: exposure 0.35, outcome (v2) 0.35
#
# set insomnia coef values 
sim_or_vector <- c(hchs_or["intercept", "insomnia"], hchs_or["bmi", "insomnia"], hchs_or["age", "insomnia"])
sim_or_vector[1] <- 0.12
# set outcome coef values
sim_or_df <- hchs_or[-4]
sim_or_df[1, 3] <- 0.0245
sim_or_df[1, 2] <- 0.00119
sim_or_df[1, 1] <- 0.000582
#
# set seed
set.seed(55)
# generate exposure
df_pop_indiv_16 <- generate_exposure(sim_or_vector, df_pop_indiv)
round(mean(df_pop_indiv_16$insomnia)*100, 1) # probability of having insomnia
# generate outcome
df_pop_indiv_16 <- generate_outcomes(sim_or_df, df_pop_indiv_16)
round(mean(df_pop_indiv_16$mci)*100, digits = 1) # probability of having mci under observed exposure status
round(mean(df_pop_indiv_16$hyp_v1)*100, digits = 1) # probability of having hypertension under observed exposure status for v1
round(mean(df_pop_indiv_16$hyp_v2)*100, digits = 1) # probability of having hypertension under observed exposure status for v2

### SAVE ---------------------------------------------------------------------
# compile all the data frames
df_pop_indiv_list <- list(
  df_pop_indiv_1, df_pop_indiv_2, df_pop_indiv_3, df_pop_indiv_4, 
  df_pop_indiv_5, df_pop_indiv_6, df_pop_indiv_7, df_pop_indiv_8, 
  df_pop_indiv_9, df_pop_indiv_10, df_pop_indiv_11, df_pop_indiv_12,
  df_pop_indiv_13, df_pop_indiv_14, df_pop_indiv_15, df_pop_indiv_16
)
# save data frames separately 
for (df_num in 1:16) {
  path <- paste0("./data/simulated_pop_data/df_pop_indiv/df_pop_indiv_", df_num, ".csv")
  write_csv(df_pop_indiv_list[[df_num]], file = path)
}

