
### SET UP ----------------------------------------------------------------------

# load packages
library(tidyverse)

# clear working environment
rm(list = ls())

# specify the data frame you want to calculate true ORs for
specify_df <- "_01"
# load individual data 
df_pop_indiv <- read_csv(paste0("./data/simulated_pop_data/df_pop_indiv", specify_df, ".csv")) # change this path to load data frame of interest

### ESTIMATE TRUE POP EFFECTS ---------------------------------------------------

# estimate true ATT and ATE (both marginal and conditional) for each outcome
# for MCI
# reformat data set so that both mci counterfactuals of an individual have their own row 
df_pop_indiv_reformat_mci <- df_pop_indiv %>%
  pivot_longer(
    cols = mci_0:mci_1, 
    names_to = "insomnia_new", 
    values_to = "mci_new"
  ) %>%
  mutate(
    insomnia_new = ifelse(insomnia_new == "mci_0", 0, 1)
  )
# estimate true marginal effects
# marginal ATE
pop_marg_ate_mci <- coefficients(
  glm(mci_new ~ 
        insomnia_new, 
      family = binomial(), 
      data = df_pop_indiv_reformat_mci)
)["insomnia_new"]
pop_marg_ate_mci
# marginal ATT
pop_marg_att_mci <- coefficients(
  glm(mci_new ~ 
        insomnia_new, 
      family = binomial(), 
      data = df_pop_indiv_reformat_mci %>% 
        filter(insomnia == 1)) # filter dataset to those who actually had insomnia
)["insomnia_new"]
pop_marg_att_mci
# another way to obtain these values
#(mean(df_pop_indiv$mci_1)/(1-mean(df_pop_indiv$mci_1)))/(mean(df_pop_indiv$mci_0)/(1-mean(df_pop_indiv$mci_0))) # odds ratio
#log((mean(df_pop_indiv$mci_1)/(1-mean(df_pop_indiv$mci_1)))/(mean(df_pop_indiv$mci_0)/(1-mean(df_pop_indiv$mci_0)))) # log odds ratio
#(mean(df_pop_indiv$mci_1[df_pop_indiv$insomnia == 1])/(1-mean(df_pop_indiv$mci_1[df_pop_indiv$insomnia == 1])))/(mean(df_pop_indiv$mci_0[df_pop_indiv$insomnia == 1])/(1-mean(df_pop_indiv$mci_0[df_pop_indiv$insomnia == 1]))) # odds ratio
#log((mean(df_pop_indiv$mci_1[df_pop_indiv$insomnia == 1])/(1-mean(df_pop_indiv$mci_1[df_pop_indiv$insomnia == 1])))/(mean(df_pop_indiv$mci_0[df_pop_indiv$insomnia == 1])/(1-mean(df_pop_indiv$mci_0[df_pop_indiv$insomnia == 1])))) # log odds ratio

# estimate true conditional effects
# conditional ATE
pop_cond_ate_mci <- coefficients(
  glm(mci_new ~ 
        insomnia_new + bmi + age, 
      family = binomial(), 
      data = df_pop_indiv_reformat_mci)
)["insomnia_new"]
pop_cond_ate_mci
# conditional ATT
pop_cond_att_mci <- coefficients(
  glm(mci_new ~ 
        insomnia_new + bmi + age, 
      family = binomial(), 
      data = df_pop_indiv_reformat_mci %>% 
        filter(insomnia == 1)) # filter data set to those actually had  insomnia
)["insomnia_new"]
pop_cond_att_mci

### TRUE POP EFFECTS FOR HYP ---------------------------------------------------

# for poisson regression with log link there is no issue with noncollapsibility
# so the marginal effect estimates can be used as the conditional effect estimates

# reformat data set so that both hypertension counterfactuals of an individual have their own row 
df_pop_indiv_reformat_hyp <- df_pop_indiv %>%
  pivot_longer(
    cols = c("hyp_0_v1", "hyp_1_v1", "hyp_0_v2", "hyp_1_v2"), 
    names_to = c("insomnia_new", ".value"), 
    names_sep = "_v"
  ) %>%
  rename(
    "hyp_new_v1" = "1",
    "hyp_new_v2" = "2"
  ) %>%
  mutate(
    insomnia_new = ifelse(insomnia_new == "hyp_0", 0, 1)
  )

# estimate true effects
# ATE
pop_marg_ate_hyp <- coefficients(
  glm(hyp_new_v2 ~ 
        insomnia_new + offset(log(visit_years)), 
      family = poisson(), 
      data = df_pop_indiv_reformat_hyp %>% 
        filter(hyp_new_v1 == 0)) # filter to keep only hypertension free people at baseline
)["insomnia_new"]
exp(pop_marg_ate_hyp)
# can also get the answer by calculating
# sum(df_pop_indiv$hyp_1_v2[df_pop_indiv$hyp_1_v1 == 0])/sum(df_pop_indiv$visit_years[df_pop_indiv$hyp_1_v1 == 0]) / 
# (sum(df_pop_indiv$hyp_0_v2[df_pop_indiv$hyp_0_v1 == 0])/sum(df_pop_indiv$visit_years[df_pop_indiv$hyp_0_v1 == 0]))
# ATT
pop_marg_att_hyp <- coefficients(
  glm(hyp_new_v2 ~ 
        insomnia_new + offset(log(visit_years)), 
      family = poisson(), 
      data = df_pop_indiv_reformat_hyp %>% 
        filter(hyp_new_v1 == 0 & insomnia == 1)) # filter to keep only hypertension free people at baseline and those who actually had insomnia
)["insomnia_new"]
pop_marg_att_hyp
#sum(df_pop_indiv$hyp_1_v2[df_pop_indiv$hyp_1_v1 == 0 & df_pop_indiv$insomnia == 1])/sum(df_pop_indiv$visit_years[df_pop_indiv$hyp_1_v1 == 0 & df_pop_indiv$insomnia == 1]) / 
#(sum(df_pop_indiv$hyp_0_v2[df_pop_indiv$hyp_0_v1 == 0 & df_pop_indiv$insomnia == 1])/sum(df_pop_indiv$visit_years[df_pop_indiv$hyp_0_v1 == 0 & df_pop_indiv$insomnia == 1]))

# estimate true conditional effects
# conditional ATE
pop_cond_ate_hyp <- coefficients(
  glm(hyp_new_v2 ~ 
        insomnia_new + bmi + age + offset(log(visit_years)), 
      family = poisson(), 
      data = df_pop_indiv_reformat_hyp %>% 
        filter(hyp_new_v1 == 0)) # filter to keep only hypertension free people at baseline
)["insomnia_new"]
pop_cond_ate_hyp
# conditional ATT
pop_cond_att_hyp <- coefficients(
  glm(hyp_new_v2 ~ 
        insomnia_new + bmi + age + offset(log(visit_years)), 
      family = poisson(), 
      data = df_pop_indiv_reformat_hyp %>% 
        filter(hyp_new_v1 == 0 & insomnia == 1)) # filter to keep only hypertension free people at baseline and those with insomnia
)["insomnia_new"]
pop_cond_att_hyp

### SAVE RESULTS ------------------------------------------------------------------

# create data frame to hold true pop coef estimates
true_pop_coef <- data.frame(rbind(
  c(pop_marg_ate_mci, pop_marg_att_mci, pop_cond_ate_mci, pop_cond_att_mci),
  c(pop_marg_ate_hyp, pop_marg_att_hyp, pop_cond_ate_hyp, pop_cond_att_hyp)
))
colnames(true_pop_coef) <- c("marg_ate", "marg_att", "cond_ate", "cond_att")
rownames(true_pop_coef) <- c("mci", "hyp")
true_pop_coef
# translate results to OR scale
true_pop_or <- exp(true_pop_coef)
true_pop_or

# save results on OR scale
write_csv(true_pop_or, file = paste0("./data/simulated_pop_data/true_pop_or", specify_df, ".csv")) # change this path to specify data frame of interest


