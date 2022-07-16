
### SET UP ----------------------------------------------------------------------

# load packages
library(tidyverse)

# clear working environment
rm(list = ls())

# specify the data frame you want to calculate true ORs for
specify_df <- "_01"
# load individual data 
df_pop_indiv <- read_csv(paste0("./data/simulated_pop_data/df_pop_indiv", specify_df, ".csv")) # change this path to load data frame of interest

### LOOK AT BMI AND AGE EFFECTS -------------------------------------------------

# reformat data sets so that both mci counterfactuals of an individual have their own row 
df_pop_indiv_reformat_mci <- df_pop_indiv %>%
  pivot_longer(
    cols = mci_0:mci_1, 
    names_to = "insomnia_new", 
    values_to = "mci_new"
  ) %>%
  mutate(
    insomnia_new = ifelse(insomnia_new == "mci_0", 0, 1)
  )
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
# mci
exp(coef(glm(mci_new ~ 
      insomnia_new + bmi + age, 
    family = binomial(), 
    data = df_pop_indiv_reformat_mci)))

# hypertension
exp(coef(glm(hyp_new_v2 ~ 
      insomnia_new + bmi + age + offset(log(visit_years)), 
    family = poisson(), 
    data = df_pop_indiv_reformat_hyp %>% 
      filter(hyp_new_v1 == 0))))
