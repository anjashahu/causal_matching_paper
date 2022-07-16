
### SET UP ----------------------------------------------------------------------

# load packages
library(tidyverse)

### ESTIMATE TRUE POP EFFECTS ---------------------------------------------------

obtain_true_effect <- function(df_pop_indiv) {
  # for mci
  #
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
  #
  # marginal ATT
  pop_marg_att_mci <- coefficients(
    glm(mci_new ~ 
          insomnia_new, 
        family = binomial(), 
        data = df_pop_indiv_reformat_mci %>% 
          filter(insomnia == 1)) # filter dataset to those who actually had insomnia
  )["insomnia_new"]
  #
  # conditional ATT
  pop_cond_att_mci <- coefficients(
    glm(mci_new ~ 
          insomnia_new + bmi + age, 
        family = binomial(), 
        data = df_pop_indiv_reformat_mci %>% 
          filter(insomnia == 1)) # filter data set to those actually had  insomnia
  )["insomnia_new"]
  
  # for hypertension
  #
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
  #
  # marginal ATT
  pop_marg_att_hyp <- coefficients(
    glm(hyp_new_v2 ~ 
          insomnia_new + offset(log(visit_years)), 
        family = poisson(), 
        data = df_pop_indiv_reformat_hyp %>% 
          filter(hyp_new_v1 == 0 & insomnia == 1)) # filter to keep only hypertension free people at baseline and those who actually had insomnia
  )["insomnia_new"]
  #
  # conditional ATT
  pop_cond_att_hyp <- coefficients(
    glm(hyp_new_v2 ~ 
          insomnia_new + bmi + age + offset(log(visit_years)), 
        family = poisson(), 
        data = df_pop_indiv_reformat_hyp %>% 
          filter(hyp_new_v1 == 0 & insomnia == 1)) # filter to keep only hypertension free people at baseline and those with insomnia
  )["insomnia_new"]
  
  # create data frame to hold true pop coef estimates
  true_pop_coef <- data.frame(rbind(c(pop_marg_att_mci, pop_cond_att_mci), c(pop_marg_att_hyp, pop_cond_att_hyp)))
  colnames(true_pop_coef) <- c("marg_att", "cond_att")
  rownames(true_pop_coef) <- c("mci", "hyp")
  true_pop_coef
  # translate results to OR scale
  true_pop_or <- exp(true_pop_coef)
  return(true_pop_or)
}

### SAVE RESULTS ------------------------------------------------------------------

# loop through all data frames and save results separately on OR scale
for (df_num in 1:16) {
  # load individual data 
  df_pop_indiv <- read_csv(paste0("./data/simulated_pop_data/df_pop_indiv/df_pop_indiv_", df_num, ".csv")) 
  # get true effect estimates for the data frame
  true_pop_or <- obtain_true_effect(df_pop_indiv)
  print(true_pop_or)
  # create path to save results
  path <- paste0("./data/simulated_pop_data/true_pop_or/true_pop_or_", df_num, ".csv")
  # save results
  write_csv(true_pop_or, file = path)
}

