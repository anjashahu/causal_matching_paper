
### SET UP ---------------------------------------------------------------------

# load packages
library(tidyverse)
library(MatchIt)

### GET INHERITED WEIGHTS -------------------------------------------------------
# create function to obtain inherited wights for propensity score matching method and add to data frame
get_inherited_weights <- function(df) {
  df %>%
    left_join(
      df %>% 
        filter(insomnia == 1) %>% 
        select(subclass, final_weights) %>%
        rename(inherited_weights = final_weights),
      by = "subclass"
    )
}

### CONDUCT MATCHING  ------------------------------------------------------

# create function to obtain matches (for psm)
conduct_one_psm <- function(df_samp) {
  # based on based on logistic regression with OSW as covariate
  # obtain matchit object  
  matchit_psm_cov <- matchit(
    insomnia ~ bmi + age + visit_years + final_weights,
    method = "nearest", # nearest neighbor matching 
    distance = "glm", # use propensity scores obtained from logistic regression
    estimand = "ATT", # controls are selected to be matched w/ treated 
    data = df_samp
  )
  # obtain match for the repetition and add to associated list
  df_matches_psm_cov <- matchit_psm_cov %>% 
    get_matches(
      weights = "matchit_weights", 
      id = "matchit_id"
    ) %>% 
    get_inherited_weights(.)
  
  return(df_matches_psm_cov)
}

