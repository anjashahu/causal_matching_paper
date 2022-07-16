
### SET UP ---------------------------------------------------------------------

# load packages
library(tidyverse)
library(MatchIt)

# note that cem package and matchit package should give us the same results 
# as we do not encounter any of the issues that can create differences in the results between the packages
# so we will use the matchit package since it's easier to work with

# note: we have removed age from this analysis

### CONDUCT CEM  ---------------------------------------------------------

# create function to conduct CEM 
conduct_cem <- function(df_samp) {
  # when CEM is based on binning of coarsened covariates + coarsened OSW
  df_matches_cem_weighted <- matchit(
    insomnia ~ bmi_cat + visit_years_cat + final_weights_cat,
    method = "cem", # CEM matching
    estimand = "ATT", 
    data = df_samp
  ) %>% 
    # obtain CEMW
    match.data(weights = "cem_weights") %>%
    # obtain CEMW x OSW
    mutate(cem_osw_weights = final_weights * cem_weights)
  
  # combine results
  return(df_matches_cem_weighted)
}


