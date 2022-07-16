
### SET UP ---------------------------------------------------------------------

# load packages
library(tidyverse)

### CALCULATE PROP SCORE -------------------------------------------------------

conduct_one_weighting <- function(df_samp) {
  # fit logistic model to obtain propensity scores
  # based on based on weighted logistic regression using OSW
  fit_ps_weighted <- glm(insomnia ~ bmi + visit_years, weights = final_weights, family = binomial(), data = df_samp)
  # based on logistic regression using OSW as covariate
  fit_ps_cov <- glm(insomnia ~ bmi + visit_years + final_weights, family = binomial(), data = df_samp)

  # obtain propensity scores weights (PSW) using both ways for calculation (IPTW and weighting by the odds)
  # and add to sample data frame and calculate the product of PSW x OSW for both
  # for weighted logistic regression using OSW
  df_samp <- df_samp %>%
    mutate(
      # IPTW
      ps_weighted = predict(fit_ps_weighted, type = "response"),
      psw_weighted_iptw_weights = (insomnia / ps_weighted) + ((1 - insomnia) / (1 - ps_weighted)), 
      psw_weighted_osw_iptw_weights = psw_weighted_iptw_weights * final_weights,
      # weighting by the odds
      psw_weighted_odds_weights = insomnia + ((1 - insomnia) * ps_weighted / (1 - ps_weighted)),
      psw_weighted_osw_odds_weights = psw_weighted_odds_weights * final_weights
    )
  # for logistic regression using OSW as covariate
  df_samp <- df_samp %>%
    mutate(
      # IPTW
      ps_cov = predict(fit_ps_cov, type = "response"),
      psw_cov_iptw_weights = (insomnia / ps_cov) + ((1 - insomnia) / (1 - ps_cov)), 
      psw_cov_osw_iptw_weights = psw_cov_iptw_weights * final_weights,
      # weighting by the odds
      psw_cov_odds_weights = insomnia + ((1 - insomnia) * ps_cov / (1 - ps_cov)),
      psw_cov_osw_odds_weights = psw_cov_odds_weights * final_weights
    )
  
  return(df_samp)
}

