
### SET UP ---------------------------------------------------------------------

# load packages
library(tidyverse)
library(survey)

### SIMULATE ALL SAMPLES  ------------------------------------------------------

# create function to simulate a sample
sim_one_sample <- function(df_pop_bg, df_pop_hh, df_pop_indiv, seed) {
  
  set.seed(seed)
  
  # sample entire BGs and create BG data frame of sample
  df_samp_bg <- df_pop_bg %>%
    # sample BGs using fixed values (i.e. number of BGs sampled per strata follow 14, 14, 5, 5, 78, 78, 100, 100 for stratas 1 through 8, respectively)
    # (we sample based on almost fixed values instead of using rbinom() to ensure that a smaller strata like 3 or 4 
    # doesn't end up with very few BG sampled by chance which will not allow survey regressions to run)
    group_by(strata) %>%
    sample_frac(size = case_when(
      strata == 1 | strata == 2 | strata == 3 | strata == 4 ~ 0.25,
      strata == 5 | strata == 6 | strata == 7 | strata == 8 ~ 0.60
    )) %>% 
    ungroup()

  # from the sampled BGs, sample entire HHs (size = 2 for every HH) and create HH data frame of sample
  # note that this method of sampling HH means that some BGs that were "sampled" will not have any HHs from them sampled 
  # so those BGs will not end up in the final sample so the fixed number of BGs sampled is not really fixed (but it's still very close to fixed)
  df_samp_hh <- df_pop_hh %>%
    # filter HH data frame of population to HHs that are in sampled BGs
    filter(bg %in% df_samp_bg$bg) %>% 
    # sample HHs using HH sampling probabilities
    mutate(sampled = rbinom(nrow(.), size = 1, prob = hh_sampling_prob)) %>%
    # filter to HHs that were sampled
    filter(sampled == 1) %>%
    # remove column indicating sampling status
    select(!sampled)
  

   # create data frame of sampled individuals and calculate final weights
  df_samp_indiv <- df_pop_indiv %>% 
    # filter individual data frame of population to individuals that are in sampled HHs
    filter(hh %in% df_samp_hh$hh) %>%
    mutate(
      # recalculate final weights by dividing by the mean of the base weights so that sum of final weights adds up to sample size
      final_weights = base_weights / mean(base_weights)
    )
  
  # return individual data frame
  return(df_samp_indiv)
}





