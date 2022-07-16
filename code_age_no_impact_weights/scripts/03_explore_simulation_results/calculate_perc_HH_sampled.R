
### SET UP ---------------------------------------------------------------------

# load packages
library(tidyverse)
library(survey)

# clear working environment
rm(list = ls())

# specify which of the individual data frames you want to use
specify_df <- "_01"
specify_df

# load data 
df_pop_bg <- read_csv("./data/simulated_pop_data/df_pop_bg.csv")
df_pop_hh <- read_csv("./data/simulated_pop_data/df_pop_hh.csv")
df_pop_indiv <- read_csv(paste0("./data/simulated_pop_data/df_pop_indiv", specify_df, ".csv")) # change this path to load the data frame of interest

# set numbers of repetitions
reps <- 1000
reps

### SIMULATE SAMPLES --------------------------------------------------------
# run R script to obtain function to simulate samples
source("./scripts/02_conduct_simulation_analyses/01_simulate_samples_function.R")

# create empty list to store samples
list_samps <- vector(mode = "list", length = reps)
# create empty list to store survey objects of samples that use OSW
list_samps_svy <- vector(mode = "list", length = reps)

# simulate all the samples and their survey design object by looping through the number of repetitions and storing each sample into list
for (rep in 1:reps) {
  if((rep %% 100) == 0) print(rep)
  # simulate sample for the repetition and add to associated list
  list_samps[[rep]] <- sim_one_sample(df_pop_bg, df_pop_hh, df_pop_indiv, rep) 
  # simulate survey design object for repetition and add to associated list
  list_samps_svy[[rep]] <- svydesign(
    ids = ~ bg + hh, 
    weights = ~ final_weights, 
    strata = ~ strata, 
    data = list_samps[[rep]]
  )
}

### CALCULATE PERCENTAGE OF HH SAMPLED IN COLLAPSED STRATA -----------------------

# create an empty vector to store the proportions of HH sampled in each collapsed strata
list_prop_HH <- vector(mode = "list", length = reps)

# loop through the samples to calculate the proportion of HH sampled in each collapsed strata
for (rep in 1:reps) {
  # create two variables in population data frame
  # one to indicate whether observation was sampled 
  # the second to indicate what collapsed strata the observation was in
  df_pop_update <- df_pop_indiv %>% 
    mutate(
      # variable of whether observation from population was sampled
      in_sample = df_pop_indiv$id %in% list_samps[[rep]]$id,
      # variable of collapsed strata
      collapse_strata = case_when(
        strata == 1 | strata == 2 ~ 1,
        strata == 3 | strata == 4 ~ 2,
        strata == 5 | strata == 6 ~ 3,
        strata == 7 | strata == 8 ~ 4
      )
    )
  
  # obtain proportion of sampled HH for each collapsed strata
  df_count_HH <- df_pop_update %>%
    # count the number of sampled and nonsampled HH for each collapsed strata
    count(collapse_strata, hh_hisp_surname, in_sample) %>%
    group_by(collapse_strata) %>%
    # create column with the total number of HHs (both sampled and nonsampled) for the collapsed strata
    # and use that column to create another column with the proportion of sampled HHs
    mutate(
      n_collapse_strata = sum(n),
      prop_in_sample = n / n_collapse_strata
    ) %>%
    # filter to just look at the counts and proportions for the sampled HHs
    filter(in_sample == TRUE) %>%
    # pull out the proportions for the sampled HHs
    pull(prop_in_sample)
  
  # add the proportions for this sample to the list containing the information for all samples
  list_prop_HH[[rep]] <- df_count_HH
}

# obtain the percentages of HH sampled in each collapsed strata for each sample
round(unlist(head(lapply(list_prop_HH, function(x) x * 100), 1)), digits = 2)



