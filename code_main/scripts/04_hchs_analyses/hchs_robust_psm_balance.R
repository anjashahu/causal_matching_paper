
### SET UP ---------------------------------------------------------------------

# load packages
library(tidyverse)
library(survey)
library(MatchIt)

# clear working environment
rm(list = ls())

# load data
hchs <- read_csv("./data/hchs_paper_data/solshahu_covariates_20210810.csv")
# remove unneeded variables and create insomnia variable
hchs <- hchs %>% 
  select(-c(METS_NCEP, SBD5, SLPDUR, DIABETES2_INDICATOR, DIABETES4_V2, METS_NCEP2_V2, CENTER)) %>%
  mutate(INSOMNIA = ifelse(WHIIRS >= 9, 1, 0), BKGRD1_C7 = BKGRD1_C7 + 1) %>%
  select(-WHIIRS) %>%
  rename_all(.funs = tolower) %>%
  mutate_at(
    .vars = c(
      "gender", "marital_status", "employed", 
      "alcohol_use", "cigarette_use", 
      "bkgrd1_c7", "education_c3"
    ), 
    .funs = as.factor
  ) %>%
  mutate(gender = ifelse(gender == "F", 0, 1)) %>%
  rename(background = bkgrd1_c7, education = education_c3, smoking = cigarette_use, employment = employed, alcohol = alcohol_use)

# function to get indicators for each level of a categorical variable
get_indicators <- function(df) {
  df %>% 
    # for marital status
    bind_cols(data.frame(model.matrix(~ marital_status - 1, data = df))) %>%
    # for employment
    bind_cols(data.frame(model.matrix(~ employment - 1, data = df))) %>%
    # for alcohol use
    bind_cols(data.frame(model.matrix(~ alcohol - 1, data = df))) %>%
    # for smoking 
    bind_cols(data.frame(model.matrix(~ smoking - 1, data = df))) %>%
    # for background
    bind_cols(data.frame(model.matrix(~ background - 1, data = df))) %>%
    # for education
    bind_cols(data.frame(model.matrix(~ education - 1, data = df)))
}

### CREATE MCI AND HYP DATAFRAMES -----------------------------------------------------------

# 1) for MCI
# filter further to just variables needed for mci analyses
# and remove missing values 
mci <- hchs %>%
  select(-c(weight_final_norm_overall, hypertension2_aha, hypertension2_aha_v2, mci3, yrs_btwn_v1v2)) %>%
  drop_na(insomnia, age, bmi, gender, marital_status, employment, alcohol, smoking, background, education, weight_norm_overall_inca) %>%
  get_indicators(.)
apply(mci, 2, function(col) sum(is.na(col)))

# 2) for hyp
hyp <- hchs %>%
  select(-c(weight_norm_overall_inca, mci, mci3)) %>%
  drop_na(insomnia, age, bmi, gender, marital_status, employment, alcohol, smoking, background, education, weight_final_norm_overall, yrs_btwn_v1v2) %>%
  get_indicators
apply(hyp, 2, function(col) sum(is.na(col)))

### CREATE FUNCTIONS TO CALC SMD ---------------------------------------------------------------

# function to calculate standardized mean difference (SMD) in a sample
calc_smd <- function(var, svy_object) {
  # mean weighted by OSW (before matching)
  mean_osw <- tibble(before_match = svyby(~ {{ var }}, by = ~ insomnia, svy_object, svymean) %>% pull({{ var }}))
  # sd in exposed group weighted by OSW (before matching)
  sd_osw <- tibble(sd_osw = svyby(~ {{ var }}, by = ~ insomnia, svy_object, svyvar) %>% pull({{ var }}) %>% sqrt(.)) %>% pull(sd_osw) %>% .[2]
  
  # calculate smd (before matching)
  # calculate (absolute) difference in means between exposure groups
  diff_mean <- abs(unlist(mean_osw[2, ] - mean_osw[1, ]))
  # calculate smd by dividing difference in means by before matching sd
  smd <- unname(diff_mean / sd_osw)
  return(smd)
}

# function to get a dataframe of SMD values for each variable
get_df_smd <- function(svy_object, analysis) {
  # calculate SMD for each variable and combine results
  df_smd <- cbind(
    # age
    age = calc_smd(age, svy_object),
    # bmi
    bmi = calc_smd(bmi, svy_object),
    # gender
    gender = calc_smd(gender, svy_object),
    # marital status 
    marital_status1 = calc_smd(marital_status1, svy_object),
    marital_status2 = calc_smd(marital_status2, svy_object),
    marital_status3 = calc_smd(marital_status3, svy_object),
    # employment 
    employment1 = calc_smd(employment1, svy_object),
    employment2 = calc_smd(employment2, svy_object),
    employment3 = calc_smd(employment3, svy_object),
    employment4 = calc_smd(employment4, svy_object),
    # alcohol use 
    alcohol1 = calc_smd(alcohol1, svy_object),
    alcohol2 = calc_smd(alcohol2, svy_object),
    alcohol3 = calc_smd(alcohol3, svy_object),
    # smoking
    smoking1 = calc_smd(smoking1, svy_object),
    smoking2 = calc_smd(smoking2, svy_object),
    smoking3 = calc_smd(smoking3, svy_object),
    # background
    background1 = calc_smd(background1, svy_object),
    background2 = calc_smd(background2, svy_object),
    background3 = calc_smd(background3, svy_object),
    background4 = calc_smd(background4, svy_object),
    background5 = calc_smd(background5, svy_object),
    background6 = calc_smd(background6, svy_object),
    background7 = calc_smd(background7, svy_object),
    # education
    education1 = calc_smd(education1, svy_object),
    education2 = calc_smd(education2, svy_object),
    education3 = calc_smd(education3, svy_object)
  ) 
  # translate to data frame with correct format 
  df_smd <- df_smd %>% 
    data.frame() %>%
    pivot_longer(cols = age:education3, names_to = "variable", values_to = analysis)
  return(df_smd)
}

### BEFORE MATCHING ---------------------------------------------------------------------

# 1) for MCI
# create survey design object
mci_svy <- svydesign(
  ids = ~ psu_id, 
  weights = ~ weight_norm_overall_inca, 
  strata = ~ strat, 
  data = mci
)
# calculate smd before matching
mci_smd_before_match <- get_df_smd(mci_svy, "before")


# 2) for hyp
# create survey design object
hyp_svy <- svydesign(
  ids = ~ psu_id, 
  weights = ~ weight_final_norm_overall, 
  strata = ~ strat, 
  data = hyp
)
# calculate smd before matching
hyp_smd_before_match <- get_df_smd(hyp_svy, "before")

### AFTER MATCHING ---------------------------------------------------------------------

# 1) for MCI
# get matched data frame after conducting robust PSM matching
# where you match based on logistic regression with OSW as covariate
mci_matchit_psm_cov <- matchit(
  insomnia ~ age + bmi + gender + marital_status + employment + alcohol + smoking + background + education + weight_norm_overall_inca,
  method = "nearest", # nearest neighbor matching 
  distance = "glm", # use propensity scores obtained from logistic regression
  estimand = "ATT", # controls are selected to be matched w/ treated 
  data = mci
)
mci_matched <- mci_matchit_psm_cov %>% 
  get_matches(
    weights = "matchit_weights", 
    id = "matchit_id"
  )
# create survey design objct
mci_matched_svy <- svydesign(
  ids = ~ psu_id, 
  weights = ~ weight_norm_overall_inca, 
  strata = ~ strat, 
  data = mci_matched
)
# calculate smd after matching
mci_smd_after_match <- get_df_smd(mci_matched_svy, "after")


# 2) for hyp
# get matched data frame after conducting robust PSM matching
# where you match based on logistic regression with OSW as covariate
hyp_matchit_psm_cov <- matchit(
  insomnia ~ age + bmi + gender + marital_status + employment + alcohol + smoking + background + education + yrs_btwn_v1v2 + weight_final_norm_overall,
  method = "nearest", # nearest neighbor matching 
  distance = "glm", # use propensity scores obtained from logistic regression
  estimand = "ATT", # controls are selected to be matched w/ treated 
  data = hyp
)
hyp_matched <- hyp_matchit_psm_cov %>% 
  get_matches(
    weights = "matchit_weights", 
    id = "matchit_id"
  )
# create survey design objct
hyp_matched_svy <-  svydesign(
  ids = ~ psu_id, 
  weights = ~ weight_final_norm_overall, 
  strata = ~ strat, 
  data = hyp_matched
)
# calculate smd after matching
hyp_smd_after_match <- get_df_smd(hyp_matched_svy, "after")


### COMBINE RESULTS AND PLOT ---------------------------------------------------------------------

# combine before and after results for each outcome
mci_smd <- mci_smd_before_match %>%
  left_join(mci_smd_after_match, by = "variable") %>%
  pivot_longer(before:after, names_to = "analysis", values_to = "smd") %>%
  mutate(outcome = rep("mci", nrow(.)))
hyp_smd <- hyp_smd_before_match %>%
  left_join(hyp_smd_after_match, by = "variable") %>%
  pivot_longer(before:after, names_to = "analysis", values_to = "smd") %>%
  mutate(outcome = rep("hypertension", nrow(.)))
# combine results from the two outcomes together
smd <- mci_smd %>%
  bind_rows(hyp_smd)

smd %>%
  ggplot(aes(x = smd, y = variable, color = analysis)) +
  geom_point() +
  facet_wrap(~outcome) +
  xlab("SMD") +
  ylab("Variable") +
  ggtitle("SMD before and after matching") +
  labs(color = "Analysis") +
  theme_bw()
ggsave("smd_plot", path = "./scripts/04_hchs_analyses/", device = "jpeg")


