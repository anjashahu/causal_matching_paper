
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
  mutate(INSOMNIA = ifelse(WHIIRS >= 9, 1, 0)) %>%
  select(-WHIIRS) %>%
  rename_all(.funs = tolower) %>%
  mutate_at(
    .vars = c(
      "gender", "marital_status", "employed", 
      "alcohol_use", "cigarette_use", 
      "bkgrd1_c7", "education_c3"
    ), 
    .funs = as.factor
  )

# covariates:
# age (cont)
# bmi (cont)
# gender (categorical)
# marital_status (categorical)
# employed (categorical)
# alcohol_use (categorical)
# cigarette_use (categorical)
# bkgrd1_c7 (categorical)
# education_c3 (categorical)
# yrs_btwn_v1v2 (cont)

# filter further to just variables needed for mci analyses
# and remove missing values 
mci <- hchs %>%
  select(-c(weight_final_norm_overall, hypertension2_aha, hypertension2_aha_v2, mci3, yrs_btwn_v1v2)) %>%
  drop_na(insomnia, age, bmi, gender, marital_status, employed, alcohol_use, cigarette_use, bkgrd1_c7, education_c3, weight_norm_overall_inca)
apply(mci, 2, function(col) sum(is.na(col)))

### PSM ----------------------------------------------------------------------

# create function get inherited weights
get_inherited_weights <- function(df, weights) {
  df %>%
    left_join(
      df %>% 
        filter(insomnia == 1) %>% 
        select(subclass, {{weights}}) %>%
        rename(inherited_weights = {{weights}}),
      by = "subclass"
    )
}

# create formula to use in matching
formula_psm_weighted <- as.formula("insomnia ~ age + bmi + gender + marital_status + employed + alcohol_use + cigarette_use + bkgrd1_c7 + education_c3")
formula_psm_cov <- as.formula("insomnia ~ age + bmi + gender + marital_status + employed + alcohol_use + cigarette_use + bkgrd1_c7 + education_c3 + weight_norm_overall_inca")

# match based on weighted logistic regression using OSW
matchit_psm_weighted <- matchit(
  formula_psm_weighted,
  method = "nearest", # nearest neighbor matching 
  distance = "glm", # use propensity scores obtained from logistic regression
  estimand = "ATT", # controls are selected to be matched w/ treated 
  s.weights = ~ weight_norm_overall_inca, # use weighted regression
  data = mci
)
df_matches_psm_weighted <- matchit_psm_weighted %>% 
  get_matches(
    weights = "matchit_weights", 
    id = "matchit_id"
  ) %>% 
  get_inherited_weights(., "weight_norm_overall_inca")

# match based on logistic regression with OSW as covariate
matchit_psm_cov <- matchit(
  formula_psm_cov,
  method = "nearest", # nearest neighbor matching 
  distance = "glm", # use propensity scores obtained from logistic regression
  estimand = "ATT", # controls are selected to be matched w/ treated 
  data = mci
)
df_matches_psm_cov <- matchit_psm_cov %>% 
  get_matches(
    weights = "matchit_weights", 
    id = "matchit_id"
  ) %>% 
  get_inherited_weights(., "weight_norm_overall_inca")

### CREATE SURVEY DESIGN OBJECT ---------------------------------------------------------------------------------

# for PSM matches using OSW
# where propensity score is calculated using weighted logistic regression 
svy_matches_psm_weighted_osw <- svydesign(
  ids = ~ psu_id, 
  weights = ~ weight_norm_overall_inca, 
  strata = ~ strat, 
  data = df_matches_psm_weighted
)
# where propensity score is calculated using logistic regression with weight as covariate
svy_matches_psm_cov_osw <- svydesign(
  ids = ~ psu_id, 
  weights = ~ weight_norm_overall_inca, 
  strata = ~ strat, 
  data = df_matches_psm_cov
)
  
# for PSM matches that use ISW
# where propensity score is calculated using weighted logistic regression 
svy_matches_psm_weighted_isw <- svydesign(
  ids = ~ psu_id,  
  weights = ~ inherited_weights, 
  strata = ~ strat, 
  data = df_matches_psm_weighted
)
# where propensity score is calculated using logistic regression with weight as covariate
svy_matches_psm_cov_isw <- svydesign(
  ids = ~ psu_id,  
  weights = ~ inherited_weights, 
  strata = ~ strat, 
  data = df_matches_psm_cov
)

### CALCULATE EFFECTS ----------------------------------------------------------

options(survey.lonely.psu = "adjust")

# list the type of analyses we will conduct for MCI for psm
psm_analyses_mci <- c(
  "psm_weighted_unadj_mci", "psm_weighted_unadj_svy_osw_mci", "psm_weighted_unadj_svy_isw_mci",
  "psm_weighted_adj_mci", "psm_weighted_adj_svy_osw_mci", "psm_weighted_adj_svy_isw_mci",
  "psm_cov_unadj_mci", "psm_cov_unadj_svy_osw_mci", "psm_cov_unadj_svy_isw_mci",
  "psm_cov_adj_mci", "psm_cov_adj_svy_osw_mci", "psm_cov_adj_svy_isw_mci"
)
# create vector to hold effect estimates from the various analyses
psm_effect_mci <- rep(NA, length(psm_analyses_mci))
names(psm_effect_mci) <- psm_analyses_mci
# create vector to hold SEs from the various analyses
psm_se_mci <- rep(NA, length(psm_analyses_mci))
names(psm_se_mci) <- psm_analyses_mci
# matrix to hold CIs
psm_ci_mci <- matrix(NA, nrow = length(psm_analyses_mci), ncol = 2, dimnames = list(psm_analyses_mci, c("lower", "upper")))
# create vector to hold number of observations 
psm_nobs_mci <- rep(NA, length(psm_analyses_mci))
names(psm_nobs_mci) <- psm_analyses_mci

# create formulas for fitting regression models 
formula_unadj <- as.formula("mci ~ insomnia")
formula_adj <- as.formula("mci ~ insomnia + age + bmi + gender + marital_status + employed + alcohol_use + cigarette_use + bkgrd1_c7 + education_c3")

# when propensity score is calculated using weighted logistic regression 
# unadjusted regression without weights
fit_psm_weighted_unadj_mci <- glm(
  formula_unadj,
  family = binomial(), 
  data = df_matches_psm_weighted, 
  subset = !is.na(mci)
)
psm_effect_mci["psm_weighted_unadj_mci"] <- summary(fit_psm_weighted_unadj_mci)$coefficients["insomnia", "Estimate"]
psm_se_mci["psm_weighted_unadj_mci"] <- summary(fit_psm_weighted_unadj_mci)$coefficients["insomnia", "Std. Error"]
psm_ci_mci["psm_weighted_unadj_mci", ] <- confint(fit_psm_weighted_unadj_mci)["insomnia",]
psm_nobs_mci["psm_weighted_unadj_mci"] <- nobs(fit_psm_weighted_unadj_mci)

# survey unadjusted regression using OSW
fit_psm_weighted_unadj_svy_osw_mci <- svyglm(
  formula_unadj, 
  family = quasibinomial(), 
  design = svy_matches_psm_weighted_osw, 
  subset = !is.na(mci)
)
psm_effect_mci["psm_weighted_unadj_svy_osw_mci"] <- summary(fit_psm_weighted_unadj_svy_osw_mci)$coefficients["insomnia", "Estimate"]
psm_se_mci["psm_weighted_unadj_svy_osw_mci"] <- summary(fit_psm_weighted_unadj_svy_osw_mci)$coefficients["insomnia", "Std. Error"]
psm_ci_mci["psm_weighted_unadj_svy_osw_mci", ] <- confint(fit_psm_weighted_unadj_svy_osw_mci)["insomnia",]
psm_nobs_mci["psm_weighted_unadj_svy_osw_mci"] <- nobs(fit_psm_weighted_unadj_svy_osw_mci)

# survey unadjusted regression using ISW
fit_psm_weighted_unadj_svy_isw_mci <- svyglm(
  formula_unadj, 
  family = quasibinomial(), 
  design = svy_matches_psm_weighted_isw,
  subset = !is.na(mci)
)
psm_effect_mci["psm_weighted_unadj_svy_isw_mci"] <- summary(fit_psm_weighted_unadj_svy_isw_mci)$coefficients["insomnia", "Estimate"]
psm_se_mci["psm_weighted_unadj_svy_isw_mci"] <- summary(fit_psm_weighted_unadj_svy_isw_mci)$coefficients["insomnia", "Std. Error"]
psm_ci_mci["psm_weighted_unadj_svy_isw_mci", ] <- confint(fit_psm_weighted_unadj_svy_isw_mci)["insomnia",]
psm_nobs_mci["psm_weighted_unadj_svy_isw_mci"] <- nobs(fit_psm_weighted_unadj_svy_isw_mci)

# adjusted regression without weights
fit_psm_weighted_adj_mci <- glm(
  formula_adj, 
  family = binomial(), 
  data = df_matches_psm_weighted,
  subset = !is.na(mci)
)
psm_effect_mci["psm_weighted_adj_mci"] <- summary(fit_psm_weighted_adj_mci)$coefficients["insomnia", "Estimate"]
psm_se_mci["psm_weighted_adj_mci"] <- summary(fit_psm_weighted_adj_mci)$coefficients["insomnia", "Std. Error"]
psm_ci_mci["psm_weighted_adj_mci", ] <- confint(fit_psm_weighted_adj_mci)["insomnia",]
psm_nobs_mci["psm_weighted_adj_mci"] <- nobs(fit_psm_weighted_adj_mci)

# survey adjusted regression using OSW
fit_psm_weighted_adj_svy_osw_mci <- svyglm(
  formula_adj, 
  family = quasibinomial(), 
  design = svy_matches_psm_weighted_osw,
  subset = !is.na(mci)
)
psm_effect_mci["psm_weighted_adj_svy_osw_mci"] <- summary(fit_psm_weighted_adj_svy_osw_mci)$coefficients["insomnia", "Estimate"]
psm_se_mci["psm_weighted_adj_svy_osw_mci"] <- summary(fit_psm_weighted_adj_svy_osw_mci)$coefficients["insomnia", "Std. Error"]
psm_ci_mci["psm_weighted_adj_svy_osw_mci", ] <- confint(fit_psm_weighted_adj_svy_osw_mci)["insomnia",]
psm_nobs_mci["psm_weighted_adj_svy_osw_mci"] <- nobs(fit_psm_weighted_adj_svy_osw_mci)

# survey adjusted regression using ISW
fit_psm_weighted_adj_svy_isw_mci <- svyglm(
  formula_adj, 
  family = quasibinomial(), 
  design = svy_matches_psm_weighted_isw,
  subset = !is.na(mci)
)
psm_effect_mci["psm_weighted_adj_svy_isw_mci"] <- summary(fit_psm_weighted_adj_svy_isw_mci)$coefficients["insomnia", "Estimate"]
psm_se_mci["psm_weighted_adj_svy_isw_mci"] <- summary(fit_psm_weighted_adj_svy_isw_mci)$coefficients["insomnia", "Std. Error"]
psm_ci_mci["psm_weighted_adj_svy_isw_mci", ] <- confint(fit_psm_weighted_adj_svy_isw_mci)["insomnia",]
psm_nobs_mci["psm_weighted_adj_svy_isw_mci"] <- nobs(fit_psm_weighted_adj_svy_isw_mci)

# when propensity score is calculated using logistic regression with weight as covariate
# unadjusted regression without weights
fit_psm_cov_unadj_mci <- glm(
  formula_unadj, 
  family = binomial(), 
  data = df_matches_psm_cov,
  subset = !is.na(mci)
)
psm_effect_mci["psm_cov_unadj_mci"] <- summary(fit_psm_cov_unadj_mci)$coefficients["insomnia", "Estimate"]
psm_se_mci["psm_cov_unadj_mci"] <- summary(fit_psm_cov_unadj_mci)$coefficients["insomnia", "Std. Error"]
psm_ci_mci["psm_cov_unadj_mci", ] <- confint(fit_psm_cov_unadj_mci)["insomnia",]
psm_nobs_mci["psm_cov_unadj_mci"] <- nobs(fit_psm_cov_unadj_mci)

# survey unadjusted regression using OSW
fit_psm_cov_unadj_svy_osw_mci <- svyglm(
  formula_unadj, 
  family = quasibinomial(), 
  design = svy_matches_psm_cov_osw,
  subset = !is.na(mci)
)
psm_effect_mci["psm_cov_unadj_svy_osw_mci"] <- summary(fit_psm_cov_unadj_svy_osw_mci)$coefficients["insomnia", "Estimate"]
psm_se_mci["psm_cov_unadj_svy_osw_mci"] <- summary(fit_psm_cov_unadj_svy_osw_mci)$coefficients["insomnia", "Std. Error"]
psm_ci_mci["psm_cov_unadj_svy_osw_mci", ] <- confint(fit_psm_cov_unadj_svy_osw_mci)["insomnia",]
psm_nobs_mci["psm_cov_unadj_svy_osw_mci"] <- nobs(fit_psm_cov_unadj_svy_osw_mci)

# survey unadjusted regression using ISW
fit_psm_cov_unadj_svy_isw_mci <- svyglm(
  formula_unadj, 
  family = quasibinomial(), 
  design = svy_matches_psm_cov_isw,
  subset = !is.na(mci)
)
psm_effect_mci["psm_cov_unadj_svy_isw_mci"] <- summary(fit_psm_cov_unadj_svy_isw_mci)$coefficients["insomnia", "Estimate"]
psm_se_mci["psm_cov_unadj_svy_isw_mci"] <- summary(fit_psm_cov_unadj_svy_isw_mci)$coefficients["insomnia", "Std. Error"]
psm_ci_mci["psm_cov_unadj_svy_isw_mci", ] <- confint(fit_psm_cov_unadj_svy_isw_mci)["insomnia",]
psm_nobs_mci["psm_cov_unadj_svy_isw_mci"] <- nobs(fit_psm_cov_unadj_svy_isw_mci)

# adjusted regression without weights
fit_psm_cov_adj_mci <- glm(
  formula_adj, 
  family = binomial(), 
  data = df_matches_psm_cov,
  subset = !is.na(mci)
)
psm_effect_mci["psm_cov_adj_mci"] <- summary(fit_psm_cov_adj_mci)$coefficients["insomnia", "Estimate"]
psm_se_mci["psm_cov_adj_mci"] <- summary(fit_psm_cov_adj_mci)$coefficients["insomnia", "Std. Error"]
psm_ci_mci["psm_cov_adj_mci", ] <- confint(fit_psm_cov_adj_mci)["insomnia",]
psm_nobs_mci["psm_cov_adj_mci"] <- nobs(fit_psm_cov_adj_mci)

# survey adjusted regression using OSW
fit_psm_cov_adj_svy_osw_mci <- svyglm(
  formula_adj, 
  family = quasibinomial(), 
  design = svy_matches_psm_cov_osw,
  subset = !is.na(mci)
)
psm_effect_mci["psm_cov_adj_svy_osw_mci"] <- summary(fit_psm_cov_adj_svy_osw_mci)$coefficients["insomnia", "Estimate"]
psm_se_mci["psm_cov_adj_svy_osw_mci"] <- summary(fit_psm_cov_adj_svy_osw_mci)$coefficients["insomnia", "Std. Error"]
psm_ci_mci["psm_cov_adj_svy_osw_mci", ] <- confint(fit_psm_cov_adj_svy_osw_mci)["insomnia",]
psm_nobs_mci["psm_cov_adj_svy_osw_mci"] <- nobs(fit_psm_cov_adj_svy_osw_mci)

# survey adjusted regression using ISW
fit_psm_cov_adj_svy_isw_mci <- svyglm(
  formula_adj, 
  family = quasibinomial(), 
  design = svy_matches_psm_cov_isw,
  subset = !is.na(mci)
)
psm_effect_mci["psm_cov_adj_svy_isw_mci"] <- summary(fit_psm_cov_adj_svy_isw_mci)$coefficients["insomnia", "Estimate"]
psm_se_mci["psm_cov_adj_svy_isw_mci"] <- summary(fit_psm_cov_adj_svy_isw_mci)$coefficients["insomnia", "Std. Error"]
psm_ci_mci["psm_cov_adj_svy_isw_mci", ] <- confint(fit_psm_cov_adj_svy_isw_mci)["insomnia",]
psm_nobs_mci["psm_cov_adj_svy_isw_mci"] <- nobs(fit_psm_cov_adj_svy_isw_mci)






