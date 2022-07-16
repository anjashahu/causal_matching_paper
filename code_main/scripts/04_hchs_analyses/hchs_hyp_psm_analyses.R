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

# filter further to just variables needed for hyp analyses
# and remove missing values 
hyp <- hchs %>%
  select(-c(weight_norm_overall_inca, mci, mci3)) %>%
  drop_na(insomnia, age, bmi, gender, marital_status, employed, alcohol_use, cigarette_use, bkgrd1_c7, education_c3, weight_final_norm_overall, yrs_btwn_v1v2)
apply(hyp, 2, function(col) sum(is.na(col)))

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
formula_psm_weighted <- as.formula("insomnia ~ age + bmi + gender + marital_status + employed + alcohol_use + cigarette_use + bkgrd1_c7 + education_c3 + yrs_btwn_v1v2")
formula_psm_cov <- as.formula("insomnia ~ age + bmi + gender + marital_status + employed + alcohol_use + cigarette_use + bkgrd1_c7 + education_c3 + yrs_btwn_v1v2 + weight_final_norm_overall")

# match based on weighted logistic regression using OSW
matchit_psm_weighted <- matchit(
  formula_psm_weighted,
  method = "nearest", # nearest neighbor matching 
  distance = "glm", # use propensity scores obtained from logistic regression
  estimand = "ATT", # controls are selected to be matched w/ treated 
  s.weights = ~ weight_final_norm_overall, # use weighted regression
  data = hyp
)
df_matches_psm_weighted <- matchit_psm_weighted %>% 
  get_matches(
    weights = "matchit_weights", 
    id = "matchit_id"
  ) %>% 
  get_inherited_weights(., "weight_final_norm_overall")

# match based on logistic regression with OSW as covariate
matchit_psm_cov <- matchit(
  formula_psm_cov,
  method = "nearest", # nearest neighbor matching 
  distance = "glm", # use propensity scores obtained from logistic regression
  estimand = "ATT", # controls are selected to be matched w/ treated 
  data = hyp
)
df_matches_psm_cov <- matchit_psm_cov %>% 
  get_matches(
    weights = "matchit_weights", 
    id = "matchit_id"
  ) %>% 
  get_inherited_weights(., "weight_final_norm_overall")

### CREATE SURVEY DESIGN OBJECT ---------------------------------------------------------------------------------

# for PSM matches using OSW
# where propensity score is calculated using weighted logistic regression 
svy_matches_psm_weighted_osw <- svydesign(
  ids = ~ psu_id, 
  weights = ~ weight_final_norm_overall, 
  strata = ~ strat, 
  data = df_matches_psm_weighted
)
# where propensity score is calculated using logistic regression with weight as covariate
svy_matches_psm_cov_osw <- svydesign(
  ids = ~ psu_id, 
  weights = ~ weight_final_norm_overall, 
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
psm_analyses_hyp <- c(
  "psm_weighted_unadj_hyp", "psm_weighted_unadj_svy_osw_hyp", "psm_weighted_unadj_svy_isw_hyp",
  "psm_weighted_adj_hyp", "psm_weighted_adj_svy_osw_hyp", "psm_weighted_adj_svy_isw_hyp",
  "psm_cov_unadj_hyp", "psm_cov_unadj_svy_osw_hyp", "psm_cov_unadj_svy_isw_hyp",
  "psm_cov_adj_hyp", "psm_cov_adj_svy_osw_hyp", "psm_cov_adj_svy_isw_hyp"
)
# create vector to hold effect estimates from the various analyses
psm_effect_hyp <- rep(NA, length(psm_analyses_hyp))
names(psm_effect_hyp) <- psm_analyses_hyp
# create vector to hold SEs from the various analyses
psm_se_hyp <- rep(NA, length(psm_analyses_hyp))
names(psm_se_hyp) <- psm_analyses_hyp
# matrix to hold CIs
psm_ci_hyp <- matrix(NA, nrow = length(psm_analyses_hyp), ncol = 2, dimnames = list(psm_analyses_hyp, c("lower", "upper")))
# create vector to hold number of observations 
psm_nobs_hyp <- rep(NA, length(psm_analyses_hyp))
names(psm_nobs_hyp) <- psm_analyses_hyp

# create formulas for fitting regression models 
formula_unadj <- as.formula("hypertension2_aha_v2 ~ insomnia + offset(log(yrs_btwn_v1v2))")
formula_adj <- as.formula("hypertension2_aha_v2 ~ insomnia + age + bmi + gender + marital_status + employed + alcohol_use + cigarette_use + bkgrd1_c7 + education_c3 + offset(log(yrs_btwn_v1v2))")

# when propensity score is calculated using weighted logistic regression 
# unadjusted regression without weights
fit_psm_weighted_unadj_hyp <- glm(
  formula_unadj,
  family = poisson(), 
  data = df_matches_psm_weighted, 
  subset = !is.na(hypertension2_aha_v2) & hypertension2_aha == 0
)
psm_effect_hyp["psm_weighted_unadj_hyp"] <- summary(fit_psm_weighted_unadj_hyp)$coefficients["insomnia", "Estimate"]
psm_se_hyp["psm_weighted_unadj_hyp"] <- summary(fit_psm_weighted_unadj_hyp)$coefficients["insomnia", "Std. Error"]
psm_ci_hyp["psm_weighted_unadj_hyp", ] <- confint(fit_psm_weighted_unadj_hyp)["insomnia",]
psm_nobs_hyp["psm_weighted_unadj_hyp"] <- nobs(fit_psm_weighted_unadj_hyp)

# survey unadjusted regression using OSW
fit_psm_weighted_unadj_svy_osw_hyp <- svyglm(
  formula_unadj, 
  family = poisson(), 
  design = svy_matches_psm_weighted_osw, 
  subset = !is.na(hypertension2_aha_v2) & hypertension2_aha == 0
)
psm_effect_hyp["psm_weighted_unadj_svy_osw_hyp"] <- summary(fit_psm_weighted_unadj_svy_osw_hyp)$coefficients["insomnia", "Estimate"]
psm_se_hyp["psm_weighted_unadj_svy_osw_hyp"] <- summary(fit_psm_weighted_unadj_svy_osw_hyp)$coefficients["insomnia", "Std. Error"]
psm_ci_hyp["psm_weighted_unadj_svy_osw_hyp", ] <- confint(fit_psm_weighted_unadj_svy_osw_hyp)["insomnia",]
psm_nobs_hyp["psm_weighted_unadj_svy_osw_hyp"] <- nobs(fit_psm_weighted_unadj_svy_osw_hyp)

# survey unadjusted regression using ISW
fit_psm_weighted_unadj_svy_isw_hyp <- svyglm(
  formula_unadj, 
  family = poisson(), 
  design = svy_matches_psm_weighted_isw,
  subset = !is.na(hypertension2_aha_v2) & hypertension2_aha == 0
)
psm_effect_hyp["psm_weighted_unadj_svy_isw_hyp"] <- summary(fit_psm_weighted_unadj_svy_isw_hyp)$coefficients["insomnia", "Estimate"]
psm_se_hyp["psm_weighted_unadj_svy_isw_hyp"] <- summary(fit_psm_weighted_unadj_svy_isw_hyp)$coefficients["insomnia", "Std. Error"]
psm_ci_hyp["psm_weighted_unadj_svy_isw_hyp", ] <- confint(fit_psm_weighted_unadj_svy_isw_hyp)["insomnia",]
psm_nobs_hyp["psm_weighted_unadj_svy_isw_hyp"] <- nobs(fit_psm_weighted_unadj_svy_isw_hyp)

# adjusted regression without weights
fit_psm_weighted_adj_hyp <- glm(
  formula_adj, 
  family = poisson(), 
  data = df_matches_psm_weighted,
  subset = !is.na(hypertension2_aha_v2) & hypertension2_aha == 0
)
psm_effect_hyp["psm_weighted_adj_hyp"] <- summary(fit_psm_weighted_adj_hyp)$coefficients["insomnia", "Estimate"]
psm_se_hyp["psm_weighted_adj_hyp"] <- summary(fit_psm_weighted_adj_hyp)$coefficients["insomnia", "Std. Error"]
psm_ci_hyp["psm_weighted_adj_hyp", ] <- confint(fit_psm_weighted_adj_hyp)["insomnia",]
psm_nobs_hyp["psm_weighted_adj_hyp"] <- nobs(fit_psm_weighted_adj_hyp)

# survey adjusted regression using OSW
fit_psm_weighted_adj_svy_osw_hyp <- svyglm(
  formula_adj, 
  family = poisson(), 
  design = svy_matches_psm_weighted_osw,
  subset = !is.na(hypertension2_aha_v2) & hypertension2_aha == 0
)
psm_effect_hyp["psm_weighted_adj_svy_osw_hyp"] <- summary(fit_psm_weighted_adj_svy_osw_hyp)$coefficients["insomnia", "Estimate"]
psm_se_hyp["psm_weighted_adj_svy_osw_hyp"] <- summary(fit_psm_weighted_adj_svy_osw_hyp)$coefficients["insomnia", "Std. Error"]
psm_ci_hyp["psm_weighted_adj_svy_osw_hyp", ] <- confint(fit_psm_weighted_adj_svy_osw_hyp)["insomnia",]
psm_nobs_hyp["psm_weighted_adj_svy_osw_hyp"] <- nobs(fit_psm_weighted_adj_svy_osw_hyp)

# survey adjusted regression using ISW
fit_psm_weighted_adj_svy_isw_hyp <- svyglm(
  formula_adj, 
  family = poisson(), 
  design = svy_matches_psm_weighted_isw,
  subset = !is.na(hypertension2_aha_v2) & hypertension2_aha == 0
)
psm_effect_hyp["psm_weighted_adj_svy_isw_hyp"] <- summary(fit_psm_weighted_adj_svy_isw_hyp)$coefficients["insomnia", "Estimate"]
psm_se_hyp["psm_weighted_adj_svy_isw_hyp"] <- summary(fit_psm_weighted_adj_svy_isw_hyp)$coefficients["insomnia", "Std. Error"]
psm_ci_hyp["psm_weighted_adj_svy_isw_hyp", ] <- confint(fit_psm_weighted_adj_svy_isw_hyp)["insomnia",]
psm_nobs_hyp["psm_weighted_adj_svy_isw_hyp"] <- nobs(fit_psm_weighted_adj_svy_isw_hyp)

# when propensity score is calculated using logistic regression with weight as covariate
# unadjusted regression without weights
fit_psm_cov_unadj_hyp <- glm(
  formula_unadj, 
  family = poisson(), 
  data = df_matches_psm_cov,
  subset = !is.na(hypertension2_aha_v2) & hypertension2_aha == 0
)
psm_effect_hyp["psm_cov_unadj_hyp"] <- summary(fit_psm_cov_unadj_hyp)$coefficients["insomnia", "Estimate"]
psm_se_hyp["psm_cov_unadj_hyp"] <- summary(fit_psm_cov_unadj_hyp)$coefficients["insomnia", "Std. Error"]
psm_ci_hyp["psm_cov_unadj_hyp", ] <- confint(fit_psm_cov_unadj_hyp)["insomnia",]
psm_nobs_hyp["psm_cov_unadj_hyp"] <- nobs(fit_psm_cov_unadj_hyp)

# survey unadjusted regression using OSW
fit_psm_cov_unadj_svy_osw_hyp <- svyglm(
  formula_unadj, 
  family = poisson(), 
  design = svy_matches_psm_cov_osw,
  subset = !is.na(hypertension2_aha_v2) & hypertension2_aha == 0
)
psm_effect_hyp["psm_cov_unadj_svy_osw_hyp"] <- summary(fit_psm_cov_unadj_svy_osw_hyp)$coefficients["insomnia", "Estimate"]
psm_se_hyp["psm_cov_unadj_svy_osw_hyp"] <- summary(fit_psm_cov_unadj_svy_osw_hyp)$coefficients["insomnia", "Std. Error"]
psm_ci_hyp["psm_cov_unadj_svy_osw_hyp", ] <- confint(fit_psm_cov_unadj_svy_osw_hyp)["insomnia",]
psm_nobs_hyp["psm_cov_unadj_svy_osw_hyp"] <- nobs(fit_psm_cov_unadj_svy_osw_hyp)

# survey unadjusted regression using ISW
fit_psm_cov_unadj_svy_isw_hyp <- svyglm(
  formula_unadj, 
  family = poisson(), 
  design = svy_matches_psm_cov_isw,
  subset = !is.na(hypertension2_aha_v2) & hypertension2_aha == 0
)
psm_effect_hyp["psm_cov_unadj_svy_isw_hyp"] <- summary(fit_psm_cov_unadj_svy_isw_hyp)$coefficients["insomnia", "Estimate"]
psm_se_hyp["psm_cov_unadj_svy_isw_hyp"] <- summary(fit_psm_cov_unadj_svy_isw_hyp)$coefficients["insomnia", "Std. Error"]
psm_ci_hyp["psm_cov_unadj_svy_isw_hyp", ] <- confint(fit_psm_cov_unadj_svy_isw_hyp)["insomnia",]
psm_nobs_hyp["psm_cov_unadj_svy_isw_hyp"] <- nobs(fit_psm_cov_unadj_svy_isw_hyp)

# adjusted regression without weights
fit_psm_cov_adj_hyp <- glm(
  formula_adj, 
  family = poisson(), 
  data = df_matches_psm_cov,
  subset = !is.na(hypertension2_aha_v2) & hypertension2_aha == 0
)
psm_effect_hyp["psm_cov_adj_hyp"] <- summary(fit_psm_cov_adj_hyp)$coefficients["insomnia", "Estimate"]
psm_se_hyp["psm_cov_adj_hyp"] <- summary(fit_psm_cov_adj_hyp)$coefficients["insomnia", "Std. Error"]
psm_ci_hyp["psm_cov_adj_hyp", ] <- confint(fit_psm_cov_adj_hyp)["insomnia",]
psm_nobs_hyp["psm_cov_adj_hyp"] <- nobs(fit_psm_cov_adj_hyp)

# survey adjusted regression using OSW
fit_psm_cov_adj_svy_osw_hyp <- svyglm(
  formula_adj, 
  family = poisson(), 
  design = svy_matches_psm_cov_osw,
  subset = !is.na(hypertension2_aha_v2) & hypertension2_aha == 0
)
psm_effect_hyp["psm_cov_adj_svy_osw_hyp"] <- summary(fit_psm_cov_adj_svy_osw_hyp)$coefficients["insomnia", "Estimate"]
psm_se_hyp["psm_cov_adj_svy_osw_hyp"] <- summary(fit_psm_cov_adj_svy_osw_hyp)$coefficients["insomnia", "Std. Error"]
psm_ci_hyp["psm_cov_adj_svy_osw_hyp", ] <- confint(fit_psm_cov_adj_svy_osw_hyp)["insomnia",]
psm_nobs_hyp["psm_cov_adj_svy_osw_hyp"] <- nobs(fit_psm_cov_adj_svy_osw_hyp)

# survey adjusted regression using ISW
fit_psm_cov_adj_svy_isw_hyp <- svyglm(
  formula_adj, 
  family = poisson(), 
  design = svy_matches_psm_cov_isw,
  subset = !is.na(hypertension2_aha_v2) & hypertension2_aha == 0
)
psm_effect_hyp["psm_cov_adj_svy_isw_hyp"] <- summary(fit_psm_cov_adj_svy_isw_hyp)$coefficients["insomnia", "Estimate"]
psm_se_hyp["psm_cov_adj_svy_isw_hyp"] <- summary(fit_psm_cov_adj_svy_isw_hyp)$coefficients["insomnia", "Std. Error"]
psm_ci_hyp["psm_cov_adj_svy_isw_hyp", ] <- confint(fit_psm_cov_adj_svy_isw_hyp)["insomnia",]
psm_nobs_hyp["psm_cov_adj_svy_isw_hyp"] <- nobs(fit_psm_cov_adj_svy_isw_hyp)



