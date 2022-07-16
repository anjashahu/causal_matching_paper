
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
hyp <- hchs %>%
  select(-c(weight_norm_overall_inca, mci, mci3)) %>%
  drop_na(insomnia, age, bmi, gender, marital_status, employed, alcohol_use, cigarette_use, bkgrd1_c7, education_c3, weight_final_norm_overall, yrs_btwn_v1v2)
apply(hyp, 2, function(col) sum(is.na(col)))

# create categorical versions of continious variables 
# set cut points to categorize continuous variables 
age_cuts <- unname(quantile(hyp$age))
bmi_cuts <- c(min(hyp$bmi), 18.5, 25, 30, max(hyp$bmi))
yrs_btwn_v1v2_cuts <- unname(quantile(hyp$yrs_btwn_v1v2))
weight_final_norm_overall_cuts <- unname(quantile(hyp$weight_final_norm_overall))
# add categorical variables to individual population data frame
hyp <- hyp %>% mutate(
  age_cat = cut(
    age, 
    age_cuts, 
    include.lowest = TRUE,
    right = FALSE
  ), 
  bmi_cat = cut(
    bmi, 
    bmi_cuts, 
    include.lowest = TRUE, 
    right = FALSE
  ), 
  weight_final_norm_overall_cat = cut(
    weight_final_norm_overall, weight_final_norm_overall_cuts, 
    include.lowest = TRUE, 
    right = FALSE
  ),
  yrs_btwn_v1v2_cat =  cut(
    yrs_btwn_v1v2, yrs_btwn_v1v2_cuts, 
    include.lowest = TRUE, 
    right = FALSE
  )
)
### CEM ----------------------------------------------------------------------

# create formula to use in matching
formula_cem_cov <- as.formula("insomnia ~ age_cat + bmi_cat + gender + marital_status + employed + alcohol_use + cigarette_use + bkgrd1_c7 + education_c3 + yrs_btwn_v1v2_cat")
formula_cem_weighted <- as.formula("insomnia ~ age_cat + bmi_cat + gender + marital_status + employed + alcohol_use + cigarette_use + bkgrd1_c7 + education_c3 + yrs_btwn_v1v2_cat + weight_final_norm_overall_cat")

# when CEM is based on binning of coarsened covariates 
df_matches_cem_cov <- matchit(
  formula_cem_cov,
  method = "cem", # CEM matching
  estimand = "ATT", 
  data = hyp
) %>% 
  # obtain CEM weights (CEMW)
  match.data(weights = "cem_weights") %>%
  # obtain CEMW x OSW
  mutate(cem_osw_weights = weight_final_norm_overall * cem_weights)

# when CEM is based on binning of coarsened covariates + coarsened OSW
df_matches_cem_weighted <- matchit(
  formula_cem_weighted,
  method = "cem", # CEM matching
  estimand = "ATT", 
  data = hyp
) %>% 
  # obtain CEMW
  match.data(weights = "cem_weights") %>%
  # obtain CEMW x OSW
  mutate(cem_osw_weights = weight_final_norm_overall * cem_weights)

### CREATE SURVEY DESIGN OBJECT ---------------------------------------------------------------------------------

# for CEM matches using CEMW
# where bins are defined by coarsened covariates
svy_matches_cem_cov_cemw <- svydesign(
  ids = ~ psu_id, 
  weights = ~ cem_weights, 
  strata = ~ strat, 
  data = df_matches_cem_cov
)
# where bins are defined by coarsened covariates + coarsened OSW
svy_matches_cem_weighted_cemw <- svydesign(
  ids = ~ psu_id, 
  weights = ~ cem_weights, 
  strata = ~ strat, 
  data = df_matches_cem_weighted
)

# for CEM matches using CEMW x OSW
svy_matches_cem_cov_cemw_osw <- svydesign(
  ids = ~ psu_id, 
  weights = ~ cem_osw_weights, 
  strata = ~ strat, 
  data = df_matches_cem_cov
)
svy_matches_cem_weighted_cemw_osw <- svydesign(
  ids = ~ psu_id, 
  weights = ~ cem_osw_weights, 
  strata = ~ strat, 
  data = df_matches_cem_weighted
)

### CALCULATE EFFECTS ----------------------------------------------------------

options(survey.lonely.psu = "adjust")

# list the type of analyses we will conduct for MCI for CEM
cem_analyses_hyp <- c(
  "cem_cov_cemw_unadj_svy_hyp", "cem_cov_cemw_osw_unadj_svy_hyp", 
  "cem_cov_cemw_adj_svy_hyp", "cem_cov_cemw_osw_adj_svy_hyp",
  "cem_weighted_cemw_unadj_svy_hyp", "cem_weighted_cemw_osw_unadj_svy_hyp", 
  "cem_weighted_cemw_adj_svy_hyp", "cem_weighted_cemw_osw_adj_svy_hyp"
)

# create vector to hold effect estimates from the various analyses
cem_effect_hyp <- rep(NA, length(cem_analyses_hyp))
names(cem_effect_hyp) <- cem_analyses_hyp
# create vector to hold SEs from the various analyses
cem_se_hyp <- rep(NA, length(cem_analyses_hyp))
names(cem_se_hyp) <- cem_analyses_hyp
# matrix to hold CIs
cem_ci_hyp <- matrix(NA, nrow = length(cem_analyses_hyp), ncol = 2, dimnames = list(cem_analyses_hyp, c("lower", "upper")))
# create vector to hold number of observations 
cem_nobs_hyp <- rep(NA, length(cem_analyses_hyp))
names(cem_nobs_hyp) <- cem_analyses_hyp

# create formulas for fitting regression models 
formula_unadj <- as.formula("hypertension2_aha_v2 ~ insomnia + offset(log(yrs_btwn_v1v2))")
formula_adj <- as.formula("hypertension2_aha_v2 ~ insomnia + age + bmi + gender + marital_status + employed + alcohol_use + cigarette_use + bkgrd1_c7 + education_c3 + offset(log(yrs_btwn_v1v2))")

# when bins are defined by coarsened covariates
# unadjusted survey regression with CEMW
fit_cem_cov_cemw_unadj_svy_hyp <- svyglm(
  formula_unadj, 
  family = poisson(), 
  design = svy_matches_cem_cov_cemw,
  subset = !is.na(hypertension2_aha_v2) & hypertension2_aha == 0
)
cem_effect_hyp["cem_cov_cemw_unadj_svy_hyp"] <- summary(fit_cem_cov_cemw_unadj_svy_hyp)$coefficients["insomnia", "Estimate"]
cem_se_hyp["cem_cov_cemw_unadj_svy_hyp"] <- summary(fit_cem_cov_cemw_unadj_svy_hyp)$coefficients["insomnia", "Std. Error"]
cem_ci_hyp["cem_cov_cemw_unadj_svy_hyp", ] <- confint(fit_cem_cov_cemw_unadj_svy_hyp)["insomnia",]
cem_nobs_hyp["cem_cov_cemw_unadj_svy_hyp"] <- nobs(fit_cem_cov_cemw_unadj_svy_hyp)

# unadjusted survey regression with CEMW x OSW
fit_cem_cov_cemw_osw_unadj_svy_hyp <- svyglm(
  formula_unadj, 
  family = poisson(), 
  design = svy_matches_cem_cov_cemw_osw,
  subset = !is.na(hypertension2_aha_v2) & hypertension2_aha == 0
)
cem_effect_hyp["cem_cov_cemw_osw_unadj_svy_hyp"] <- summary(fit_cem_cov_cemw_osw_unadj_svy_hyp)$coefficients["insomnia", "Estimate"]
cem_se_hyp["cem_cov_cemw_osw_unadj_svy_hyp"] <- summary(fit_cem_cov_cemw_osw_unadj_svy_hyp)$coefficients["insomnia", "Std. Error"]
cem_ci_hyp["cem_cov_cemw_osw_unadj_svy_hyp", ] <- confint(fit_cem_cov_cemw_osw_unadj_svy_hyp)["insomnia",]
cem_nobs_hyp["cem_cov_cemw_osw_unadj_svy_hyp"] <- nobs(fit_cem_cov_cemw_osw_unadj_svy_hyp)

# adjusted survey regression with CEMW
fit_cem_cov_cemw_adj_svy_hyp <- svyglm(
  formula_adj, 
  family = poisson(), 
  design = svy_matches_cem_cov_cemw,
  subset = !is.na(hypertension2_aha_v2) & hypertension2_aha == 0
)
cem_effect_hyp["cem_cov_cemw_adj_svy_hyp"] <- summary(fit_cem_cov_cemw_adj_svy_hyp)$coefficients["insomnia", "Estimate"]
cem_se_hyp["cem_cov_cemw_adj_svy_hyp"] <- summary(fit_cem_cov_cemw_adj_svy_hyp)$coefficients["insomnia", "Std. Error"]
cem_ci_hyp["cem_cov_cemw_adj_svy_hyp", ] <- confint(fit_cem_cov_cemw_adj_svy_hyp)["insomnia",]
cem_nobs_hyp["cem_cov_cemw_adj_svy_hyp"] <- nobs(fit_cem_cov_cemw_adj_svy_hyp)

# adjusted survey regression with CEMW x OSW
fit_cem_cov_cemw_osw_adj_svy_hyp <- svyglm(
  formula_adj, 
  family = poisson(), 
  design = svy_matches_cem_cov_cemw_osw,
  subset = !is.na(hypertension2_aha_v2) & hypertension2_aha == 0
)
cem_effect_hyp["cem_cov_cemw_osw_adj_svy_hyp"] <- summary(fit_cem_cov_cemw_osw_adj_svy_hyp)$coefficients["insomnia", "Estimate"]
cem_se_hyp["cem_cov_cemw_osw_adj_svy_hyp"] <- summary(fit_cem_cov_cemw_osw_adj_svy_hyp)$coefficients["insomnia", "Std. Error"]
cem_ci_hyp["cem_cov_cemw_osw_adj_svy_hyp", ] <- confint(fit_cem_cov_cemw_osw_adj_svy_hyp)["insomnia",]
cem_nobs_hyp["cem_cov_cemw_osw_adj_svy_hyp"] <- nobs(fit_cem_cov_cemw_osw_adj_svy_hyp)

# when bins are defined by coarsened covariates + coarsened OSW
# unadjusted survey regression with CEMW
fit_cem_weighted_cemw_unadj_svy_hyp <- svyglm(
  formula_unadj, 
  family = poisson(), 
  design = svy_matches_cem_weighted_cemw,
  subset = !is.na(hypertension2_aha_v2) & hypertension2_aha == 0
)
cem_effect_hyp["cem_weighted_cemw_unadj_svy_hyp"] <- summary(fit_cem_weighted_cemw_unadj_svy_hyp)$coefficients["insomnia", "Estimate"]
cem_se_hyp["cem_weighted_cemw_unadj_svy_hyp"] <- summary(fit_cem_weighted_cemw_unadj_svy_hyp)$coefficients["insomnia", "Std. Error"]
cem_ci_hyp["cem_weighted_cemw_unadj_svy_hyp", ] <- confint(fit_cem_weighted_cemw_unadj_svy_hyp)["insomnia",]
cem_nobs_hyp["cem_weighted_cemw_unadj_svy_hyp"] <- nobs(fit_cem_weighted_cemw_unadj_svy_hyp)

# unadjusted survey regression with CEMW x OSW
fit_cem_weighted_cemw_osw_unadj_svy_hyp <- svyglm(
  formula_unadj, 
  family = poisson(), 
  design = svy_matches_cem_weighted_cemw_osw,
  subset = !is.na(hypertension2_aha_v2) & hypertension2_aha == 0
)
cem_effect_hyp["cem_weighted_cemw_osw_unadj_svy_hyp"] <- summary(fit_cem_weighted_cemw_osw_unadj_svy_hyp)$coefficients["insomnia", "Estimate"]
cem_se_hyp["cem_weighted_cemw_osw_unadj_svy_hyp"] <- summary(fit_cem_weighted_cemw_osw_unadj_svy_hyp)$coefficients["insomnia", "Std. Error"]
cem_ci_hyp["cem_weighted_cemw_osw_unadj_svy_hyp", ] <- confint(fit_cem_weighted_cemw_osw_unadj_svy_hyp)["insomnia",]
cem_nobs_hyp["cem_weighted_cemw_osw_unadj_svy_hyp"] <- nobs(fit_cem_weighted_cemw_osw_unadj_svy_hyp)

# adjusted survey regression with CEMW
fit_cem_weighted_cemw_adj_svy_hyp <- svyglm(
  formula_adj, 
  family = poisson(), 
  design = svy_matches_cem_weighted_cemw,
  subset = !is.na(hypertension2_aha_v2) & hypertension2_aha == 0
)
cem_effect_hyp["cem_weighted_cemw_adj_svy_hyp"] <- summary(fit_cem_weighted_cemw_adj_svy_hyp)$coefficients["insomnia", "Estimate"]
cem_se_hyp["cem_weighted_cemw_adj_svy_hyp"] <- summary(fit_cem_weighted_cemw_adj_svy_hyp)$coefficients["insomnia", "Std. Error"]
cem_ci_hyp["cem_weighted_cemw_adj_svy_hyp", ] <- confint(fit_cem_weighted_cemw_adj_svy_hyp)["insomnia",]
cem_nobs_hyp["cem_weighted_cemw_adj_svy_hyp"] <- nobs(fit_cem_weighted_cemw_adj_svy_hyp)

# adjusted survey regression with CEMW x OSW
fit_cem_weighted_cemw_osw_adj_svy_hyp <- svyglm(
  formula_adj, 
  family = poisson(), 
  design = svy_matches_cem_weighted_cemw_osw,
  subset = !is.na(hypertension2_aha_v2) & hypertension2_aha == 0
)
cem_effect_hyp["cem_weighted_cemw_osw_adj_svy_hyp"] <- summary(fit_cem_weighted_cemw_osw_adj_svy_hyp)$coefficients["insomnia", "Estimate"]
cem_se_hyp["cem_weighted_cemw_osw_adj_svy_hyp"] <- summary(fit_cem_weighted_cemw_osw_adj_svy_hyp)$coefficients["insomnia", "Std. Error"]
cem_ci_hyp["cem_weighted_cemw_osw_adj_svy_hyp", ] <- confint(fit_cem_weighted_cemw_osw_adj_svy_hyp)["insomnia",]
cem_nobs_hyp["cem_weighted_cemw_osw_adj_svy_hyp"] <- nobs(fit_cem_weighted_cemw_osw_adj_svy_hyp)



