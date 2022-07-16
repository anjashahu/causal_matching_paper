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

### WEIGHTING ----------------------------------------------------------------------

# create formula to use in propensity score calculation 
formula_psm_weighted <- as.formula("insomnia ~ age + bmi + gender + marital_status + employed + alcohol_use + cigarette_use + bkgrd1_c7 + education_c3")
formula_psm_cov <- as.formula("insomnia ~ age + bmi + gender + marital_status + employed + alcohol_use + cigarette_use + bkgrd1_c7 + education_c3 + weight_norm_overall_inca")

# fit logistic model to obtain propensity scores
# based on based on weighted logistic regression using OSW
fit_ps_weighted <- glm(formula_psm_weighted, weights = weight_norm_overall_inca, family = binomial(), data = mci)
# based on logistic regression using OSW as covariate
fit_ps_cov <- glm(formula_psm_cov, family = binomial(), data = mci)

# obtain propensity scores weights (PSW) using both ways for calculation (IPTW and weighting by the odds)
# and add to sample data frame and calculate the product of PSW x OSW for both
# for weighted logistic regression using OSW
mci <- mci %>%
  mutate(
    # IPTW
    ps_weighted = predict(fit_ps_weighted, type = "response"),
    psw_weighted_iptw_weights = (insomnia / ps_weighted) + ((1 - insomnia) / (1 - ps_weighted)), 
    psw_weighted_osw_iptw_weights = psw_weighted_iptw_weights * weight_norm_overall_inca,
    # weighting by the odds
    psw_weighted_odds_weights = insomnia + ((1 - insomnia) * ps_weighted / (1 - ps_weighted)),
    psw_weighted_osw_odds_weights = psw_weighted_odds_weights * weight_norm_overall_inca
  )
# for logistic regression using OSW as covariate
mci <- mci %>%
  mutate(
    # IPTW
    ps_cov = predict(fit_ps_cov, type = "response"),
    psw_cov_iptw_weights = (insomnia / ps_cov) + ((1 - insomnia) / (1 - ps_cov)), 
    psw_cov_osw_iptw_weights = psw_cov_iptw_weights * weight_norm_overall_inca,
    # weighting by the odds
    psw_cov_odds_weights = insomnia + ((1 - insomnia) * ps_cov / (1 - ps_cov)),
    psw_cov_osw_odds_weights = psw_cov_odds_weights * weight_norm_overall_inca
  )

### CREATE SURVEY DESIGN OBJECT ---------------------------------------------------------------------------------

# for weighting using PSW

# where propensity score is calculated using weighted logistic regression 
# for IPTW
svy_psw_weighted_iptw <- svydesign(
  ids = ~ psu_id, 
  weights = ~ psw_weighted_iptw_weights, 
  strata = ~ strat, 
  data = mci
)  
# for weighting by the odds
svy_psw_weighted_odds <- svydesign(
  ids = ~ psu_id, 
  weights = ~ psw_weighted_odds_weights, 
  strata = ~ strat, 
  data = mci
)

# where propensity score is calculated using logistic regression with weight as covariate
# for IPTW
svy_psw_cov_iptw <- svydesign(
  ids = ~ psu_id, 
  weights = ~ psw_cov_iptw_weights, 
  strata = ~ strat, 
  data = mci
) 
# for weighting by the odds
svy_psw_cov_odds <- svydesign(
  ids = ~ psu_id, 
  weights = ~ psw_cov_odds_weights, 
  strata = ~ strat, 
  data = mci
) 


# for weighting using PSW x OSW

# where propensity score is calculated using weighted logistic regression 
# for IPTW
svy_psw_weighted_osw_iptw <- svydesign(
  ids = ~ psu_id, 
  weights = ~ psw_weighted_osw_iptw_weights, 
  strata = ~ strat, 
  data = mci
) 
# for weighting by the odds
svy_psw_weighted_osw_odds <- svydesign(
  ids = ~ psu_id, 
  weights = ~ psw_weighted_osw_odds_weights, 
  strata = ~ strat, 
  data = mci
) 

# where propensity score is calculated using logistic regression with weight as covariate
# for IPTW
svy_psw_cov_osw_iptw <- svydesign(
  ids = ~ psu_id,
  weights = ~ psw_cov_osw_iptw_weights, 
  strata = ~ strat, 
  data = mci
) 
# for weighting by the odds
svy_psw_cov_osw_odds <- svydesign(
  ids = ~ psu_id, 
  weights = ~ psw_cov_osw_odds_weights, 
  strata = ~ strat, 
  data = mci
) 

### CALCULATE EFFECTS ----------------------------------------------------------

options(survey.lonely.psu = "adjust")


# list the type of analyses we will conduct for MCI for weighting
weighting_analyses_mci <- c(
  "psw_weighted_iptw_unadj_svy_mci", "psw_weighted_iptw_osw_unadj_svy_mci", 
  "psw_weighted_odds_unadj_svy_mci", "psw_weighted_odds_osw_unadj_svy_mci",
  "psw_weighted_iptw_adj_svy_mci", "psw_weighted_iptw_osw_adj_svy_mci", 
  "psw_weighted_odds_adj_svy_mci", "psw_weighted_odds_osw_adj_svy_mci",
  "psw_cov_iptw_unadj_svy_mci", "psw_cov_iptw_osw_unadj_svy_mci", 
  "psw_cov_odds_unadj_svy_mci", "psw_cov_odds_osw_unadj_svy_mci",
  "psw_cov_iptw_adj_svy_mci", "psw_cov_iptw_osw_adj_svy_mci", 
  "psw_cov_odds_adj_svy_mci", "psw_cov_odds_osw_adj_svy_mci"
)

# create vector to hold effect estimates from the various analyses
weighting_effect_mci <- rep(NA, length(weighting_analyses_mci))
names(weighting_effect_mci) <- weighting_analyses_mci
# create vector to hold SEs from the various analyses
weighting_se_mci <- rep(NA, length(weighting_analyses_mci))
names(weighting_se_mci) <- weighting_analyses_mci
# matrix to hold CIs
weighting_ci_mci <- matrix(NA, nrow = length(weighting_analyses_mci), ncol = 2, dimnames = list(weighting_analyses_mci, c("lower", "upper")))
# create vector to hold number of observations 
weighting_nobs_mci <- rep(NA, length(weighting_analyses_mci))
names(weighting_nobs_mci) <- weighting_analyses_mci

# create formulas for fitting regression models 
formula_unadj <- as.formula("mci ~ insomnia")
formula_adj <- as.formula("mci ~ insomnia + age + bmi + gender + marital_status + employed + alcohol_use + cigarette_use + bkgrd1_c7 + education_c3")


# fit the various regressions

# when propensity score is calculated using weighted logistic regression 
# unadjusted survey regression with PSW using IPTW
fit_psw_weighted_iptw_unadj_svy_mci <- svyglm(
  formula_unadj, 
  family = quasibinomial(), 
  design = svy_psw_weighted_iptw,
  subset = !is.na(mci)
)
weighting_effect_mci["psw_weighted_iptw_unadj_svy_mci"] <- summary(fit_psw_weighted_iptw_unadj_svy_mci)$coefficients["insomnia", "Estimate"]
weighting_se_mci["psw_weighted_iptw_unadj_svy_mci"] <- summary(fit_psw_weighted_iptw_unadj_svy_mci)$coefficients["insomnia", "Std. Error"]
weighting_ci_mci["psw_weighted_iptw_unadj_svy_mci", ] <- confint(fit_psw_weighted_iptw_unadj_svy_mci)["insomnia",]
weighting_nobs_mci["psw_weighted_iptw_unadj_svy_mci"] <- nobs(fit_psw_weighted_iptw_unadj_svy_mci)

# unadjusted survey regression with PSW x OSW using IPTW
fit_psw_weighted_iptw_osw_unadj_svy_mci <- svyglm(
  formula_unadj, 
  family = quasibinomial(), 
  design = svy_psw_weighted_osw_iptw,
  subset = !is.na(mci)
)
weighting_effect_mci["psw_weighted_iptw_osw_unadj_svy_mci"] <- summary(fit_psw_weighted_iptw_osw_unadj_svy_mci)$coefficients["insomnia", "Estimate"]
weighting_se_mci["psw_weighted_iptw_osw_unadj_svy_mci"] <- summary(fit_psw_weighted_iptw_osw_unadj_svy_mci)$coefficients["insomnia", "Std. Error"]
weighting_ci_mci["psw_weighted_iptw_osw_unadj_svy_mci", ] <- confint(fit_psw_weighted_iptw_osw_unadj_svy_mci)["insomnia",]
weighting_nobs_mci["psw_weighted_iptw_osw_unadj_svy_mci"] <- nobs(fit_psw_weighted_iptw_osw_unadj_svy_mci)

# unadjusted survey regression with PSW using weighting by the odds
fit_psw_weighted_odds_unadj_svy_mci <-  svyglm(
  formula_unadj, 
  family = quasibinomial(), 
  design = svy_psw_weighted_odds,
  subset = !is.na(mci)
)
weighting_effect_mci["psw_weighted_odds_unadj_svy_mci"] <- summary(fit_psw_weighted_odds_unadj_svy_mci)$coefficients["insomnia", "Estimate"]
weighting_se_mci["psw_weighted_odds_unadj_svy_mci"] <- summary(fit_psw_weighted_odds_unadj_svy_mci)$coefficients["insomnia", "Std. Error"]
weighting_ci_mci["psw_weighted_odds_unadj_svy_mci", ] <- confint(fit_psw_weighted_odds_unadj_svy_mci)["insomnia",]
weighting_nobs_mci["psw_weighted_odds_unadj_svy_mci"] <- nobs(fit_psw_weighted_odds_unadj_svy_mci)

# unadjusted survey regression with PSW x OSW using weighting by the odds
fit_psw_weighted_odds_osw_unadj_svy_mci <-  svyglm(
  formula_unadj, 
  family = quasibinomial(), 
  design = svy_psw_weighted_osw_odds,
  subset = !is.na(mci)
)
weighting_effect_mci["psw_weighted_odds_osw_unadj_svy_mci"] <- summary(fit_psw_weighted_odds_osw_unadj_svy_mci)$coefficients["insomnia", "Estimate"]
weighting_se_mci["psw_weighted_odds_osw_unadj_svy_mci"] <- summary(fit_psw_weighted_odds_osw_unadj_svy_mci)$coefficients["insomnia", "Std. Error"]
weighting_ci_mci["psw_weighted_odds_osw_unadj_svy_mci", ] <- confint(fit_psw_weighted_odds_osw_unadj_svy_mci)["insomnia",]
weighting_nobs_mci["psw_weighted_odds_osw_unadj_svy_mci"] <- nobs(fit_psw_weighted_odds_osw_unadj_svy_mci)

# adjusted survey regression with PSW using IPTW
fit_psw_weighted_iptw_adj_svy_mci <- svyglm(
  formula_adj, 
  family = quasibinomial(), 
  design = svy_psw_weighted_iptw,
  subset = !is.na(mci)
)
weighting_effect_mci["psw_weighted_iptw_adj_svy_mci"] <- summary(fit_psw_weighted_iptw_adj_svy_mci)$coefficients["insomnia", "Estimate"]
weighting_se_mci["psw_weighted_iptw_adj_svy_mci"] <- summary(fit_psw_weighted_iptw_adj_svy_mci)$coefficients["insomnia", "Std. Error"]
weighting_ci_mci["psw_weighted_iptw_adj_svy_mci", ] <- confint(fit_psw_weighted_iptw_adj_svy_mci)["insomnia",]
weighting_nobs_mci["psw_weighted_iptw_adj_svy_mci"] <- nobs(fit_psw_weighted_iptw_adj_svy_mci)

# adjusted survey regression with PSW x OSW using IPTW
fit_psw_weighted_iptw_osw_adj_svy_mci <- svyglm(
  formula_adj, 
  family = quasibinomial(), 
  design = svy_psw_weighted_osw_iptw,
  subset = !is.na(mci)
)
weighting_effect_mci["psw_weighted_iptw_osw_adj_svy_mci"] <- summary(fit_psw_weighted_iptw_osw_adj_svy_mci)$coefficients["insomnia", "Estimate"]
weighting_se_mci["psw_weighted_iptw_osw_adj_svy_mci"] <- summary(fit_psw_weighted_iptw_osw_adj_svy_mci)$coefficients["insomnia", "Std. Error"]
weighting_ci_mci["psw_weighted_iptw_osw_adj_svy_mci", ] <- confint(fit_psw_weighted_iptw_osw_adj_svy_mci)["insomnia",]
weighting_nobs_mci["psw_weighted_iptw_osw_adj_svy_mci"] <- nobs(fit_psw_weighted_iptw_osw_adj_svy_mci)

# adjusted survey regression with PSW using weighting by the odds
fit_psw_weighted_odds_adj_svy_mci <-  svyglm(
  formula_adj, 
  family = quasibinomial(), 
  design = svy_psw_weighted_odds,
  subset = !is.na(mci)
)
weighting_effect_mci["psw_weighted_odds_adj_svy_mci"] <- summary(fit_psw_weighted_odds_adj_svy_mci)$coefficients["insomnia", "Estimate"]
weighting_se_mci["psw_weighted_odds_adj_svy_mci"] <- summary(fit_psw_weighted_odds_adj_svy_mci)$coefficients["insomnia", "Std. Error"]
weighting_ci_mci["psw_weighted_odds_adj_svy_mci", ] <- confint(fit_psw_weighted_odds_adj_svy_mci)["insomnia",]
weighting_nobs_mci["psw_weighted_odds_adj_svy_mci"] <- nobs(fit_psw_weighted_odds_adj_svy_mci)

# adjusted survey regression with PSW x OSW using weighting by the odds
fit_psw_weighted_odds_osw_adj_svy_mci <-  svyglm(
  formula_adj, 
  family = quasibinomial(), 
  design = svy_psw_weighted_osw_odds,
  subset = !is.na(mci)
)
weighting_effect_mci["psw_weighted_odds_osw_adj_svy_mci"] <- summary(fit_psw_weighted_odds_osw_adj_svy_mci)$coefficients["insomnia", "Estimate"]
weighting_se_mci["psw_weighted_odds_osw_adj_svy_mci"] <- summary(fit_psw_weighted_odds_osw_adj_svy_mci)$coefficients["insomnia", "Std. Error"]
weighting_ci_mci["psw_weighted_odds_osw_adj_svy_mci", ] <- confint(fit_psw_weighted_odds_osw_adj_svy_mci)["insomnia",]
weighting_nobs_mci["psw_weighted_odds_osw_adj_svy_mci"] <- nobs(fit_psw_weighted_odds_osw_adj_svy_mci)

# when propensity score is calculated using logistic regression with weight as covariate
# unadjusted survey regression with PSW using IPTW
fit_psw_cov_iptw_unadj_svy_mci <- svyglm(
  formula_unadj, 
  family = quasibinomial(), 
  design = svy_psw_cov_iptw,
  subset = !is.na(mci)
)
weighting_effect_mci["psw_cov_iptw_unadj_svy_mci"] <- summary(fit_psw_cov_iptw_unadj_svy_mci)$coefficients["insomnia", "Estimate"]
weighting_se_mci["psw_cov_iptw_unadj_svy_mci"] <- summary(fit_psw_cov_iptw_unadj_svy_mci)$coefficients["insomnia", "Std. Error"]
weighting_ci_mci["psw_cov_iptw_unadj_svy_mci", ] <- confint(fit_psw_cov_iptw_unadj_svy_mci)["insomnia",]
weighting_nobs_mci["psw_cov_iptw_unadj_svy_mci"] <- nobs(fit_psw_cov_iptw_unadj_svy_mci)

# unadjusted survey regression with PSW x OSW using IPTW
fit_psw_cov_iptw_osw_unadj_svy_mci <- svyglm(
  formula_unadj, 
  family = quasibinomial(), 
  design = svy_psw_cov_osw_iptw,
  subset = !is.na(mci)
)
weighting_effect_mci["psw_cov_iptw_osw_unadj_svy_mci"] <- summary(fit_psw_cov_iptw_osw_unadj_svy_mci)$coefficients["insomnia", "Estimate"]
weighting_se_mci["psw_cov_iptw_osw_unadj_svy_mci"] <- summary(fit_psw_cov_iptw_osw_unadj_svy_mci)$coefficients["insomnia", "Std. Error"]
weighting_ci_mci["psw_cov_iptw_osw_unadj_svy_mci", ] <- confint(fit_psw_cov_iptw_osw_unadj_svy_mci)["insomnia",]
weighting_nobs_mci["psw_cov_iptw_osw_unadj_svy_mci"] <- nobs(fit_psw_cov_iptw_osw_unadj_svy_mci)

# unadjusted survey regression with PSW using weighting by the odds
fit_psw_cov_odds_unadj_svy_mci <-  svyglm(
  formula_unadj, 
  family = quasibinomial(), 
  design = svy_psw_cov_odds,
  subset = !is.na(mci)
)
weighting_effect_mci["psw_cov_odds_unadj_svy_mci"] <- summary(fit_psw_cov_odds_unadj_svy_mci)$coefficients["insomnia", "Estimate"]
weighting_se_mci["psw_cov_odds_unadj_svy_mci"] <- summary(fit_psw_cov_odds_unadj_svy_mci)$coefficients["insomnia", "Std. Error"]
weighting_ci_mci["psw_cov_odds_unadj_svy_mci", ] <- confint(fit_psw_cov_odds_unadj_svy_mci)["insomnia",]
weighting_nobs_mci["psw_cov_odds_unadj_svy_mci"] <- nobs(fit_psw_cov_odds_unadj_svy_mci)

# unadjusted survey regression with PSW x OSW using weighting by the odds
fit_psw_cov_odds_osw_unadj_svy_mci <-  svyglm(
  formula_unadj, 
  family = quasibinomial(), 
  design = svy_psw_cov_osw_odds,
  subset = !is.na(mci)
)
weighting_effect_mci["psw_cov_odds_osw_unadj_svy_mci"] <- summary(fit_psw_cov_odds_osw_unadj_svy_mci)$coefficients["insomnia", "Estimate"]
weighting_se_mci["psw_cov_odds_osw_unadj_svy_mci"] <- summary(fit_psw_cov_odds_osw_unadj_svy_mci)$coefficients["insomnia", "Std. Error"]
weighting_ci_mci["psw_cov_odds_osw_unadj_svy_mci", ] <- confint(fit_psw_cov_odds_osw_unadj_svy_mci)["insomnia",]
weighting_nobs_mci["psw_cov_odds_osw_unadj_svy_mci"] <- nobs(fit_psw_cov_odds_osw_unadj_svy_mci)

# adjusted survey regression with PSW using IPTW
fit_psw_cov_iptw_adj_svy_mci <- svyglm(
  formula_adj, 
  family = quasibinomial(), 
  design = svy_psw_cov_iptw,
  subset = !is.na(mci)
)
weighting_effect_mci["psw_cov_iptw_adj_svy_mci"] <- summary(fit_psw_cov_iptw_adj_svy_mci)$coefficients["insomnia", "Estimate"]
weighting_se_mci["psw_cov_iptw_adj_svy_mci"] <- summary(fit_psw_cov_iptw_adj_svy_mci)$coefficients["insomnia", "Std. Error"]
weighting_ci_mci["psw_cov_iptw_adj_svy_mci", ] <- confint(fit_psw_cov_iptw_adj_svy_mci)["insomnia",]
weighting_nobs_mci["psw_cov_iptw_adj_svy_mci"] <- nobs(fit_psw_cov_iptw_adj_svy_mci)

# adjusted survey regression with PSW x OSW using IPTW
fit_psw_cov_iptw_osw_adj_svy_mci <- svyglm(
  formula_adj, 
  family = quasibinomial(), 
  design = svy_psw_cov_osw_iptw,
  subset = !is.na(mci)
)
weighting_effect_mci["psw_cov_iptw_osw_adj_svy_mci"] <- summary(fit_psw_cov_iptw_osw_adj_svy_mci)$coefficients["insomnia", "Estimate"]
weighting_se_mci["psw_cov_iptw_osw_adj_svy_mci"] <- summary(fit_psw_cov_iptw_osw_adj_svy_mci)$coefficients["insomnia", "Std. Error"]
weighting_ci_mci["psw_cov_iptw_osw_adj_svy_mci", ] <- confint(fit_psw_cov_iptw_osw_adj_svy_mci)["insomnia",]
weighting_nobs_mci["psw_cov_iptw_osw_adj_svy_mci"] <- nobs(fit_psw_cov_iptw_osw_adj_svy_mci)

# adjusted survey regression with PSW using weighting by the odds
fit_psw_cov_odds_adj_svy_mci <-  svyglm(
  formula_adj, 
  family = quasibinomial(), 
  design = svy_psw_cov_odds,
  subset = !is.na(mci)
)
weighting_effect_mci["psw_cov_odds_adj_svy_mci"] <- summary(fit_psw_cov_odds_adj_svy_mci)$coefficients["insomnia", "Estimate"]
weighting_se_mci["psw_cov_odds_adj_svy_mci"] <- summary(fit_psw_cov_odds_adj_svy_mci)$coefficients["insomnia", "Std. Error"]
weighting_ci_mci["psw_cov_odds_adj_svy_mci", ] <- confint(fit_psw_cov_odds_adj_svy_mci)["insomnia",]
weighting_nobs_mci["psw_cov_odds_adj_svy_mci"] <- nobs(fit_psw_cov_odds_adj_svy_mci)

# adjusted survey regression with PSW x OSW using weighting by the odds
fit_psw_cov_odds_osw_adj_svy_mci <-  svyglm(
  formula_adj, 
  family = quasibinomial(), 
  design = svy_psw_cov_osw_odds,
  subset = !is.na(mci)
)
weighting_effect_mci["psw_cov_odds_osw_adj_svy_mci"] <- summary(fit_psw_cov_odds_osw_adj_svy_mci)$coefficients["insomnia", "Estimate"]
weighting_se_mci["psw_cov_odds_osw_adj_svy_mci"] <- summary(fit_psw_cov_odds_osw_adj_svy_mci)$coefficients["insomnia", "Std. Error"]
weighting_ci_mci["psw_cov_odds_osw_adj_svy_mci", ] <- confint(fit_psw_cov_odds_osw_adj_svy_mci)["insomnia",]
weighting_nobs_mci["psw_cov_odds_osw_adj_svy_mci"] <- nobs(fit_psw_cov_odds_osw_adj_svy_mci)



