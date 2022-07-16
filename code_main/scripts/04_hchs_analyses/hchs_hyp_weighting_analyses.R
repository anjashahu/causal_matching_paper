
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


### WEIGHTING ----------------------------------------------------------------------

# create formula to use in propensity score calculation 
formula_psm_weighted <- as.formula("insomnia ~ age + bmi + gender + marital_status + employed + alcohol_use + cigarette_use + bkgrd1_c7 + education_c3 + yrs_btwn_v1v2")
formula_psm_cov <- as.formula("insomnia ~ age + bmi + gender + marital_status + employed + alcohol_use + cigarette_use + bkgrd1_c7 + education_c3 + yrs_btwn_v1v2 + weight_final_norm_overall")

# fit logistic model to obtain propensity scores
# based on based on weighted logistic regression using OSW
fit_ps_weighted <- glm(formula_psm_weighted, weights = weight_final_norm_overall, family = binomial(), data = hyp)
# based on logistic regression using OSW as covariate
fit_ps_cov <- glm(formula_psm_cov, family = binomial(), data = hyp)

# obtain propensity scores weights (PSW) using both ways for calculation (IPTW and weighting by the odds)
# and add to sample data frame and calculate the product of PSW x OSW for both
# for weighted logistic regression using OSW
hyp <- hyp %>%
  mutate(
    # IPTW
    ps_weighted = predict(fit_ps_weighted, type = "response"),
    psw_weighted_iptw_weights = (insomnia / ps_weighted) + ((1 - insomnia) / (1 - ps_weighted)), 
    psw_weighted_osw_iptw_weights = psw_weighted_iptw_weights * weight_final_norm_overall,
    # weighting by the odds
    psw_weighted_odds_weights = insomnia + ((1 - insomnia) * ps_weighted / (1 - ps_weighted)),
    psw_weighted_osw_odds_weights = psw_weighted_odds_weights * weight_final_norm_overall
  )
# for logistic regression using OSW as covariate
hyp <- hyp %>%
  mutate(
    # IPTW
    ps_cov = predict(fit_ps_cov, type = "response"),
    psw_cov_iptw_weights = (insomnia / ps_cov) + ((1 - insomnia) / (1 - ps_cov)), 
    psw_cov_osw_iptw_weights = psw_cov_iptw_weights * weight_final_norm_overall,
    # weighting by the odds
    psw_cov_odds_weights = insomnia + ((1 - insomnia) * ps_cov / (1 - ps_cov)),
    psw_cov_osw_odds_weights = psw_cov_odds_weights * weight_final_norm_overall
  )

### CREATE SURVEY DESIGN OBJECT ---------------------------------------------------------------------------------

# for weighting using PSW

# where propensity score is calculated using weighted logistic regression 
# for IPTW
svy_psw_weighted_iptw <- svydesign(
  ids = ~ psu_id, 
  weights = ~ psw_weighted_iptw_weights, 
  strata = ~ strat, 
  data = hyp
)  
# for weighting by the odds
svy_psw_weighted_odds <- svydesign(
  ids = ~ psu_id, 
  weights = ~ psw_weighted_odds_weights, 
  strata = ~ strat, 
  data = hyp
)

# where propensity score is calculated using logistic regression with weight as covariate
# for IPTW
svy_psw_cov_iptw <- svydesign(
  ids = ~ psu_id, 
  weights = ~ psw_cov_iptw_weights, 
  strata = ~ strat, 
  data = hyp
) 
# for weighting by the odds
svy_psw_cov_odds <- svydesign(
  ids = ~ psu_id, 
  weights = ~ psw_cov_odds_weights, 
  strata = ~ strat, 
  data = hyp
) 


# for weighting using PSW x OSW

# where propensity score is calculated using weighted logistic regression 
# for IPTW
svy_psw_weighted_osw_iptw <- svydesign(
  ids = ~ psu_id, 
  weights = ~ psw_weighted_osw_iptw_weights, 
  strata = ~ strat, 
  data = hyp
) 
# for weighting by the odds
svy_psw_weighted_osw_odds <- svydesign(
  ids = ~ psu_id, 
  weights = ~ psw_weighted_osw_odds_weights, 
  strata = ~ strat, 
  data = hyp
) 

# where propensity score is calculated using logistic regression with weight as covariate
# for IPTW
svy_psw_cov_osw_iptw <- svydesign(
  ids = ~ psu_id,
  weights = ~ psw_cov_osw_iptw_weights, 
  strata = ~ strat, 
  data = hyp
) 
# for weighting by the odds
svy_psw_cov_osw_odds <- svydesign(
  ids = ~ psu_id, 
  weights = ~ psw_cov_osw_odds_weights, 
  strata = ~ strat, 
  data = hyp
) 

### CALCULATE EFFECTS ----------------------------------------------------------

options(survey.lonely.psu = "adjust")

# list the type of analyses we will conduct for hypertension for weighting
weighting_analyses_hyp <- c(
  "psw_weighted_iptw_unadj_svy_hyp", "psw_weighted_iptw_osw_unadj_svy_hyp", 
  "psw_weighted_odds_unadj_svy_hyp", "psw_weighted_odds_osw_unadj_svy_hyp",
  "psw_weighted_iptw_adj_svy_hyp", "psw_weighted_iptw_osw_adj_svy_hyp", 
  "psw_weighted_odds_adj_svy_hyp", "psw_weighted_odds_osw_adj_svy_hyp",
  "psw_cov_iptw_unadj_svy_hyp", "psw_cov_iptw_osw_unadj_svy_hyp", 
  "psw_cov_odds_unadj_svy_hyp", "psw_cov_odds_osw_unadj_svy_hyp",
  "psw_cov_iptw_adj_svy_hyp", "psw_cov_iptw_osw_adj_svy_hyp", 
  "psw_cov_odds_adj_svy_hyp", "psw_cov_odds_osw_adj_svy_hyp"
)

# create vector to hold effect estimates from the various analyses
weighting_effect_hyp <- rep(NA, length(weighting_analyses_hyp))
names(weighting_effect_hyp) <- weighting_analyses_hyp
# create vector to hold SEs from the various analyses
weighting_se_hyp <- rep(NA, length(weighting_analyses_hyp))
names(weighting_se_hyp) <- weighting_analyses_hyp
# matrix to hold CIs
weighting_ci_hyp <- matrix(NA, nrow = length(weighting_analyses_hyp), ncol = 2, dimnames = list(weighting_analyses_hyp, c("lower", "upper")))
# create vector to hold number of observations 
weighting_nobs_hyp <- rep(NA, length(weighting_analyses_hyp))
names(weighting_nobs_hyp) <- weighting_analyses_hyp


# create formulas for fitting regression models 
formula_unadj <- as.formula("hypertension2_aha_v2 ~ insomnia + offset(log(yrs_btwn_v1v2))")
formula_adj <- as.formula("hypertension2_aha_v2 ~ insomnia + age + bmi + gender + marital_status + employed + alcohol_use + cigarette_use + bkgrd1_c7 + education_c3 + offset(log(yrs_btwn_v1v2))")


# fit the various regressions

# when propensity score is calculated using weighted logistic regression 
# unadjusted survey regression with PSW using IPTW
fit_psw_weighted_iptw_unadj_svy_hyp <- svyglm(
  formula_unadj, 
  family = poisson(),  
  design = svy_psw_weighted_iptw,
  subset = !is.na(hypertension2_aha_v2) & hypertension2_aha == 0
)
weighting_effect_hyp["psw_weighted_iptw_unadj_svy_hyp"] <- summary(fit_psw_weighted_iptw_unadj_svy_hyp)$coefficients["insomnia", "Estimate"]
weighting_se_hyp["psw_weighted_iptw_unadj_svy_hyp"] <- summary(fit_psw_weighted_iptw_unadj_svy_hyp)$coefficients["insomnia", "Std. Error"]
weighting_ci_hyp["psw_weighted_iptw_unadj_svy_hyp", ] <- confint(fit_psw_weighted_iptw_unadj_svy_hyp)["insomnia",]
weighting_nobs_hyp["psw_weighted_iptw_unadj_svy_hyp"] <- nobs(fit_psw_weighted_iptw_unadj_svy_hyp)

# unadjusted survey regression with PSW x OSW using IPTW
fit_psw_weighted_iptw_osw_unadj_svy_hyp <- svyglm(
  formula_unadj, 
  family = poisson(), 
  design = svy_psw_weighted_osw_iptw,
  subset = !is.na(hypertension2_aha_v2) & hypertension2_aha == 0
)
weighting_effect_hyp["psw_weighted_iptw_osw_unadj_svy_hyp"] <- summary(fit_psw_weighted_iptw_osw_unadj_svy_hyp)$coefficients["insomnia", "Estimate"]
weighting_se_hyp["psw_weighted_iptw_osw_unadj_svy_hyp"] <- summary(fit_psw_weighted_iptw_osw_unadj_svy_hyp)$coefficients["insomnia", "Std. Error"]
weighting_ci_hyp["psw_weighted_iptw_osw_unadj_svy_hyp", ] <- confint(fit_psw_weighted_iptw_osw_unadj_svy_hyp)["insomnia",]
weighting_nobs_hyp["psw_weighted_iptw_osw_unadj_svy_hyp"] <- nobs(fit_psw_weighted_iptw_osw_unadj_svy_hyp)

# unadjusted survey regression with PSW using weighting by the odds
fit_psw_weighted_odds_unadj_svy_hyp <-  svyglm(
  formula_unadj, 
  family = poisson(), 
  design = svy_psw_weighted_odds,
  subset = !is.na(hypertension2_aha_v2) & hypertension2_aha == 0
)
weighting_effect_hyp["psw_weighted_odds_unadj_svy_hyp"] <- summary(fit_psw_weighted_odds_unadj_svy_hyp)$coefficients["insomnia", "Estimate"]
weighting_se_hyp["psw_weighted_odds_unadj_svy_hyp"] <- summary(fit_psw_weighted_odds_unadj_svy_hyp)$coefficients["insomnia", "Std. Error"]
weighting_ci_hyp["psw_weighted_odds_unadj_svy_hyp", ] <- confint(fit_psw_weighted_odds_unadj_svy_hyp)["insomnia",]
weighting_nobs_hyp["psw_weighted_odds_unadj_svy_hyp"] <- nobs(fit_psw_weighted_odds_unadj_svy_hyp)

# unadjusted survey regression with PSW x OSW using weighting by the odds
fit_psw_weighted_odds_osw_unadj_svy_hyp <-  svyglm(
  formula_unadj, 
  family = poisson(), 
  design = svy_psw_weighted_osw_odds,
  subset = !is.na(hypertension2_aha_v2) & hypertension2_aha == 0
)
weighting_effect_hyp["psw_weighted_odds_osw_unadj_svy_hyp"] <- summary(fit_psw_weighted_odds_osw_unadj_svy_hyp)$coefficients["insomnia", "Estimate"]
weighting_se_hyp["psw_weighted_odds_osw_unadj_svy_hyp"] <- summary(fit_psw_weighted_odds_osw_unadj_svy_hyp)$coefficients["insomnia", "Std. Error"]
weighting_ci_hyp["psw_weighted_odds_osw_unadj_svy_hyp", ] <- confint(fit_psw_weighted_odds_osw_unadj_svy_hyp)["insomnia",]
weighting_nobs_hyp["psw_weighted_odds_osw_unadj_svy_hyp"] <- nobs(fit_psw_weighted_odds_osw_unadj_svy_hyp)

# adjusted survey regression with PSW using IPTW
fit_psw_weighted_iptw_adj_svy_hyp <- svyglm(
  formula_adj, 
  family = poisson(), 
  design = svy_psw_weighted_iptw,
  subset = !is.na(hypertension2_aha_v2) & hypertension2_aha == 0
)
weighting_effect_hyp["psw_weighted_iptw_adj_svy_hyp"] <- summary(fit_psw_weighted_iptw_adj_svy_hyp)$coefficients["insomnia", "Estimate"]
weighting_se_hyp["psw_weighted_iptw_adj_svy_hyp"] <- summary(fit_psw_weighted_iptw_adj_svy_hyp)$coefficients["insomnia", "Std. Error"]
weighting_ci_hyp["psw_weighted_iptw_adj_svy_hyp", ] <- confint(fit_psw_weighted_iptw_adj_svy_hyp)["insomnia",]
weighting_nobs_hyp["psw_weighted_iptw_adj_svy_hyp"] <- nobs(fit_psw_weighted_iptw_adj_svy_hyp)

# adjusted survey regression with PSW x OSW using IPTW
fit_psw_weighted_iptw_osw_adj_svy_hyp <- svyglm(
  formula_adj, 
  family = poisson(), 
  design = svy_psw_weighted_osw_iptw,
  subset = !is.na(hypertension2_aha_v2) & hypertension2_aha == 0
)
weighting_effect_hyp["psw_weighted_iptw_osw_adj_svy_hyp"] <- summary(fit_psw_weighted_iptw_osw_adj_svy_hyp)$coefficients["insomnia", "Estimate"]
weighting_se_hyp["psw_weighted_iptw_osw_adj_svy_hyp"] <- summary(fit_psw_weighted_iptw_osw_adj_svy_hyp)$coefficients["insomnia", "Std. Error"]
weighting_ci_hyp["psw_weighted_iptw_osw_adj_svy_hyp", ] <- confint(fit_psw_weighted_iptw_osw_adj_svy_hyp)["insomnia",]
weighting_nobs_hyp["psw_weighted_iptw_osw_adj_svy_hyp"] <- nobs(fit_psw_weighted_iptw_osw_adj_svy_hyp)

# adjusted survey regression with PSW using weighting by the odds
fit_psw_weighted_odds_adj_svy_hyp <-  svyglm(
  formula_adj, 
  family = poisson(), 
  design = svy_psw_weighted_odds,
  subset = !is.na(hypertension2_aha_v2) & hypertension2_aha == 0
)
weighting_effect_hyp["psw_weighted_odds_adj_svy_hyp"] <- summary(fit_psw_weighted_odds_adj_svy_hyp)$coefficients["insomnia", "Estimate"]
weighting_se_hyp["psw_weighted_odds_adj_svy_hyp"] <- summary(fit_psw_weighted_odds_adj_svy_hyp)$coefficients["insomnia", "Std. Error"]
weighting_ci_hyp["psw_weighted_odds_adj_svy_hyp", ] <- confint(fit_psw_weighted_odds_adj_svy_hyp)["insomnia",]
weighting_nobs_hyp["psw_weighted_odds_adj_svy_hyp"] <- nobs(fit_psw_weighted_odds_adj_svy_hyp)

# adjusted survey regression with PSW x OSW using weighting by the odds
fit_psw_weighted_odds_osw_adj_svy_hyp <-  svyglm(
  formula_adj, 
  family = poisson(), 
  design = svy_psw_weighted_osw_odds,
  subset = !is.na(hypertension2_aha_v2) & hypertension2_aha == 0
)
weighting_effect_hyp["psw_weighted_odds_osw_adj_svy_hyp"] <- summary(fit_psw_weighted_odds_osw_adj_svy_hyp)$coefficients["insomnia", "Estimate"]
weighting_se_hyp["psw_weighted_odds_osw_adj_svy_hyp"] <- summary(fit_psw_weighted_odds_osw_adj_svy_hyp)$coefficients["insomnia", "Std. Error"]
weighting_ci_hyp["psw_weighted_odds_osw_adj_svy_hyp", ] <- confint(fit_psw_weighted_odds_osw_adj_svy_hyp)["insomnia",]
weighting_nobs_hyp["psw_weighted_odds_osw_adj_svy_hyp"] <- nobs(fit_psw_weighted_odds_osw_adj_svy_hyp)

# when propensity score is calculated using logistic regression with weight as covariate
# unadjusted survey regression with PSW using IPTW
fit_psw_cov_iptw_unadj_svy_hyp <- svyglm(
  formula_unadj, 
  family = poisson(), 
  design = svy_psw_cov_iptw,
  subset = !is.na(hypertension2_aha_v2) & hypertension2_aha == 0
)
weighting_effect_hyp["psw_cov_iptw_unadj_svy_hyp"] <- summary(fit_psw_cov_iptw_unadj_svy_hyp)$coefficients["insomnia", "Estimate"]
weighting_se_hyp["psw_cov_iptw_unadj_svy_hyp"] <- summary(fit_psw_cov_iptw_unadj_svy_hyp)$coefficients["insomnia", "Std. Error"]
weighting_ci_hyp["psw_cov_iptw_unadj_svy_hyp", ] <- confint(fit_psw_cov_iptw_unadj_svy_hyp)["insomnia",]
weighting_nobs_hyp["psw_cov_iptw_unadj_svy_hyp"] <- nobs(fit_psw_cov_iptw_unadj_svy_hyp)

# unadjusted survey regression with PSW x OSW using IPTW
fit_psw_cov_iptw_osw_unadj_svy_hyp <- svyglm(
  formula_unadj, 
  family = poisson(), 
  design = svy_psw_cov_osw_iptw,
  subset = !is.na(hypertension2_aha_v2) & hypertension2_aha == 0
)
weighting_effect_hyp["psw_cov_iptw_osw_unadj_svy_hyp"] <- summary(fit_psw_cov_iptw_osw_unadj_svy_hyp)$coefficients["insomnia", "Estimate"]
weighting_se_hyp["psw_cov_iptw_osw_unadj_svy_hyp"] <- summary(fit_psw_cov_iptw_osw_unadj_svy_hyp)$coefficients["insomnia", "Std. Error"]
weighting_ci_hyp["psw_cov_iptw_osw_unadj_svy_hyp", ] <- confint(fit_psw_cov_iptw_osw_unadj_svy_hyp)["insomnia",]
weighting_nobs_hyp["psw_cov_iptw_osw_unadj_svy_hyp"] <- nobs(fit_psw_cov_iptw_osw_unadj_svy_hyp)

# unadjusted survey regression with PSW using weighting by the odds
fit_psw_cov_odds_unadj_svy_hyp <-  svyglm(
  formula_unadj, 
  family = poisson(), 
  design = svy_psw_cov_odds,
  subset = !is.na(hypertension2_aha_v2) & hypertension2_aha == 0
)
weighting_effect_hyp["psw_cov_odds_unadj_svy_hyp"] <- summary(fit_psw_cov_odds_unadj_svy_hyp)$coefficients["insomnia", "Estimate"]
weighting_se_hyp["psw_cov_odds_unadj_svy_hyp"] <- summary(fit_psw_cov_odds_unadj_svy_hyp)$coefficients["insomnia", "Std. Error"]
weighting_ci_hyp["psw_cov_odds_unadj_svy_hyp", ] <- confint(fit_psw_cov_odds_unadj_svy_hyp)["insomnia",]
weighting_nobs_hyp["psw_cov_odds_unadj_svy_hyp"] <- nobs(fit_psw_cov_odds_unadj_svy_hyp)

# unadjusted survey regression with PSW x OSW using weighting by the odds
fit_psw_cov_odds_osw_unadj_svy_hyp <-  svyglm(
  formula_unadj, 
  family = poisson(), 
  design = svy_psw_cov_osw_odds,
  subset = !is.na(hypertension2_aha_v2) & hypertension2_aha == 0
)
weighting_effect_hyp["psw_cov_odds_osw_unadj_svy_hyp"] <- summary(fit_psw_cov_odds_osw_unadj_svy_hyp)$coefficients["insomnia", "Estimate"]
weighting_se_hyp["psw_cov_odds_osw_unadj_svy_hyp"] <- summary(fit_psw_cov_odds_osw_unadj_svy_hyp)$coefficients["insomnia", "Std. Error"]
weighting_ci_hyp["psw_cov_odds_osw_unadj_svy_hyp", ] <- confint(fit_psw_cov_odds_osw_unadj_svy_hyp)["insomnia",]
weighting_nobs_hyp["psw_cov_odds_osw_unadj_svy_hyp"] <- nobs(fit_psw_cov_odds_osw_unadj_svy_hyp)

# adjusted survey regression with PSW using IPTW
fit_psw_cov_iptw_adj_svy_hyp <- svyglm(
  formula_adj, 
  family = poisson(), 
  design = svy_psw_cov_iptw,
  subset = !is.na(hypertension2_aha_v2) & hypertension2_aha == 0
)
weighting_effect_hyp["psw_cov_iptw_adj_svy_hyp"] <- summary(fit_psw_cov_iptw_adj_svy_hyp)$coefficients["insomnia", "Estimate"]
weighting_se_hyp["psw_cov_iptw_adj_svy_hyp"] <- summary(fit_psw_cov_iptw_adj_svy_hyp)$coefficients["insomnia", "Std. Error"]
weighting_ci_hyp["psw_cov_iptw_adj_svy_hyp", ] <- confint(fit_psw_cov_iptw_adj_svy_hyp)["insomnia",]
weighting_nobs_hyp["psw_cov_iptw_adj_svy_hyp"] <- nobs(fit_psw_cov_iptw_adj_svy_hyp)

# adjusted survey regression with PSW x OSW using IPTW
fit_psw_cov_iptw_osw_adj_svy_hyp <- svyglm(
  formula_adj, 
  family = poisson(), 
  design = svy_psw_cov_osw_iptw,
  subset = !is.na(hypertension2_aha_v2) & hypertension2_aha == 0
)
weighting_effect_hyp["psw_cov_iptw_osw_adj_svy_hyp"] <- summary(fit_psw_cov_iptw_osw_adj_svy_hyp)$coefficients["insomnia", "Estimate"]
weighting_se_hyp["psw_cov_iptw_osw_adj_svy_hyp"] <- summary(fit_psw_cov_iptw_osw_adj_svy_hyp)$coefficients["insomnia", "Std. Error"]
weighting_ci_hyp["psw_cov_iptw_osw_adj_svy_hyp", ] <- confint(fit_psw_cov_iptw_osw_adj_svy_hyp)["insomnia",]
weighting_nobs_hyp["psw_cov_iptw_osw_adj_svy_hyp"] <- nobs(fit_psw_cov_iptw_osw_adj_svy_hyp)

# adjusted survey regression with PSW using weighting by the odds
fit_psw_cov_odds_adj_svy_hyp <-  svyglm(
  formula_adj, 
  family = poisson(), 
  design = svy_psw_cov_odds,
  subset = !is.na(hypertension2_aha_v2) & hypertension2_aha == 0
)
weighting_effect_hyp["psw_cov_odds_adj_svy_hyp"] <- summary(fit_psw_cov_odds_adj_svy_hyp)$coefficients["insomnia", "Estimate"]
weighting_se_hyp["psw_cov_odds_adj_svy_hyp"] <- summary(fit_psw_cov_odds_adj_svy_hyp)$coefficients["insomnia", "Std. Error"]
weighting_ci_hyp["psw_cov_odds_adj_svy_hyp", ] <- confint(fit_psw_cov_odds_adj_svy_hyp)["insomnia",]
weighting_nobs_hyp["psw_cov_odds_adj_svy_hyp"] <- nobs(fit_psw_cov_odds_adj_svy_hyp)

# adjusted survey regression with PSW x OSW using weighting by the odds
fit_psw_cov_odds_osw_adj_svy_hyp <-  svyglm(
  formula_adj, 
  family = poisson(), 
  design = svy_psw_cov_osw_odds,
  subset = !is.na(hypertension2_aha_v2) & hypertension2_aha == 0
)
weighting_effect_hyp["psw_cov_odds_osw_adj_svy_hyp"] <- summary(fit_psw_cov_odds_osw_adj_svy_hyp)$coefficients["insomnia", "Estimate"]
weighting_se_hyp["psw_cov_odds_osw_adj_svy_hyp"] <- summary(fit_psw_cov_odds_osw_adj_svy_hyp)$coefficients["insomnia", "Std. Error"]
weighting_ci_hyp["psw_cov_odds_osw_adj_svy_hyp", ] <- confint(fit_psw_cov_odds_osw_adj_svy_hyp)["insomnia",]
weighting_nobs_hyp["psw_cov_odds_osw_adj_svy_hyp"] <- nobs(fit_psw_cov_odds_osw_adj_svy_hyp)






