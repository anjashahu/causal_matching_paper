
### SET UP ----------------------------------------------------------------------

# load packages
library(tidyverse)
library(survey)

# clear working environment
rm(list = ls())

### LOAD IN THE REAL DATA -------------------------------------------------------

# load data
hchs <- read.csv("./data/hchs_paper_data/solshahu_covariates_20210810.csv")
# select needed variables and create insomnia variable
hchs <- hchs %>% 
  select(
    STRAT, 
    PSU_ID, 
    WEIGHT_FINAL_NORM_OVERALL, 
    WEIGHT_NORM_OVERALL_INCA,
    AGE, 
    BMI, 
    MCI,
    HYPERTENSION2_AHA, 
    HYPERTENSION2_AHA_V2,
    YRS_BTWN_V1V2,
    WHIIRS
  ) %>%
  # make all column names lower case
  rename_all(.funs = tolower) %>%
  # create insomnia variable
  mutate(insomnia = ifelse(whiirs >= 9, 1, 0))

### OBTAIN COEFFICIENTS FOR INSOMNIA ----------------------------------------------

# for insomnia
fit_svy_insomnia <- svyglm(
  insomnia ~ bmi + age, 
  family = quasibinomial(), 
  design = svydesign(
    ids = ~ psu_id, 
    weights = ~ weight_final_norm_overall, 
    strata = ~ strat, 
    data = hchs
  )
)
summary(fit_svy_insomnia)

### OBTAIN COEFFICIENTS FOR OUTCOMES ----------------------------------------------

# when generating the outcome in the simulated population, we will use
# coefficient values that are based on the real HCHS/SOL data
# so we are obtaining those values here

# for mci (at visit 2)
fit_svy_mci <- svyglm(
  mci ~ insomnia + bmi + age, 
  family = quasibinomial(), 
  design = svydesign(
    ids = ~ psu_id, 
    weights = ~ weight_norm_overall_inca, 
    strata = ~ strat, 
    data = hchs %>% filter(!is.na(weight_norm_overall_inca))
  )
)
summary(fit_svy_mci)

# for hypertension (at visit 1)
fit_svy_hyp_v1 <- svyglm(
  hypertension2_aha ~ insomnia + bmi + age, 
  family = quasibinomial(), 
  design = svydesign(
    ids = ~ psu_id, 
    weights = ~ weight_final_norm_overall, 
    strata = ~ strat, 
    data = hchs
  )
)
summary(fit_svy_hyp_v1)

# for hypertension (at visit 2)
fit_svy_hyp_v2 <- svyglm(
  hypertension2_aha_v2 ~ insomnia + bmi + age + yrs_btwn_v1v2, 
  family = quasibinomial(), 
  design = svydesign(
    ids = ~ psu_id, 
    weights = ~ weight_final_norm_overall, 
    strata = ~ strat, 
    data = hchs
  )
)
summary(fit_svy_hyp_v2)

### SAVE RESULTS -----------------------------------------------------------------

# combine results into one data frame
hchs_coef <- data.frame(
  hyp_v2 = coef(fit_svy_hyp_v2),
  hyp_v1 = c(coef(fit_svy_hyp_v1), NA),
  mci = c(coef(fit_svy_mci), NA),
  insomnia = c(coef(fit_svy_insomnia)["(Intercept)"], NA, coef(fit_svy_insomnia)[c("bmi", "age")], NA)
) 
# calculate the odds ratios and round to 3 digits
hchs_or <- round(exp(hchs_coef), digits = 3)
# fix row names of data frame with results
row.names(hchs_or) <- c("intercept", "insomnia", "bmi", "age", "visit_years")
# save outcome results
save(hchs_or, file = "./data/simulated_pop_data/hchs_or.RData")
