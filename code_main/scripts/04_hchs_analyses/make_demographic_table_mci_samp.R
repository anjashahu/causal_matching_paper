

### SET UP ---------------------------------------------------------------------

# load packages
library(tidyverse)
library(survey)
library(MatchIt)
library(tableone)

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
# and remove observations with missing values from predictors
mci <- hchs %>%
  select(-c(weight_final_norm_overall, hypertension2_aha, hypertension2_aha_v2, mci3, yrs_btwn_v1v2)) %>%
  drop_na(insomnia, age, bmi, gender, marital_status, employed, alcohol_use, cigarette_use, bkgrd1_c7, education_c3)
apply(mci, 2, function(col) sum(is.na(col)))

# remove observations with missing values from outcome
mci <- mci %>%
  drop_na(mci, weight_norm_overall_inca) 
apply(mci, 2, function(col) sum(is.na(col)))

###  TABLE ----------------------------------------------------------------------

# get unweighted totals first
mci %>%
  count(insomnia) %>%
  bind_rows(colSums(.[,-1]))

# recode categorical variables 
# background
mci$bkgrd1_c7_cat <- 
  factor(
    mci$bkgrd1_c7,
    levels = 0:6,
    labels = c("Dominican", "Central American", "Cuban", 
    "Mexican", "Puerto Rican", "South American", "More than one/Other heritage")
  )

# alcohol
mci$alcohol_use_cat <- 
  factor(
    mci$alcohol_use,
    levels = 1:3,
    labels = c("Never", "Former", "Current")
  )

# smoking
mci$cigarette_use_cat <- 
  factor(
    mci$cigarette_use,
    levels = 1:3,
    labels = c("Never", "Former", "Current")
  )

# marital status
mci$marital_status_cat <- 
  factor(
    mci$marital_status,
    levels = 1:3,
    labels = c("Single", "Married or living with partner", "Separated, divorced, or widow(er)")
  )

# education
mci$education_c3_cat <- 
  factor(
    mci$education_c3,
    levels = 1:3,
    labels = c("No high school diploma or GED", "At most a High school diploma/GED", "Greater than high school (or GED) education")
  )

# employment
mci$employed_cat <- 
  factor(
    mci$employed,
    levels = 1:4,
    labels = c("Retired and not currently employed", "Not retired", "Employed part-time (<35 hours/week)", "Employed full-time (>35 hours/week)")
  )


# create survey design object of data
svy_mci <- svydesign(
  ids = ~ psu_id, 
  weights = ~ weight_norm_overall_inca, 
  strata = ~ strat, 
  data = mci
)

# get weighted table
table_mci_svy <- svyCreateTableOne(
  vars =  c("bkgrd1_c7_cat", "alcohol_use_cat", "cigarette_use_cat", "age", "gender",
            "marital_status_cat", "education_c3_cat", "bmi", "employed_cat"), 
  factorVars = c("bkgrd1_c7_cat", "alcohol_use_cat", "cigarette_use_cat", 
                 "gender", "marital_status_cat", "education_c3_cat", "employed_cat"),
  strata = "insomnia", addOverall = TRUE, test = FALSE,
  data = svy_mci
)
print(table_mci_svy, varLabels = TRUE, contDigits = 2, catDigits = 1, showAllLevels = TRUE)


