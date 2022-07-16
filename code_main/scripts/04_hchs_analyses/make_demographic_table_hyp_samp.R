
### SET UP ---------------------------------------------------------------------

# load packages
library(tidyverse)
library(survey)
library(MatchIt)
library(tableone)
library(knitr)

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
# and remove observations with missing values from predictors
hyp <- hchs %>%
  select(-c(weight_norm_overall_inca, mci, mci3)) %>%
  drop_na(insomnia, age, bmi, gender, marital_status, employed, alcohol_use, cigarette_use, bkgrd1_c7, education_c3)
apply(hyp, 2, function(col) sum(is.na(col)))

# remove observations with missing values from outcomee
hyp <- hyp %>%
  drop_na(hypertension2_aha_v2, yrs_btwn_v1v2) 
apply(hyp, 2, function(col) sum(is.na(col)))

# remove observations with hypertension at baseline 
hyp <- hyp %>% 
  filter(hypertension2_aha == 0)

###  TABLE ----------------------------------------------------------------------

# get unweighted totals first
hyp %>%
  count(insomnia) %>%
  bind_rows(colSums(.[,-1]))

# recode categorical variables 
# background
hyp$bkgrd1_c7_cat <- 
  factor(
    hyp$bkgrd1_c7,
    levels = 0:6,
    labels = c("Dominican", "Central American", "Cuban", 
               "Mexican", "Puerto Rican", "South American", "More than one/Other heritage")
  )

# alcohol
hyp$alcohol_use_cat <- 
  factor(
    hyp$alcohol_use,
    levels = 1:3,
    labels = c("Never", "Former", "Current")
  )

# smoking
hyp$cigarette_use_cat <- 
  factor(
    hyp$cigarette_use,
    levels = 1:3,
    labels = c("Never", "Former", "Current")
  )

# marital status
hyp$marital_status_cat <- 
  factor(
    hyp$marital_status,
    levels = 1:3,
    labels = c("Single", "Married or living with partner", "Separated, divorced, or widow(er)")
  )

# education
hyp$education_c3_cat <- 
  factor(
    hyp$education_c3,
    levels = 1:3,
    labels = c("No high school diploma or GED", "At most a High school diploma/GED", "Greater than high school (or GED) education")
  )

# employment
hyp$employed_cat <- 
  factor(
    hyp$employed,
    levels = 1:4,
    labels = c("Retired and not currently employed", "Not retired", "Employed part-time (<35 hours/week)", "Employed full-time (>35 hours/week)")
  )


# create survey design object of data
svy_hyp <- svydesign(
  ids = ~ psu_id, 
  weights = ~ weight_final_norm_overall, 
  strata = ~ strat, 
  data = hyp
)

# get weighted table
table_hyp_svy <- svyCreateTableOne(
  vars =  c("yrs_btwn_v1v2", "bkgrd1_c7_cat", "alcohol_use_cat", "cigarette_use_cat", "age", "gender",
            "marital_status_cat", "education_c3_cat", "bmi", "employed_cat"), 
  factorVars = c("bkgrd1_c7_cat", "alcohol_use_cat", "cigarette_use_cat", 
                 "gender", "marital_status_cat", "education_c3_cat", "employed_cat"),
  strata = "insomnia", addOverall = TRUE, test = FALSE,
  data = svy_hyp
)
print(table_hyp_svy, varLabels = TRUE, contDigits = 2, catDigits = 1, showAllLevels = TRUE)




