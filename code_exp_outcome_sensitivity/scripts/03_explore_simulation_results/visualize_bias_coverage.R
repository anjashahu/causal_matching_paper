
### SET UP ---------------------------------------------------------------------

# load packages
library(tidyverse)

# clear working environment
rm(list = ls())

# load  results 
df_mci_main <- readRDS("./data/simulation_results/bias_coverage_results/effects_mci_main.RData") %>%
  mutate(specification = "correct")
df_hyp_main <- readRDS("./data/simulation_results/bias_coverage_results/effects_hyp_main.RData") %>%
  mutate(specification = "correct")
df_mci_confounding <- readRDS("./data/simulation_results/bias_coverage_results/effects_mci_confounding.RData") %>%
  mutate(specification = "under (no age)")
df_hyp_confounding <- readRDS("./data/simulation_results/bias_coverage_results/effects_hyp_confounding.RData") %>%
  mutate(specification = "under (no age)")

# bind results for each outcome
df_mci <- bind_rows(df_mci_main, df_mci_confounding) 
df_hyp <- bind_rows(df_hyp_main, df_hyp_confounding)

### CREATE PLOTS -------------------------------------------------------------


df_hyp_confounding %>% 
  bind_rows(df_hyp_main, df_mci_confounding, df_mci_main) %>%
  pull(coverage) %>% 
  summary()

df_hyp_confounding %>% 
  bind_rows(df_hyp_main, df_mci_confounding, df_mci_main) %>%
  pull(bias) %>% 
  summary()

create_plot <- function(df, comparison) {
  if (comparison == "coverage") {
    plot <- df %>%
      mutate_at(vars(effect, true, bias, coverage), list(~ round(., 3))) %>%
      ggplot(aes(x = as.factor(exposure_prev), y = as.factor(outcome_prev), fill = coverage)) +
      geom_tile(color = "black") +
      facet_grid(method ~ specification) +
      geom_text(aes(label = coverage), color = "black", size = 2) +
      scale_fill_gradient2(low = "#075AFF", mid = "#FFFFCC", high = "#FF0000", midpoint = 0.95, limits = c(0.85, 1)) +
      coord_fixed() +
      xlab("Exposure prevalence") +
      ylab("Outcome prevalence") +
      theme_minimal() +
      theme(panel.grid.major = element_blank(), legend.position = "bottom")
  } else if (comparison == "bias") {
    plot <- df %>%
      mutate_at(vars(effect, true, bias, coverage), list(~ round(., 3))) %>%
      ggplot(aes(x = as.factor(exposure_prev), y = as.factor(outcome_prev), fill = bias)) +
      geom_tile(color = "black") +
      facet_grid(method ~ specification) +
      geom_text(aes(label = bias), color = "black", size = 2) +
      scale_fill_gradient2(low = "#075AFF", mid = "#FFFFCC", high = "#FF0000", midpoint = 0, limits = c(-0.05, 0.15)) +
      coord_fixed() +
      xlab("Exposure prevalence") +
      ylab("Outcome prevalence") +
      theme_minimal() +
      theme(panel.grid.major = element_blank(), legend.position = "bottom")
  }
  
  return(plot)
}


# look at plots
# mci analyses
create_plot(df_mci, "coverage") + 
  ggtitle("Coverage, MCI")
ggsave("cov_mci", path = "./scripts/03_explore_simulation_results/", device = "jpeg")
create_plot(df_mci, "bias") +
  ggtitle("Bias, MCI")
ggsave("bias_mci", path = "./scripts/03_explore_simulation_results/", device = "jpeg")
# hypertension analyses
create_plot(df_hyp, "coverage") +
  ggtitle("Coverage, hypertension")
ggsave("cov_hyp", path = "./scripts/03_explore_simulation_results/", device = "jpeg")
create_plot(df_hyp, "bias") +
  ggtitle("Bias, hypertension")
ggsave("bias_hyp", path = "./scripts/03_explore_simulation_results/", device = "jpeg")







