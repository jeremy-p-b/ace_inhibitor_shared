# import globals
prepare_clusters <- function() {
  source("scripts/load_globals.R")
  source("~/glp1/R/analysis_helpers_weighting.R")
  source("~/glp1/R/analysis_helpers.R")
  blas_set_num_threads(1)
}
prepare_clusters()

# load data -----------------------------------------------------------------------------------
YEAR <- 2025
source(here::here("scripts", "load_pregnancies.R"))
aceis <- read_parquet(glue("{COMP_DATA}/ace_inhibitor{YEAR}.parquet"))

# Define covariates -----------------------------------------------------------------------------------

source(here::here("scripts", "define_covariates.R"))

# Prepare data -----------------------------------------------------------------------------------

acei_pregnancies <- pregnancies %>% filter(prior_ace_inhibitor==1)
acei_pregnancies %<>% filter(last_ace_inhibitor+ddays(last_ace_inhibitor_daysupp) >= lmp)
acei_pregnancies %<>% filter(last_ace_inhibitor_daysupp %in% c(30,90))
acei_pregnancies %<>% filter(emprel!=3)

# Create intervention probabilities -----------------------------------------------------------

# inspect distribution of time between medications in year prior
dispensing_patterns <- calculate_dispensing_patterns(acei_pregnancies, aceis, lmp + ddays(13), lmp - ddays(365-14))
chosen_pints <- calculate_pints(dispensing_patterns, c(30,90))

# compare to distribution during pregnancy
clone_and_calculate_weights_rxstrategy(acei_pregnancies, aceis, chosen_pints, overall_rhs_formula, weeks_extra=4, lmp_offset=14) %>%
  filter(sim_continuer == 1 & drug == 1 & num_dispensed==1 & nonlivebirth==0 & week <= 13) %>% tabyl(weeks_since_last, last_weeksupp) %>%
  adorn_percentages(denominator="col") %>% adorn_pct_formatting() %>% adorn_ns(position="front")

plotted_pints <- chosen_pints %>%
  group_by(last_weeksupp) %>%
  mutate(cumpint = if_else(weeks_since_last==15, sum(as.integer(weeks_since_last >= 15)*cumpint), cumpint)) %>%
  filter(weeks_since_last <= 15)

CairoPDF(glue("output/acei_weeks_until_next{YEAR}.pdf"), width=9, height=6)
ggplot(data=chosen_pints %>% filter(weeks_since_last < 20) %>% mutate(last_weeksupp = factor(last_weeksupp)),
       aes(x=weeks_since_last)) + geom_col(aes(y=cumpint, fill=last_weeksupp), position="identity", alpha=0.7) +
  theme_minimal() +
  labs(fill="Weeks' supply", x="Weeks until next dispensation", y="Proportion of dispensations") +
  theme(
    axis.title = element_text(size=12),
    axis.text = element_text(size=12),
    legend.title = element_text(size=12),
    legend.text = element_text(size=12),
  )
dev.off()

# calculate effective sample size
acei_pregnancies %>%
  clone_and_calculate_weights_rxstrategy(aceis, chosen_pints, overall_rhs_formula, weeks_extra=4, lmp_offset=14) %>%
  filter(week == 42 & sim_continuer == 1 & censor == 0 & nonlivebirth == 0) %>%
  summarise(across(starts_with("ipcw"), ~ sum(.)^2/sum(.^2)))

# Estimates -----------------------------------------------------------------------------------

# examine week of dispensing
clone_and_calculate_weights_rxstrategy(acei_pregnancies, aceis, chosen_pints, overall_rhs_formula, weeks_extra=4, lmp_offset=14) %>%
  filter(sim_continuer == 1 & drug == 1 & num_dispensed==1 & nonlivebirth==0 & censor==0) %>% tabyl(week, last_weeksupp) %>%
  adorn_percentages(denominator="col") %>% adorn_pct_formatting() %>% adorn_ns(position="front")

acei_pregnancies %>% calculate_weighted_estimates_rxstrategy(aceis, chosen_pints, overall_rhs_formula, weeks_extra=4, lmp_offset=14) %>%
  filter(week==42) %>% select(analysis, starts_with("risk_ratio"))

acei_pregnancies %>% calculate_weighted_estimates_rxstrategy(aceis, chosen_pints, overall_rhs_formula, weeks_extra=4, lmp_offset=14, gee_livebirth=TRUE) %>%
  select(analysis, starts_with("risk_ratio")) %>% View()

# Analysis of ACE inhibitors---------------------------------------------------

# main analysis
clust <- makeCluster(5)
clusterCall(clust, prepare_clusters)
acei_adjusted_estimates <- pblapply(1:500, sample_compute_weighted_rxstrategy, pregnancies=acei_pregnancies, drugs=aceis, pints=chosen_pints, rhs_formula=overall_rhs_formula, weeks_extra=4, lmp_offset=14, cl=clust)
acei_adjusted_estimates %<>% summarise_boot_weighted_rxstrategy(pregnancies=acei_pregnancies, drugs=aceis, pints=chosen_pints, rhs_formula=overall_rhs_formula, weeks_extra=4, lmp_offset=14)
acei_adjusted_estimates %>% write_csv(glue("output/acei_adjusted_results{YEAR}.csv"))
stopCluster(clust)

# no risk related behaviour
acei_pregnancies_noriskrelated <- acei_pregnancies %>%
  filter(prior_teratogen==0 & alcohol_abuse==0 & drug_abuse==0 & tobacco_use==0)

clust <- makeCluster(5)
clusterCall(clust, prepare_clusters)
acei_adjusted_estimates_noriskrelated <- pblapply(1:500, sample_compute_weighted_rxstrategy, pregnancies=acei_pregnancies_noriskrelated, drugs=aceis, pints=chosen_pints, rhs_formula=overall_rhs_formula_nolifestyle, weeks_extra=4, lmp_offset=14, cl=clust)
acei_adjusted_estimates_noriskrelated %<>% summarise_boot_weighted_rxstrategy(pregnancies=acei_pregnancies_noriskrelated, drugs=aceis, pints=chosen_pints, rhs_formula=overall_rhs_formula_nolifestyle, weeks_extra=4, lmp_offset=14)
acei_adjusted_estimates_noriskrelated %>% write_csv(glue("output/acei_adjusted_results_noriskrelated{YEAR}.csv"))
stopCluster(clust)

# less than 5 weeks supply remaining at conception
acei_pregnancies_less5weeks <- acei_pregnancies %>%
  filter(last_ace_inhibitor+ddays(last_ace_inhibitor_daysupp) < lmp + ddays(49))

clust <- makeCluster(5)
clusterCall(clust, prepare_clusters)
acei_adjusted_estimates_less5weeks <- pblapply(1:500, sample_compute_weighted_rxstrategy, pregnancies=acei_pregnancies_less5weeks, drugs=aceis, pints=chosen_pints, rhs_formula=overall_rhs_formula, weeks_extra=4, lmp_offset=14, cl=clust)
acei_adjusted_estimates_less5weeks %<>% summarise_boot_weighted_rxstrategy(pregnancies=acei_pregnancies_less5weeks, drugs=aceis, pints=chosen_pints, rhs_formula=overall_rhs_formula, weeks_extra=4, lmp_offset=14)
acei_adjusted_estimates_less5weeks %>% write_csv(glue("output/acei_adjusted_results_less5weeks{YEAR}.csv"))
stopCluster(clust)

# after oct 2015
acei_pregnancies_afteroct2015 <- acei_pregnancies %>%
  filter(lmp >= ymd("2015-10-01"))

clust <- makeCluster(5)
clusterCall(clust, prepare_clusters)
acei_adjusted_estimates_afteroct2015 <- pblapply(1:500, sample_compute_weighted_rxstrategy, pregnancies=acei_pregnancies_afteroct2015, drugs=aceis, pints=chosen_pints, rhs_formula=overall_rhs_formula, weeks_extra=4, lmp_offset=14, cl=clust)
acei_adjusted_estimates_afteroct2015 %<>% summarise_boot_weighted_rxstrategy(pregnancies=acei_pregnancies_afteroct2015, drugs=aceis, pints=chosen_pints, rhs_formula=overall_rhs_formula, weeks_extra=4, lmp_offset=14)
acei_adjusted_estimates_afteroct2015 %>% write_csv(glue("output/acei_adjusted_results_afteroct2015{YEAR}.csv"))
stopCluster(clust)

# 90 day supply
acei_pregnancies_30daysupply <- acei_pregnancies %>%
  filter(last_ace_inhibitor_daysupp==30)

clust <- makeCluster(5)
clusterCall(clust, prepare_clusters)
acei_adjusted_estimates_30daysupply <- pblapply(1:500, sample_compute_weighted_rxstrategy, pregnancies=acei_pregnancies_30daysupply, drugs=aceis, pints=chosen_pints, rhs_formula=overall_rhs_formula, weeks_extra=4, lmp_offset=14, cl=clust)
acei_adjusted_estimates_30daysupply %<>% summarise_boot_weighted_rxstrategy(pregnancies=acei_pregnancies_30daysupply, drugs=aceis, pints=chosen_pints, rhs_formula=overall_rhs_formula, weeks_extra=4, lmp_offset=14)
acei_adjusted_estimates_30daysupply %>% write_csv(glue("output/acei_adjusted_results_30daysupply{YEAR}.csv"))
stopCluster(clust)

# 90 day supply
acei_pregnancies_90daysupply <- acei_pregnancies %>%
  filter(last_ace_inhibitor_daysupp==90)

clust <- makeCluster(5)
clusterCall(clust, prepare_clusters)
acei_adjusted_estimates_90daysupply <- pblapply(1:500, sample_compute_weighted_rxstrategy, pregnancies=acei_pregnancies_90daysupply, drugs=aceis, pints=chosen_pints, rhs_formula=overall_rhs_formula, weeks_extra=4, lmp_offset=14, cl=clust)
acei_adjusted_estimates_90daysupply %<>% summarise_boot_weighted_rxstrategy(pregnancies=acei_pregnancies_90daysupply, drugs=aceis, pints=chosen_pints, rhs_formula=overall_rhs_formula, weeks_extra=4, lmp_offset=14)
acei_adjusted_estimates_90daysupply %>% write_csv(glue("output/acei_adjusted_results_90daysupply{YEAR}.csv"))
stopCluster(clust)


