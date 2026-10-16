# import globals
source("scripts/load_globals.R")
source("~/glp1/R/table_helpers.R")
source("~/glp1/R/plot_km_survival_ci.R")
source("~/glp1/R/plot_km_survival_ci_by_analysis.R")
source("~/glp1/R/analysis_helpers_weighting.R")
source("~/glp1/R/analysis_helpers.R")

YEAR <- 2025

source("scripts/load_pregnancies.R")
source(here::here("scripts", "define_covariates.R"))
aceis <- read_parquet(glue("{COMP_DATA}/ace_inhibitor{YEAR}.parquet"))

calculate_denoms <- function(pregnancies, expanded=TRUE) {
  if (expanded) {
    pregnancies %<>% group_by(newid, sim_continuer) %>%
      filter(week == max(week) & censor == 0) %>%
      group_by(pregid) %>%
      filter(row_number()==1) %>%
      ungroup()
  }
  linked_livebirths <- pregnancies %>% filter(nonlivebirth == 0 & infant_enrol30to90==1)

  denoms <- tribble(~Outcome, ~`No. outcomes`,
                    "Non-live birth",  pregnancies %>% filter(!(final_type %in% c(1,2))) %>% nrow(),
                    "Spontaneous abortion", pregnancies %>%  filter(final_type == 3) %>% nrow(),
                    "Stillbirth", pregnancies %>% filter(final_type == 5) %>% nrow(),
                    "Elective termination", pregnancies %>% filter(final_type ==4) %>% nrow(),
                    "Linked live births", linked_livebirths %>% nrow(),
                    "Live birth with SGA", linked_livebirths %>% filter(sga==1) %>% nrow(),
                    "SGA among live births", linked_livebirths %>% filter(sga==1) %>% nrow(),
                    "Live birth with LGA", linked_livebirths %>% filter(lga==1) %>% nrow(),
                    "LGA among live births", linked_livebirths %>% filter(lga==1) %>% nrow(),
                    "Live birth with MCM", linked_livebirths %>% filter(mcm_no_chrom==1) %>% nrow(),
                    "MCM among live births", linked_livebirths %>% filter(mcm_no_chrom==1) %>% nrow())
  return(denoms)
}

# add labels to malformations
variable_lbls <- read_excel("~/glp1/variable_labels.xlsx")
for (varname in variable_lbls$varname) {
  if (paste0(varname, "_no_chrom") %in% colnames(pregnancies)) {
    var_label(pregnancies[[paste0(varname, "_no_chrom")]]) <- capitalize(variable_lbls[variable_lbls[["varname"]] == varname, "label"][[1,1]])
  }
}

acei_pregnancies <- pregnancies %>% filter(prior_ace_inhibitor==1)
acei_pregnancies %<>% filter(last_ace_inhibitor+ddays(last_ace_inhibitor_daysupp) >= lmp)
acei_pregnancies %<>% filter(last_ace_inhibitor_daysupp %in% c(30,90))
acei_pregnancies %<>% filter(emprel!=3)

dispensing_patterns <- calculate_dispensing_patterns(acei_pregnancies, aceis, lmp + ddays(13), lmp - ddays(365-14))
chosen_pints <- calculate_pints(dispensing_patterns, c(30,90))

acei_cloned_expanded_pregnancies <- acei_pregnancies %>%
  clone_and_calculate_weights_rxstrategy(aceis, chosen_pints, overall_rhs_formula, weeks_extra=4, lmp_offset=14)

acei_cloned_expanded_pregnancies_noriskrelated <- acei_pregnancies %>%
  clone_and_calculate_weights_rxstrategy(aceis, chosen_pints, overall_rhs_formula, weeks_extra=4, lmp_offset=14)

acei_cloned_expanded_pregnancies_less5weeks <- acei_pregnancies %>%
  filter(last_ace_inhibitor+ddays(last_ace_inhibitor_daysupp) < lmp + ddays(49)) %>%
  clone_and_calculate_weights_rxstrategy(aceis, chosen_pints, overall_rhs_formula, weeks_extra=4, lmp_offset=14)

acei_cloned_expanded_pregnancies_afteroct2015 <- acei_pregnancies %>%
  filter(lmp >= ymd("2015-10-01")) %>%
  clone_and_calculate_weights_rxstrategy(aceis, chosen_pints, overall_rhs_formula, weeks_extra=4, lmp_offset=14)

# Primary results -----------------------------------------------------------------------------

acei_overall_adjusted_results_weighted <- read_csv(glue("output/acei_adjusted_results{YEAR}.csv"))
risks_overall_adjusted <- create_table_by_analysis(acei_overall_adjusted_results_weighted)

out_nos <- calculate_denoms(acei_cloned_expanded_pregnancies)

outcome_list <- c("Non-live birth", "Spontaneous abortion", "Elective termination",
                  "MCM among live births", "Live birth with MCM")

forest_rr_crude <- risks_overall_adjusted %>%
  filter(Analysis == "Unweighted - nonparametric outcome model") %>%
  create_risk_forest_plot(out_nos, outcome_list)

forest_rr_adjusted <- risks_overall_adjusted %>%
  filter(Analysis == "IPCW user specified - nonparametric outcome model") %>%
  create_risk_forest_plot(out_nos, outcome_list)


CairoPDF(glue("output/acei_forest_plot_crude{YEAR}.pdf"), width=10, height=5)
forest_rr_crude
dev.off()

CairoPDF(glue("output/acei_forest_plot_adjusted{YEAR}.pdf"), width=10, height=5)
forest_rr_adjusted %>% fp_set_style(box="black", line="black", txt_gp=fpTxtGp(summary = gpar(fontfamily = "", cex = 1),
                                                                              label = gpar(fontfamily = "", cex = 1),ticks = gpar(fontfamily = "", cex = 1)))
dev.off()



# SE - no risk-related  -----------------------------------------------------------------------------

acei_overall_adjusted_results_weighted_noriskrelated <- read_csv(glue("output/acei_adjusted_results_noriskrelated{YEAR}.csv"))
risks_overall_adjusted_noriskrelated <- create_table_by_analysis(acei_overall_adjusted_results_weighted_noriskrelated)

out_nos_noriskrelated <- calculate_denoms(acei_cloned_expanded_pregnancies %>% filter(prior_teratogen==0 & alcohol_abuse==0 & drug_abuse==0 & tobacco_use==0))

outcome_list <- c("Non-live birth", "Spontaneous abortion", "Elective termination",
                  "MCM among live births", "Live birth with MCM")

forest_rr_adjusted_noriskrelated <- risks_overall_adjusted_noriskrelated %>%
  filter(Analysis == "IPCW user specified - nonparametric outcome model") %>%
  create_risk_forest_plot(out_nos_noriskrelated, outcome_list)

CairoPDF(glue("output/acei_forest_plot_adjusted_noriskrelated{YEAR}.pdf"), width=10, height=5)
forest_rr_adjusted_noriskrelated
dev.off()

# SE - less than 5 weeks -----------------------------------------------------------------------------

acei_overall_adjusted_results_weighted_less5weeks <- read_csv(glue("output/acei_adjusted_results_less5weeks{YEAR}.csv"))
risks_overall_adjusted_less5weeks <- create_table_by_analysis(acei_overall_adjusted_results_weighted_less5weeks)

out_nos_less5weeks <- calculate_denoms(acei_cloned_expanded_pregnancies_less5weeks)

outcome_list <- c("Non-live birth", "Spontaneous abortion", "Elective termination",
                  "MCM among live births", "Live birth with MCM")

forest_rr_adjusted_less5weeks <- risks_overall_adjusted_less5weeks %>%
  filter(Analysis == "IPCW user specified - nonparametric outcome model") %>%
  create_risk_forest_plot(out_nos_less5weeks, outcome_list)

CairoPDF(glue("output/acei_forest_plot_adjusted_less5weeks{YEAR}.pdf"), width=10, height=5)
forest_rr_adjusted_less5weeks
dev.off()

# SE - after 2015 -----------------------------------------------------------------------------

acei_overall_adjusted_results_weighted_afteroct2015 <- read_csv(glue("output/acei_adjusted_results_afteroct2015{YEAR}.csv"))
risks_overall_adjusted_afteroct2015 <- create_table_by_analysis(acei_overall_adjusted_results_weighted_afteroct2015)

out_nos_afteroct2015 <- calculate_denoms(acei_cloned_expanded_pregnancies_afteroct2015)

outcome_list <- c("Non-live birth", "Spontaneous abortion", "Elective termination",
                  "MCM among live births", "Live birth with MCM")

forest_rr_adjusted_afteroct2015 <- risks_overall_adjusted_afteroct2015 %>%
  filter(Analysis == "IPCW user specified - nonparametric outcome model") %>%
  create_risk_forest_plot(out_nos_afteroct2015, outcome_list)

CairoPDF(glue("output/acei_forest_plot_adjusted_afteroct2015{YEAR}.pdf"), width=10, height=5)
forest_rr_adjusted_afteroct2015
dev.off()

# Denominators -----------------------------------------------------------------------

denoms <- tribble(
  ~ Group, ~Subgroup, ~No.,
  "Overall", "Overall", acei_pregnancies %>% nrow(),
  "Overall", "Uncensored linked live births", acei_cloned_expanded_pregnancies %>% filter(week==42 & infant_enrol30to90==1 & censor==0) %>% nrow(),
  "No risk-related", "No risk-related factors", acei_pregnancies %>%
    filter(prior_teratogen==0 & alcohol_abuse==0 & drug_abuse==0 & tobacco_use==0) %>%
    nrow(),
  "No risk-related", "Uncensored linked live births", acei_cloned_expanded_pregnancies %>%
    filter(prior_teratogen==0 & alcohol_abuse==0 & drug_abuse==0 & tobacco_use==0) %>% filter(week==42 & infant_enrol30to90==1 & censor==0) %>% nrow(),
  "Less 5 weeks", "No risk-related factors", acei_pregnancies %>%
    filter(last_ace_inhibitor+ddays(last_ace_inhibitor_daysupp) < lmp + ddays(49))  %>% nrow(),
  "Less 5 weeks", "Uncensored linked live births", acei_cloned_expanded_pregnancies_less5weeks %>% filter(week==42 & infant_enrol30to90==1 & censor==0) %>% nrow(),
  "After Oct 2015", "No risk-related factors", acei_pregnancies %>%
    filter(lmp >= ymd("2015-10-01")) %>% nrow(),
  "After Oct 2015", "Uncensored linked live births", acei_cloned_expanded_pregnancies_afteroct2015 %>% filter(week==42 & infant_enrol30to90==1 & censor==0) %>% nrow())
denoms %>% write_csv(glue("output/acei_denoms{YEAR}.csv"))

# Malformations -----------------------------------------------------------------------

malformations_list <- c("mcm_no_chrom", "cardio_no_chrom", "oral_no_chrom", "cns_no_chrom", "eye_no_chrom", "ear_no_chrom",
                        "vas_no_chrom", "resp_no_chrom", "gastro_no_chrom", "genit_overall_no_chrom", "urinary_no_chrom",
                        "muscu_no_chrom", "limb_no_chrom", "abdom_no_chrom", "other_no_chrom")

pretty_malformations_list <- c("mcm_no_chrom" = "Any major", "cardio_no_chrom" = "Cardiovascular", "oral_no_chrom" = "Oral", "cns_no_chrom" = "CNS", "eye_no_chrom" = "Eye", "ear_no_chrom" ="Ear",
                               "vas_no_chrom" = "Other vascular", "resp_no_chrom" = "Respiratory", "gastro_no_chrom" = "Gastrointestinal", "genit_overall_no_chrom"="Genital", "urinary_no_chrom" = "Urinary",
                               "muscu_no_chrom" = "Musculoskeletal", "limb_no_chrom"="Limb", "abdom_no_chrom"="Abdominal wall", "other_no_chrom"="Other")

acei_uncensored_linked_livebirths <- acei_cloned_expanded_pregnancies %>% filter(infant_enrol30to90==1 & censor==0 & week==42)

# crude
res <- tibble_row(Group = "Denominator", `ACEi continuer`=acei_uncensored_linked_livebirths %>% filter(sim_continuer==1) %>% nrow(),
                  `ACEi discontinuer`=acei_uncensored_linked_livebirths %>% filter(sim_continuer==0) %>% nrow())

for (malformation in malformations_list) {
  res %<>% bind_rows(
    tibble_row(Group=pretty_malformations_list[[malformation]],
               `ACEi continuer`= sum(acei_uncensored_linked_livebirths[acei_uncensored_linked_livebirths["sim_continuer"]==1, malformation]),
               `ACEi discontinuer`= sum(acei_uncensored_linked_livebirths[acei_uncensored_linked_livebirths["sim_continuer"]==0, malformation])
    ))
}

res %<>% mutate(`ACEi continuer` = format_perc(`ACEi continuer`, acei_uncensored_linked_livebirths %>% filter(sim_continuer==1) %>% nrow()))
res %<>% mutate(`ACEi discontinuer` = format_perc(`ACEi discontinuer`, acei_uncensored_linked_livebirths %>% filter(sim_continuer==0) %>% nrow()))

res %>% write_csv(glue("output/acei_malformation_counts{YEAR}.csv"))

# weighted
res <- tibble_row(Group = "Denominator", `ACEi continuer`=acei_uncensored_linked_livebirths %>% filter(sim_continuer==1) %>% {sum(.[["ipcw_user"]])},
                  `ACEi discontinuer`=acei_uncensored_linked_livebirths %>% filter(sim_continuer==0) %>% {sum(.[["ipcw_user"]])})

for (malformation in malformations_list) {
  res %<>% bind_rows(
    tibble_row(Group=pretty_malformations_list[[malformation]],
               `ACEi continuer`= sum(acei_uncensored_linked_livebirths[acei_uncensored_linked_livebirths["sim_continuer"]==1, malformation] * acei_uncensored_linked_livebirths[acei_uncensored_linked_livebirths["sim_continuer"]==1, "ipcw_user"]),
               `ACEi discontinuer`= sum(acei_uncensored_linked_livebirths[acei_uncensored_linked_livebirths["sim_continuer"]==0, malformation] * acei_uncensored_linked_livebirths[acei_uncensored_linked_livebirths["sim_continuer"]==0, "ipcw_user"])
    ))
}
res %<>% mutate(across(c(`ACEi continuer`, `ACEi discontinuer`), round))

res %<>% mutate(`ACEi continuer` = format_perc(`ACEi continuer`, acei_uncensored_linked_livebirths %>% filter(sim_continuer==1) %>% {sum(.[["ipcw_user"]])}))
res %<>% mutate(`ACEi discontinuer` = format_perc(`ACEi discontinuer`, acei_uncensored_linked_livebirths %>% filter(sim_continuer==0) %>% {sum(.[["ipcw_user"]])}))

res %>% write_csv(glue("output/acei_malformation_counts_weighted{YEAR}.csv"))
