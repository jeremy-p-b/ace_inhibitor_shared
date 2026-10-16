# import globals
source("scripts/load_globals.R")
source("~/glp1/R/analysis_helpers.R")
source("~/glp1/R/analysis_helpers_weighting.R")

# load and process data -----------------------------------------------------------------------
YEAR <- 2025
source(here::here("scripts", "load_pregnancies.R"))
aceis <- read_parquet(glue("{COMP_DATA}/ace_inhibitor{YEAR}.parquet"))

# add covariates
pregnancies %<>%
  mutate(overweight_or_obese = as.integer(bmi_cat %in% c("Overweight", "Obese", "Severe obese", "Unspecified obesity")),
         obese = as.integer(bmi_cat %in% c("Obese", "Severe obese", "Unspecified obesity")),
         severe_obese = as.integer(bmi_cat %in% c("Severe obese")),
         last_ace_inhibitor_dosage = if_else(!is.na(last_ace_inhibitor_gennme), paste0(last_ace_inhibitor_gennme, " - ", last_ace_inhibitor_dosage), NA),
         age_lmp_squared = age_lmp^2,
         multiple_linked_infant = as.integer(!is.na(enrolid_inf2) & infant_enrol30to90 ==1),
         num_ace_inhibitor_pregnancy_cat = factor(case_when(num_ace_inhibitor_pregnancy == 0 ~ "0",
                                         num_ace_inhibitor_pregnancy == 1 ~ "1",
                                         num_ace_inhibitor_pregnancy >= 2 ~ "2+",
         ), levels=c("0", "1", "2+")),
         year_lmp = year(lmp),
         REGION = factor(REGION, levels=c(1,2,3,4,5), labels=c("North east", "North central", "South", "West", "Unknown")),
         year_lmp_cat = case_when(year_lmp <= 2016 ~ "2011-2016",
                                  .default = "2017-2023"),
         weeks_since_last_prior=-calculate_week(last_ace_inhibitor, lmp),
         last_prior_ace_inhibitor_weeksupp=floor(as.numeric(last_ace_inhibitor_daysupp/7))
  )


# load variable labels
variable_lbls <- read_excel("~/glp1/variable_labels.xlsx")
for (varname in variable_lbls$varname) {
  if (varname %in% colnames(pregnancies)) {
    var_label(pregnancies[[varname]]) <- capitalize(variable_lbls[variable_lbls[["varname"]] == varname, "label"][[1,1]])
  }
  if (paste0("prior_", varname) %in% colnames(pregnancies)) {
    var_label(pregnancies[[paste0("prior_", varname)]]) <- paste0("Prior ", variable_lbls[variable_lbls[["varname"]] == varname, "label"][[1,1]])
  }
  if (paste0("ft_", varname) %in% colnames(pregnancies)) {
    var_label(pregnancies[[paste0("ft_", varname)]]) <- paste0("First trimester ", variable_lbls[variable_lbls[["varname"]] == varname, "label"][[1,1]])
  }
}


source(here::here("scripts", "define_covariates.R"))

# define cohort -------------------------------------------------------------------------------

# select pregnancies of interest
acei_pregnancies <- pregnancies %>% filter(prior_ace_inhibitor==1)
acei_pregnancies %<>% filter(last_ace_inhibitor+ddays(last_ace_inhibitor_daysupp) >= lmp)
acei_pregnancies %<>% filter(last_ace_inhibitor_daysupp %in% c(30,90))
acei_pregnancies %<>% filter(emprel != 3)

dispensing_patterns <- calculate_dispensing_patterns(acei_pregnancies, aceis, lmp + ddays(13), lmp - ddays(365-14))
chosen_pints <- calculate_pints(dispensing_patterns, c(30,90), keepdenom=TRUE)


add_vars <- function(pregnancies) {
  pregnancies %<>%
    group_by(newid, sim_continuer) %>%
    mutate(week_ace_inhibitor = if_else(sim_continuer==1, max(week*as.integer(drug == 1 & censor==0 & prior_nonlivebirth==0)), min(as.integer(week==0)*(-weeks_since_last-1)))) %>%
    mutate(last_ace_inhibitor_dispensation = case_when(week_ace_inhibitor < 0 ~ "Prior to LMP",
                                              week_ace_inhibitor < 11 ~ "First trimester",
                                              week_ace_inhibitor < 25 ~ "Second trimester",
                                              week_ace_inhibitor >= 25 ~ "Third trimester")) %>%
    mutate(last_ace_inhibitor_dispensation = factor(last_ace_inhibitor_dispensation, levels=c("Prior to LMP", "First trimester", "Second trimester", "Third trimester"))) %>%
    ungroup()
  return(pregnancies)
}

acei_cloned_pregnancies <- acei_pregnancies %>%
  clone_and_calculate_weights_rxstrategy(aceis, chosen_pints, overall_rhs_formula, weeks_extra=4, lmp_offset=14) %>%
  mutate(remaining_supply = as.integer(drug == 1| weeks_since_last <= last_weeksupp))

for (chosen_week in c(2:12)) {
acei_cloned_pregnancies %<>% group_by(newid, sim_continuer) %>%
    arrange(newid, sim_continuer, week) %>%
    mutate(weeks_since_last_previous = first(weeks_since_last)+week) %>%
    mutate(remaining_prior_supply = as.integer(weeks_since_last_previous <= last_prior_ace_inhibitor_weeksupp)) %>%
    mutate("remaining_prior_weeks_supply_week{chosen_week}" := max(as.integer(week == (chosen_week-2) & remaining_prior_supply == 1))) %>%
    mutate("remaining_weeks_supply_week{chosen_week}" := max(as.integer(week == (chosen_week-2) & remaining_supply == 1))) %>%
    ungroup()
}

acei_cloned_pregnancies %<>% group_by(newid, sim_continuer) %>%
  mutate(weeks_supply_at_conception = max(as.integer(week == 0)*(last_prior_ace_inhibitor_weeksupp-weeks_since_last_previous))) %>%
  ungroup()

acei_livebirths <- acei_cloned_pregnancies %>%
  add_vars() %>%
  filter(week == 42 & censor==0 & nonlivebirth==0)

acei_week42 <- acei_cloned_pregnancies %>%
  add_vars() %>%
  filter(week == 42 & censor==0)

var_label(acei_livebirths[["last_ace_inhibitor_dispensation"]]) <- capitalize(variable_lbls[variable_lbls[["varname"]] == "last_ace_inhibitor_dispensation", "label"][[1,1]])
var_label(acei_livebirths[["weeks_supply_at_conception"]]) <- capitalize(variable_lbls[variable_lbls[["varname"]] == "weeks_supply_at_conception", "label"][[1,1]])

acei_cloned_pregnancies_w6 <- acei_pregnancies %>%
  filter(last_ace_inhibitor+ddays(last_ace_inhibitor_daysupp) < lmp + ddays(49)) %>%
  clone_and_calculate_weights_rxstrategy(aceis, chosen_pints, overall_rhs_formula, weeks_extra=4, lmp_offset=14) %>%
  mutate(remaining_supply = as.integer(drug == 1| weeks_since_last <= last_weeksupp))

for (chosen_week in c(2:12)) {
  acei_cloned_pregnancies_w6 %<>% group_by(newid, sim_continuer) %>%
    arrange(newid, sim_continuer, week) %>%
    mutate(weeks_since_last_previous = first(weeks_since_last)+week) %>%
    mutate(remaining_prior_supply = as.integer(weeks_since_last_previous <= last_prior_ace_inhibitor_weeksupp)) %>%
    mutate("remaining_prior_weeks_supply_week{chosen_week}" := max(as.integer(week == (chosen_week-2) & remaining_prior_supply == 1))) %>%
    mutate("remaining_weeks_supply_week{chosen_week}" := max(as.integer(week == (chosen_week-2) & remaining_supply == 1))) %>%
    ungroup()
}

acei_livebirths_w6 <- acei_cloned_pregnancies_w6 %>%
  add_vars() %>%
  filter(week == 42 & censor==0 & nonlivebirth==0)


# For publication -----------------------------------------------------------------------------

# create spreadsheet for publication outputs
wb2 <- createWorkbook()
set_wb_options(wb2)

# Flow chart numbers
inclusions <- read_csv(glue("{COMP_DATA}/acei_pregnancy_cohort_2025_080825.csv"))
inclusions %<>% bind_rows(tibble_row(description="Without child/other",
                                     no_pregnancies=nrow(pregnancies %>%
                                                           filter(emprel!=3))))
inclusions %<>% bind_rows(tibble_row(description="With ACE inhibitor dispensation",
                                     no_pregnancies=nrow(pregnancies %>%
                                                           filter(emprel!=3) %>%
                                                           filter(prior_ace_inhibitor==1))))
inclusions %<>% bind_rows(tibble_row(description="With 30 or 90 days' supply",
                                     no_pregnancies=nrow(pregnancies %>%
                                                           filter(emprel!=3) %>%
                                                           filter(prior_ace_inhibitor==1) %>%
                                                           filter(last_ace_inhibitor_daysupp %in% c(30,90)))))
inclusions %<>% bind_rows(tibble_row(description="With last dispensation + days' supply on or after LMP",
                                     no_pregnancies=nrow(pregnancies %>%
                                                           filter(emprel!=3) %>%
                                                           filter(prior_ace_inhibitor==1) %>%
                                                           filter(last_ace_inhibitor_daysupp %in% c(30,90)) %>%
                                                           filter(last_ace_inhibitor+ddays(last_ace_inhibitor_daysupp) >= lmp))))

inclusions %<>% mutate(excluded = lag(no_pregnancies) - no_pregnancies)

addWorksheet(wb2, sheetName = "Flow chart numbers", gridLines = FALSE)
writeDescriptiveDataTableAuto(wb2, sheet = 1, x = inclusions, description="Flow chart numbers")

# Numbers for text
numbers_for_text <- tribble(
  ~Group, ~`No. pregnancies`,
  "ACEI cohort", nrow(acei_pregnancies),
  "Pregnancy loss", nrow(acei_pregnancies %>% filter(nonlivebirth==1)),
  "No pregnancy loss", nrow(acei_pregnancies %>% filter(nonlivebirth==0)),
  "SAB", nrow(acei_pregnancies %>% filter(spontaneous_abortion==1)),
  "Elective termination", nrow(acei_pregnancies %>% filter(termination==1)),
  "Stillbirth", nrow(acei_pregnancies %>% filter(stillbirth==1)),
  "Uncensored live births", nrow(acei_livebirths),
  "Linked uncensored live births", nrow(acei_livebirths %>% filter(linked_inf==1)),
  "Linked uncensored live births w infant enrolment", nrow(acei_livebirths %>% filter(infant_enrol30to90 == 1)),
  "SGA", nrow(acei_livebirths %>% filter(infant_enrol30to90 == 1 & sga==1)),
  "LGA", nrow(acei_livebirths %>% filter(infant_enrol30to90 == 1 & lga==1)),
  "MCM", nrow(acei_livebirths %>% filter(infant_enrol30to90 == 1 & mcm_no_chrom==1)),
)

addWorksheet(wb2, sheetName = "Numbers for text", gridLines = FALSE)
writeDescriptiveDataTableAuto(wb2, sheet = 2, x = numbers_for_text, description="Numbers for text")

# proportion of dispensations in year prior to LMP by week
distribution_treatment <- chosen_pints %>% group_by(last_weeksupp) %>%
  arrange(last_weeksupp, weeks_since_last) %>%
  mutate(weeks_since_last = weeks_since_last + 1) %>%
  mutate(across(c(cumpint, num), ~ if_else(weeks_since_last == 20, sum(.*as.integer(weeks_since_last >= 20)), .))) %>%
  mutate(normalised_cumpint = if_else(weeks_since_last <= last_weeksupp+4, cumpint/sum(as.integer(weeks_since_last <= last_weeksupp+4)*cumpint), 0)) %>%
  mutate(pint = if_else(weeks_since_last <= last_weeksupp+4, normalised_cumpint/(1-lag(cumsum(normalised_cumpint), default=0)), 0)) %>%
  ungroup() %>%
  filter(weeks_since_last <= 20) %>%
  mutate(weeks_since_last = if_else(weeks_since_last==20, "20+", as.character(weeks_since_last))) %>%
  mutate(weeks_since_last = factor(weeks_since_last, levels=c(as.character(c(0:19)), "20+"))) %>%
  pivot_wider(names_from=last_weeksupp, values_from=c(cumpint, num, normalised_cumpint, pint))
distribution_treatment %<>% rename(`Weeks to next dispensation` = weeks_since_last, `Proportion among four weeks' supply`=cumpint_4,
                                  `Proportion among twelve weeks' supply`=cumpint_12,
                                  `Normalised among four weeks' supply`=normalised_cumpint_4,
                                  `Normalised among twelve weeks' supply`=normalised_cumpint_12,
                                  `Intervention probability among four weeks' supply`=pint_4,
                                  `Intervention probability among twelve weeks' supply`=pint_12,
                                  `Number dispensations among four weeks' supply`=num_4,
                                  `Number dispensations among six weeks' supply`=num_12)

addWorksheet(wb2, sheetName = "Rx distribution", gridLines = FALSE)
writeDescriptiveDataTableAuto(wb2, sheet = 3, x = distribution_treatment, description="Intervention distribution of treatment")


# Table 1
supp_table3_vars <- c("age_lmp", "year_lmp_cat", "alcohol_abuse", "drug_abuse", "tobacco_use", "obese", "overweight_or_obese", "severe_obese", "chronic_hypertension", "type2_diabetes", "hyperglycemia", "dyslipidemia",
                      "infertility", "pcos", "congestive_heart_failure", "ischemic_heart_disease", "renal_disease", "prior_cholesterol_lowering_medication", "prior_fertility_medication", "prior_insulin", "prior_metformin",
                      "prior_alpha_agonist", "prior_beta_blocker", "prior_cc_blocker", "prior_thiazide_diuretic",
                      "prior_proc_hba1c_test", "prior_proc_uacr_test", "prior_proc_lipid_panel", "num_outpatient",
                      "last_ace_inhibitor_gennme", "last_ace_inhibitor_daysupp")


supp_table3 <- acei_pregnancies  %>%
  select(all_of(supp_table3_vars)) %>%
  mutate(last_ace_inhibitor_gennme = str_extract(last_ace_inhibitor_gennme, "([A-Za-z]+pril)")) %>%
  mutate(across(everything(), convert_numeric_to_binary)) %>%
  tbl_summary(include = everything(), statistic=list(all_continuous() ~ "{mean} ({sd})"),
              digits = list(all_categorical() ~ c(0,1), all_continuous() ~ c(1,1))) %>%
  as_pretty_tibble()
supp_table3 %>% flextable()

addWorksheet(wb2, sheetName = "Supp Table 3", gridLines = FALSE)
writeDescriptiveDataTableAuto(wb2, sheet = 4, x = supp_table3, description="Characteristics of study population")

table1_vars <- c("age_lmp", "year_lmp_cat", "alcohol_abuse", "drug_abuse", "tobacco_use", "obese", "overweight_or_obese", "severe_obese", "chronic_hypertension", "type2_diabetes", "hyperglycemia", "dyslipidemia",
                 "infertility", "pcos", "congestive_heart_failure", "ischemic_heart_disease", "renal_disease", "prior_cholesterol_lowering_medication", "prior_fertility_medication", "prior_insulin", "prior_metformin",
                 "prior_alpha_agonist", "prior_beta_blocker", "prior_cc_blocker", "prior_thiazide_diuretic",
                 "prior_proc_hba1c_test", "prior_proc_uacr_test", "prior_proc_lipid_panel", "num_outpatient", "prior_teratogen", "prior_statin",
                 "ft_alpha_agonist", "ft_beta_blocker", "ft_cc_blocker", "ft_thiazide_diuretic", "ft_teratogen", "ft_statin", "ft_metformin", "ft_insulin",
                "last_ace_inhibitor_daysupp", "weeks_supply_at_conception",
                "num_ace_inhibitor_pregnancy_cat", "infant_enrol30to90", "multiple_linked_infant", "last_ace_inhibitor_dispensation"
                 )

table1 <- acei_livebirths %>%
  mutate(pretty_sim_continuer = factor(if_else(sim_continuer == 0, "Non-continuer", "Continuer"), levels=c("Continuer","Non-continuer"))) %>%
  select(pretty_sim_continuer, all_of(table1_vars)) %>%
  tbl_summary(by=pretty_sim_continuer, include = everything(), statistic = list(all_continuous() ~ "{mean} ({sd})"),
              digits = list(all_categorical() ~ c(0,1), all_continuous() ~ c(1,1))) %>%
  add_difference(everything() ~ "smd") %>%
  as_pretty_tibble() %>% select(!`95% CI`) %>%
  rename(ASMD = `Difference`) %>% mutate(ASMD = abs(as.double(ASMD))) %>%
  mutate(ASMD = if_else(str_detect(Characteristic, "Last GLP"), NA_real_, ASMD))
table1 %>% flextable()


addWorksheet(wb2, sheetName = "Table 1", gridLines = FALSE)
writeDescriptiveDataTableAuto(wb2, sheet = 5, x = table1, description="Characteristics of live birth pregnancies")

weighted_table1_vars <- table1_vars[!(table1_vars %in% c("num_ace_inhibitor_pregnancy_cat"))]

table1_weighted <- acei_livebirths %>%
  mutate(pretty_sim_continuer = factor(if_else(sim_continuer == 0, "Non-continuer", "Continuer"), levels=c("Continuer","Non-continuer"))) %>%
  select(newid, pretty_sim_continuer, all_of(weighted_table1_vars), starts_with("ipcw")) %>%
  {survey::svydesign(id=~newid, weights=~ipcw_user, data=.)} %>%
  tbl_svysummary(by=pretty_sim_continuer, include = all_of(weighted_table1_vars), statistic = list(all_continuous() ~ "{mean} ({sd})"),
                 digits = list(all_categorical() ~ c(0,1), all_continuous() ~ c(1,1))) %>%
  add_difference(everything() ~ "smd") %>%
  as_pretty_tibble() %>% select(!`95% CI`) %>%
  rename(ASMD = `Difference`) %>% mutate(ASMD = abs(as.double(ASMD)))
table1_weighted %>% flextable()

table1_weighted_all <- acei_week42 %>%
  mutate(pretty_sim_continuer = factor(if_else(sim_continuer == 0, "Non-continuer", "Continuer"), levels=c("Continuer","Non-continuer"))) %>%
  select(newid, pretty_sim_continuer, all_of(weighted_table1_vars), starts_with("ipcw")) %>%
  {survey::svydesign(id=~newid, weights=~ipcw_user, data=.)} %>%
  tbl_svysummary(by=pretty_sim_continuer, include = all_of(weighted_table1_vars), statistic = list(all_continuous() ~ "{mean} ({sd})"),
                 digits = list(all_categorical() ~ c(0,1), all_continuous() ~ c(1,1))) %>%
  add_difference(everything() ~ "smd") %>%
  as_pretty_tibble() %>% select(!`95% CI`) %>%
  rename(ASMD = `Difference`) %>% mutate(ASMD = abs(as.double(ASMD)))
table1_weighted_all %>% flextable()


addWorksheet(wb2, sheetName = "Table 1 weighted", gridLines = FALSE)
writeDescriptiveDataTableAuto(wb2, sheet = 6, x = table1_weighted, description="Weighted characteristics of live birth pregnancies")

addWorksheet(wb2, sheetName = "Table 1 weighted all", gridLines = FALSE)
writeDescriptiveDataTableAuto(wb2, sheet = 7, x = table1_weighted_all, description="Weighted characteristics of all pregnancies")


# characteristics of linked pregnancies with and without linkage and adequate enrolment
linkage_characteristics <- acei_livebirths %>%
  mutate(pretty_infant_enrol30to90 = if_else(infant_enrol30to90==0, "Not linked with adequate enrolment", "Linked with adequate enrolment")) %>%
  select(pretty_infant_enrol30to90, all_of(table1_vars)) %>%
  tbl_summary(by=pretty_infant_enrol30to90, include = all_of(table1_vars), statistic = list(all_continuous() ~ "{mean} ({sd})"),
              digits = list(all_categorical() ~ c(0,1), all_continuous() ~ c(1,1))) %>%
  add_difference(everything() ~ "smd") %>%
  as_pretty_tibble() %>%
  select(!`95% CI`) %>%
  rename(ASMD = `Difference`) %>% mutate(ASMD = abs(as.double(ASMD)))
linkage_characteristics %>% flextable()

addWorksheet(wb2, sheetName = "Linkage characteristics", gridLines = FALSE)
writeDescriptiveDataTableAuto(wb2, sheet = 8, x = linkage_characteristics, description="Characteristic of pregnancies with and without linkage and adequate infant enrolment")


# Save workbook -------------------------------------------------------------------------------

saveWorkbook(wb2, glue("output/acei_descriptive_publication{YEAR}_{get_formatted_date()}.xlsx"), overwrite = TRUE)

# love plot of GLP-1RA pregnancies by year (not for publication) -------------------------------------------------------------------------------

loveplot_vars <- c("age_lmp", "year_lmp", "weeks_supply_at_conception", overall_covariates)
var_label(acei_livebirths$weeks_supply_at_conception) <- "Weeks' remaining supply at conception"
var_label(acei_livebirths$year_lmp) <- "Year of LMP"
var_label(acei_livebirths$pcos) <- "PMOS"

CairoPDF(glue("output/acei_loveplot{YEAR}.pdf"), width=8, height=8)
love.plot(acei_livebirths %>% select(all_of(loveplot_vars)), treat=acei_livebirths$sim_continuer, weights=acei_livebirths$ipcw_user, binary="std", limits=c(0,0.8), abs=TRUE, thresholds=0.1, var.names=get_variable_labels(acei_livebirths %>% select(all_of(loveplot_vars)))) +
  theme(
    axis.title = element_text(size=14),
    axis.text = element_text(size=14),
    legend.title = element_text(size=14),
    legend.text = element_text(size=14)
  )
  dev.off()

# Create remaining supply plot ----------------------------------------------------------------

create_remaining_supply_plot <- function(livebirths, prefix, weighted) {
  formatted_data <- livebirths %>%
    mutate(pretty_sim_continuer = factor(if_else(sim_continuer == 0, "Non-continuer", "Continuer"), levels=c("Continuer","Discontinuer"))) %>%
    group_by(pretty_sim_continuer)

  if (weighted) {
    formatted_data %<>% summarise(across(starts_with(prefix), ~ sum(ipcw_user*.)), total=sum(ipcw_user)) %>% ungroup()
  } else {
    formatted_data %<>% summarise(across(starts_with(prefix), ~ sum(.)), total=n()) %>% ungroup()
  }

  remaining_supply_plot <- formatted_data %>%
  pivot_longer(starts_with(prefix), names_prefix=prefix, values_to="count", names_to="gestage") %>%
  mutate(gestage=factor(gestage, levels=c(2:12))) %>%
    mutate(perc=100*count/total) %>%
    ggplot(aes(x=gestage, y=perc, fill=pretty_sim_continuer)) + geom_bar(position="dodge", stat="identity") +
    theme_minimal() +
    labs(x="Gestational week", y="% with remaining supply", fill="Group") +
    scale_y_continuous(limits=c(0,100)) +
    theme(
      axis.title = element_text(size=14),
      axis.text = element_text(size=14),
      legend.title = element_text(size=12),
      legend.text = element_text(size=12),
    )
  return(remaining_supply_plot)
}

remaining_prior_unweighted_plot <- acei_livebirths %>%
  create_remaining_supply_plot("remaining_prior_weeks_supply_week", weighted=FALSE)

remaining_prior_weighted_plot <- acei_livebirths %>%
  create_remaining_supply_plot("remaining_prior_weeks_supply_week", weighted=TRUE)

remaining_unweighted_plot <- acei_livebirths %>%
  create_remaining_supply_plot("remaining_weeks_supply_week", weighted=FALSE)

remaining_weighted_plot <- acei_livebirths %>%
  create_remaining_supply_plot("remaining_weeks_supply_week", weighted=TRUE)

CairoPDF(glue("output/remaining_supply_ace_inhibitor{YEAR}.pdf"), width=12, height=8)
plot_grid(
  plot_grid(remaining_prior_unweighted_plot + theme(legend.position = "none"),
          remaining_prior_weighted_plot + theme(legend.position = "none"),
          remaining_unweighted_plot + theme(legend.position = "none"),
          remaining_weighted_plot + theme(legend.position = "none"),
          labels=c("A", "B", "C", "D"),nrow=2, ncol=2),
  get_legend(remaining_prior_unweighted_plot),
  ncol=2,
  rel_widths = c(3,0.6)
  )
dev.off()


remaining_prior_unweighted_plot_w6 <- acei_livebirths_w6 %>%
  create_remaining_supply_plot("remaining_prior_weeks_supply_week", weighted=FALSE)

remaining_prior_weighted_plot_w6 <- acei_livebirths_w6 %>%
  create_remaining_supply_plot("remaining_prior_weeks_supply_week", weighted=TRUE)

remaining_unweighted_plot_w6 <- acei_livebirths_w6 %>%
  create_remaining_supply_plot("remaining_weeks_supply_week", weighted=FALSE)

remaining_weighted_plot_w6 <- acei_livebirths_w6 %>%
  create_remaining_supply_plot("remaining_weeks_supply_week", weighted=TRUE)

CairoPDF(glue("output/remaining_supply_w6_ace_inhibitor{YEAR}.pdf"), width=12, height=8)
plot_grid(
  plot_grid(remaining_prior_unweighted_plot_w6 + theme(legend.position = "none"),
          remaining_prior_weighted_plot_w6 + theme(legend.position = "none"),
          remaining_unweighted_plot_w6 + theme(legend.position = "none"),
          remaining_weighted_plot_w6 + theme(legend.position = "none"),
          labels=c("A", "B", "C", "D"),nrow=2, ncol=2),
  get_legend(remaining_prior_unweighted_plot_w6),
  ncol=2,
  rel_widths = c(3,0.6))
dev.off()


# Create alluvial plot ------------------------------------------------------------------------

create_alluvial_counts <- function(pregnancies) {
  pregnancies %<>% mutate(tt_drug_type=if_else(enddate <= lmp + ddays(188), "ST pregnancy loss", tt_drug_type))
  pregnancies %<>% mutate(st_drug_type=if_else(enddate <= lmp + ddays(90), "FT pregnancy loss", st_drug_type))
  alluvial_counts <- pregnancies %>%
    group_by(prior_drug_type, ft_drug_type, st_drug_type) %>%
    summarise(count=n(), count_livebirth=sum(final_type < 3), count_linked_livebirth=sum(infant_enrol30to90, na.rm=TRUE)) %>%
    arrange(desc(count)) %>%
    ungroup()
  return(alluvial_counts)
}
alluvial_counts <- create_alluvial_counts(acei_pregnancies)

alluvial_plot <- alluvial_counts %>%
  filter(count_livebirth > 75) %>%
  mutate(across(c(prior_drug_type, ft_drug_type, st_drug_type), ~ str_replace_all(., "&", "&\n"))) %>%
  mutate(across(c(prior_drug_type, ft_drug_type, st_drug_type), ~ factor(., levels = sort(unique(c(prior_drug_type, ft_drug_type, st_drug_type)))))) %>%
  complete(prior_drug_type, ft_drug_type, st_drug_type) %>% mutate(count_livebirth = replace_na(count_livebirth, 0)) %>%
  mutate(id = row_number()) %>%
  pivot_longer(cols = c(prior_drug_type, ft_drug_type, st_drug_type), names_to = "time", values_to="drug_type") %>%
  mutate(drug_type = factor(drug_type), time=factor(time, levels=c("prior_drug_type", "ft_drug_type", "st_drug_type"))) %>%
  mutate(time = fct_recode(time, "Prior" = "prior_drug_type", "First trimester" = "ft_drug_type", "Second trimester" = "st_drug_type")) %>%
  ggplot(aes(x = time, stratum=drug_type, y = count_livebirth, fill=drug_type, label=drug_type, alluvium = id)) +
  geom_flow() +
  geom_stratum(alpha = .5) +
  #  geom_text(stat = "stratum", size = 3) +
  xlab("Period") +
  ylab("No. pregnancies") +
  labs(fill="Antihypertensive") +
  theme_minimal() +
  theme(
    axis.title = element_text(size=14),
    axis.text = element_text(size=14),
    legend.title = element_text(size=12),
    legend.text = element_text(size=12),
  )
CairoPDF(glue("output/alluvial_ace_inhibitor{YEAR}_pptx.pdf"), width=9, height=6)
alluvial_plot
dev.off()


# Inspect teratogens --------------------------------------------------------------------------

pregnancies_teratogen <- read_parquet(glue("{COMP_DATA}/pregnancies2025_w_teratogens.parquet"))
acei_pregnancies %>%
  select(pregid) %>%
  left_join(pregnancies_teratogen, by="pregid") %>%
  summarise(across(starts_with("prior"), sum)) %>%
  transpose_tibble("teratogen") %>%
  arrange(desc(`0`))
write_csv(glue("output/acei_teratogen{YEAR}.csv"))

# Boxplot of weights --------------------------------------------------------------------------

CairoPDF(glue("output/acei_weight_distribution{YEAR}.pdf"), width=6, height=6)
ggplot(acei_week42 %>% mutate(ln_ipcw_user = log(ipcw_user)) %>%
         mutate(pretty_sim_continuer = factor(if_else(sim_continuer == 0, "Non-continuer", "Continuer"), levels=c("Continuer","Non-continuer"))),
       aes(pretty_sim_continuer, ln_ipcw_user)) +
  geom_boxplot() +
  xlab("Continuation status") +
  ylab("Log weight") +
  stat_boxplot(geom = "errorbar", width=0.2) +
  theme_classic() +
  scale_y_continuous(breaks=c(-3:3)) +
  theme(
    axis.title = element_text(size=14),
    axis.text = element_text(size=14),
    legend.title = element_text(size=14),
    legend.text = element_text(size=14)
  )
dev.off()
