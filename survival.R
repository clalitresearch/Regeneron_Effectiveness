library(survival)
library(broom)
library(survminer)
library(dplyr)
library(officer)
     
## -----------------------------------------------------------------------
##                  calc_survival_variables
## -----------------------------------------------------------------------
# This function is given a raw dataframe,
#                        an indication of an outcome column (Hospitalization / Severe / Death)
#                        and indication of the required analysis type (full or one of 2 subgroup analysis)
# and returns a new df with the required columns for survival analysis
# ------------------------------------------------------------------------------------------------------------------------------
calc_survival_variables = function(df_raw, outcome_date_col, sub_grp_anlsys, params) {
  administrative_fwp_end_date =  max(as.Date((df_raw$retrieval_date)))
  fwp_days = params$get('fwp_days')
  
  
  df =  
    df_raw %>%
    group_by(subclass) %>%
    mutate (infection_to_regeneron_days = max(as.numeric(difftime(regeneron_date, index_date, units = "days")), na.rm = T)) %>%
    ungroup() %>%
    # virtual treatment - so can be defined both for the treated an the untreated
    mutate(virt_treatment_date = index_date + infection_to_regeneron_days) %>%
    rename(curr_outcome = {{ outcome_date_col }}) %>% 
    mutate(days_to_event = as.numeric(difftime(curr_outcome, virt_treatment_date, units = "days"))) %>% 
    # Note days_to_censoring ignored the possibility of an event before censoring
    mutate (days_to_censoring  = as.numeric(difftime(administrative_fwp_end_date,    virt_treatment_date, units = "days"))) %>%
    # we limit the followup to fwp_days (30) days
    mutate (days_to_censoring  = pmin(days_to_censoring, fwp_days)) %>% 
    mutate (days_of_fwp        = pmin(days_to_censoring, days_to_event, na.rm = T)) %>%
    mutate (is_event           = case_when(days_of_fwp == days_to_event ~ 1, T ~ 0)) %>%
    rename(rowid = subclass)

  result_df = df %>%
    mutate(sw_heart = pmax(psw_cdc_possible_cereb_vasc, psw_heart)) %>% 
    mutate(sw_diabetes = pmax(psw_dm_type1, psw_dm_type2)) %>% 
    mutate(sw_respiratory = pmax(psw_cdc_possible_asthma, psw_cdc_certain_copd, psw_cdc_possible_other_resp)) %>% 
    select(rowid, days_of_fwp, is_event, exposure_cat, num_vacc_doses,
           age, sex, sector_cat, SES,
           pn_flu_vacc_5yr, pn_pcr_tests_9m, bmi_cat, pc_smoke_desc, pn_pack_years, num_cdc_risk_factors, sw_imunosuppressed,
           psw_preg,  simplified_vacc1_ww, sw_waned_immunity, sw_heart, sw_diabetes, sw_respiratory,psw_cdc_possible_htn,
           psw_cdc_certain_ckd, psw_cdc_certain_cancer, psw_cdc_possible_neuro, psw_hepatic)
  
  
  if (sub_grp_anlsys == '60_and_above'){
    result_df = result_df %>%
      filter(age >= 60)
  }
  if (sub_grp_anlsys == 'below_60'){
    result_df = result_df %>%
      filter(age < 60)
  }
 
  return (result_df)
}


## -----------------------------------------------------------------------
get_cox_output = function(df, formula_str, ph_test=F, filename = NA) {
## -----------------------------------------------------------------------
  formula = as.formula(formula_str)
  coxmodel =  coxph(formula, data=df )
  
  summary = broom::tidy(coxmodel, conf.int=TRUE) %>%
  mutate(across(c("estimate","conf.low","conf.high"), exp)) %>% 
  select(term, estimate, std.error, conf.low, conf.high, p.value) %>% 
  rename(HR = estimate) %>% 

  mutate(across(where(is.numeric), ~round(.x, 3)))
  
  outcome = deparse(formula)[1]

  if (!is.na(filename)) {
  read_docx() %>%  # a new, empty document
    body_add_table(summary, style = "table_template") %>% 
    print(target= glue("{filename}"))
  }

  print(summary %>% kable())
  y =  Surv(df$days_of_fwp, df$is_event) 
  with(df, plot(survfit(y ~ exposure_cat ), col=c("red", "blue"),
                fun="cloglog",
       xlab = "Time (days)",
       ylab = "log-log of S(t)",
       main = "Complementary log-log plot"))

  plt = ggcoxdiagnostics(coxmodel, type = "deviance")
  plot(plt)
  if (ph_test) {
    test.ph = cox.zph(coxmodel)
      ggcoxzph(test.ph)
  }
}
