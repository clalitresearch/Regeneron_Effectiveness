---
title: "Table_1"
format: docx
editor: visual
---

```{r}
pacman::p_load(tidyverse, cri.utils, skimr, gtsummary)
cri.utils::activate_keytab("tommus")
con <- cri.utils::cri_connection()
Data <- connect_table(con, "dbo.[scheme]", "table_name") %>%
  dplyr::collect() %>%
  mutate(pn_flu_vacc_5yr = as.numeric(as.character(pn_flu_vacc_5yr)),
        historic_DT2_diag = if_else(as.numeric(historic_DT2_diag) == 2, 1, as.numeric(historic_DT2_diag)),
        Condition = as_factor(Condition)) %>%
  mutate(Condition = fct_relevel(Condition, c('not_active', 'Active'))) %>%
  mutate(non_smoker = if_else(Smoking == 1, 1,
                                        if_else(Smoking > 1, 0, NA_real_)))


```

# Table 1 by condition general pop

```{r, message=FALSE,warning=FALSE}

# started with 9867 in active, finished before matching 6667
table_1 <- Data %>%
  select(age_at_index_date, Gender, Sector, bmi_value, BMI_CAT, 
         historicCVD, historic_DT2_diag, historic_diab_complication, 
         kod_machoz, SUG_YISHUV, HbA1c_Base_HbA1c_CAT, HbA1c_base_BT,
          Glucose_base_BT, HDL_base_BT, LDL_base_BT, Trig_base_BT,
         corona_COM, corona_COM_CAT, pstr_charlson_with_age, Smoking, SES,
          
          HTN_diag, COPD_diag, astma_diag, Bronchiectasis_diag, Pulmonary_fibrosis_diag, CKD_diag, hyperlipidemia_diag,
           pn_flu_vacc_5yr, pn_gp_visits_5yr, Condition) %>%
   mutate(Smoking = as.factor(if_else(Smoking == 1, 'not_smoker',
                           if_else(Smoking == 2, 'past_smoker',
                                   if_else(Smoking == 3, 'Smoker', NA_character_)))),
         SES = as.factor(if_else(SES == 1, 'Low',
                          if_else(SES == 2, 'Moderate', 
                               if_else(SES == 3, 'High', NA_character_)))),
         corona_COM_CAT = as.factor(if_else(corona_COM_CAT == 0, '0',
                                  if_else(corona_COM_CAT == 1, '1_2',
                                          if_else(corona_COM_CAT == 2, '3_4',
                                                  if_else(corona_COM_CAT == 3, '5', NA_character_))))),
         Sector = as.factor(if_else(Sector == 1, 'Arab',
                          if_else(Sector == 2, 'Ultra_orthodox',
                                  if_else(Sector == 9, 'General', NA_character_))))) %>% 
  mutate(BMI_CAT        = fct_relevel(BMI_CAT, c('UNDER', 'NORMAL', 'OVER', 'OBESE', 'MORBID_OBESE')),
         Smoking        = fct_relevel(Smoking, c('not_smoker', 'past_smoker', 'Smoker')),
         SES            = fct_relevel(SES, c('Low', 'Moderate', 'High')),
         corona_COM_CAT = fct_relevel(corona_COM_CAT, c('0', '1_2', '3_4', '5')),
         Sector         = fct_relevel(Sector, c('General', 'Ultra_orthodox', 'Arab')),
         kod_machoz = as.factor(as.numeric(kod_machoz))) %>% 
  tbl_summary(by = 'Condition',
            type = c(pn_flu_vacc_5yr = "continuous"),
            statistic  = list(all_continuous() ~ "{mean} ({sd}), {median}, {IQR}"))

table_1

```

```{r, message=FALSE,warning=FALSE}

# started with 9867 in active, finished before matching 6667
table_2 <- Data %>%
  filter(Condition == 'Active') %>% 
  mutate(did_leave = if_else(is.na(left_app_date), 'not_leave', 'leave')) %>% 
  select(age_at_index_date, Gender, Sector, bmi_value, BMI_CAT, 
         historicCVD, historic_DT2_diag, historic_diab_complication, 
          SUG_YISHUV, HbA1c_Base_HbA1c_CAT, HbA1c_base_BT,
          Glucose_base_BT, HDL_base_BT, LDL_base_BT, Trig_base_BT,
         corona_COM, corona_COM_CAT, pstr_charlson_with_age, Smoking, SES,
          
          HTN_diag, COPD_diag, astma_diag, Bronchiectasis_diag, Pulmonary_fibrosis_diag, CKD_diag, hyperlipidemia_diag,
           pn_flu_vacc_5yr, pn_gp_visits_5yr, did_leave) %>%
   mutate(Smoking = as.factor(if_else(Smoking == 1, 'not_smoker',
                           if_else(Smoking == 2, 'past_smoker',
                                   if_else(Smoking == 3, 'Smoker', NA_character_)))),
         SES = as.factor(if_else(SES == 1, 'Low',
                          if_else(SES == 2, 'Moderate', 
                               if_else(SES == 3, 'High', NA_character_)))),
         corona_COM_CAT = as.factor(if_else(corona_COM_CAT == 0, '0',
                                  if_else(corona_COM_CAT == 1, '1_2',
                                          if_else(corona_COM_CAT == 2, '3_4',
                                                  if_else(corona_COM_CAT == 3, '5', NA_character_))))),
         Sector = as.factor(if_else(Sector == 1, 'Arab',
                          if_else(Sector == 2, 'Ultra_orthodox',
                                  if_else(Sector == 9, 'General', NA_character_))))) %>% 
  mutate(BMI_CAT        = fct_relevel(BMI_CAT, c('UNDER', 'NORMAL', 'OVER', 'OBESE', 'MORBID_OBESE')),
         Smoking        = fct_relevel(Smoking, c('not_smoker', 'past_smoker', 'Smoker')),
         SES            = fct_relevel(SES, c('Low', 'Moderate', 'High')),
         corona_COM_CAT = fct_relevel(corona_COM_CAT, c('0', '1_2', '3_4', '5')),
         Sector         = fct_relevel(Sector, c('General', 'Ultra_orthodox', 'Arab'))
         # kod_machoz = as.factor(as.numeric(kod_machoz))
         ) %>% 
  tbl_summary(by = 'did_leave',
            type = c(pn_flu_vacc_5yr = "continuous"),
            statistic  = list(all_continuous() ~ "{mean} ({sd}), {median}, {IQR}")) %>% 
   add_p()

table_2

```
