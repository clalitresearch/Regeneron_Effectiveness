---
title: "Survival by activity groups"
format: html
editor: visual
---

```{r}
pacman::p_load(tidyverse, cri.utils, skimr, gtsummary, ggfortify, ggplot2, gmodels, survival, grDevices, patchwork, smoothHR)
cri.utils::activate_keytab("tommuspim") # innovationprod
con <- cri.utils::cri_connection() 

Data <- connect_table(con, "db.[scheme]", "table_name") %>%
  collect() %>%
      mutate(imputed_pn_gp_visits_5yr = ifelse(is.na(pn_gp_visits_5yr), round(median(pn_gp_visits_5yr, na.rm = T)), pn_gp_visits_5yr)
             , historic_DT2_diag = if_else(historic_DT2_diag == 2, 1, 0))


Data_steps <- Data %>% 
  drop_na(activity_level) %>% 
  mutate(activity_groups = if_else(activity_level == 'Low', 'Mild', activity_level)) %>% 
  mutate(activity_groups_num = as_factor(case_when(
    Condition == 'not_active' ~ 0
    , activity_groups == 'Mild' ~ 1
    , activity_groups == 'Moderate' ~ 2
    , T ~ 3
    )))

```

```{r}


remainPaires_surv <- function(data){
  
  Complic <- data %>% 
    group_by(subclass) %>% 
    mutate(n = length(subclass)) %>% 
    filter(n < 2) %>% 
    distinct(subclass) %>% 
    pull(subclass)
  
return(Complic)
  
}

prepareData <- function(data
                        , history
                        , occuranceDate
                        , titleName){ 
  
        data <- data %>%
            filter({{history}} == 0)
    
     if( titleName %in% c('CVA')){

      data <- data %>%
            group_by(subclass) %>%
            mutate(lead_cvd = lead({{history}})
             , lag_cvd = lag({{history}})) %>%
             mutate(summing = if_else(({{history}} + lead_cvd == 1) | ({{history}} + lag_cvd == 1), 'remove', 'ok')) %>%
            filter(is.na(summing)) %>%
            ungroup()

    }

  invalid_paires <- remainPaires_surv(data = data)
  data <- data %>%
      filter(! subclass %in% invalid_paires)
    
  data <- data %>%
       mutate(newOccurance = if_else(is.na({{occuranceDate}}), 0, 1)) %>% # 
       mutate(Duration     = if_else(newOccurance == 0, as.numeric(loss_to_FU - index_date), as.numeric({{occuranceDate}} - index_date))) %>%    
       mutate(newOccurance = if_else(Duration > 1000, 0, newOccurance)
              , Duration     = if_else(Duration > 1000, 1000, Duration))
   
    return(data)
}
 


give_HR <- function(data) {
  
HR <- data %>%
  broom::tidy( conf.int = TRUE) %>%
    mutate(across(c("estimate", "conf.low", "conf.high"), exp)) %>%
    mutate(across(where(is.numeric), ~ round(.x, 2))) %>%
    select(-statistic, -'std.error') %>%
    as.data.frame() 

return(HR)
  
}



create_figure <- function(the_model, accuracy, high_limit, title, colors, legend.position) {

plot <- survminer::ggsurvplot(the_model, fun = 'cumhaz',xlab = "", ylab = c(""), legend.title = "",conf.int = TRUE, legend = "right", legend.labs = c('Non Active users', 'Mild', 'Moderate', 'High'), size = 0.05, censor.size = 0.01, xlim = c(0, 1000), font.main = c(20, "bold"), surv.scale = c("percent") ,  risk.table = TRUE)

plot_ch <- plot$plot 

plt_1 <- 
plot_ch$data %>% 
  mutate(lt = case_when(
    strata ==  "Non-app users" ~ 'solid'
    , strata ==  "Mild" ~ 'Dotdash '
    , strata ==  "Moderate" ~ 'dotted'
    , T  ~ 'dashed')
  ) %>% 
  ggplot(aes(x = time, y = surv, group = strata, colour = strata, linetype = lt)) + 
  theme_classic() +
  geom_line() +
     scale_color_manual(values =  c('black', 'gray20', 'gray60', 'gray90')) +
  ggtitle(title) +
  theme(axis.text.x    = element_text(size = 20)
        , axis.text.y   = element_text(size = 20)
        , axis.title.y  = element_text(size = 20)
        , plot.title    = element_text(size = 20, face = 'bold')
        , legend.text   = element_text(size = 20)) + 
  xlab('') + ylab('') +
  theme(legend.position = legend.position) +
  scale_y_continuous(labels = scales::percent_format(accuracy = accuracy), limits = c(0, high_limit)) + 
  scale_x_continuous(breaks = c(0,  200, 400, 600, 800, 1000))  + guides(color = guide_legend(title = NULL))  +
  guides(linetype = 'none')
 
  return(plt_1)
  
}



```

## CVD

```{r}

activeCVD_steps <- prepareData(Data_steps
                         , history = historicCVD
                         , occuranceDate = CVDoccuranceDate
                         , titleName = 'CVD')

activeCVD_steps %>% select(activity_groups_num, newOccurance) %>% table()
 un_model_cvd_steps <- coxph(Surv(Duration, newOccurance) ~ activity_groups_num, data = activeCVD_steps)
 give_HR(un_model_cvd_steps)
 
 
model_cvd_steps <- coxph(Surv(Duration, newOccurance) ~ activity_groups_num +
                       age_at_index_date +
                       I(age_at_index_date * age_at_index_date) +
                       bmi_value +
                       Gender +
                       log(imputed_pn_gp_visits_5yr) +
                       as.factor(SES) +
                       historic_DT2_diag  +
                       HTN_diag +
                       hyperlipidemia_diag ,
                   data = activeCVD_steps)
 
S_R(model_cvd_steps)

hr_model_cvd_steps <- give_HR(model_cvd_steps) %>% 
   mutate(case = 'CVD') %>% 
   head(3)

hr_model_cvd_steps

activeCVD_steps <- activeCVD_steps %>% 
  mutate(activity_groups = if_else(is.na(activity_groups), 'Non Active', activity_groups)) %>% 
  mutate(activity_groups = as.factor(activity_groups)) %>% 
  mutate(activity_groups = fct_relevel(activity_groups, c('Non Active', 'Mild', 'Moderate', 'High')))

model_cvd_steps_fit <- survival::survfit(survival::Surv(Duration, newOccurance) ~ activity_groups,
                               data = activeCVD_steps, id = guid_tz)


```

```{r}
activeCVD_steps <- activeCVD_steps %>% 
  mutate(activity_groups = as.character(activity_groups)
    , new_groups = as.factor(ifelse(activity_groups_num == 0, 'No Active', activity_groups))) %>% 
  mutate(new_groups = fct_relevel(new_groups, c('No Active', 'Mild', 'Moderate', 'High')))


ch_cvd <- survival::survfit(survival::Surv(Duration, newOccurance) ~ new_groups,
                               data = activeCVD_steps, id = guid_tz)
 
plt_1 = create_figure(the_model = ch_cvd, accuracy = 1, high_limit = 0.05, title = "C. Cardiovascular Disease", legend.position = "none")

```

## diabetes

```{r}
activeDiab_steps <- prepareData(Data_steps
                          , history = historic_DT2_diag
                          , occuranceDate = diab_occurance_date
                          , 'diabetes')
activeDiab_steps %>% select(activity_groups_num, newOccurance) %>% table()

 model_diab_steps <- coxph(Surv(Duration, newOccurance) ~ activity_groups_num , data = activeDiab_steps)
  give_HR(model_diab_steps)                         

 model_diab_steps <- coxph(Surv(Duration, newOccurance) ~ activity_groups_num +
                       age_at_index_date +
                       I(age_at_index_date * age_at_index_date) +
                       bmi_value +
                       Gender +
                       log(imputed_pn_gp_visits_5yr) +
                       as.factor(SES) +
                       historicCVD  +
                       HTN_diag +
                       hyperlipidemia_diag ,
                   data = activeDiab_steps)
 
S_R(model_diab_steps)

hr_model_diab_steps <- give_HR(model_diab_steps)%>% 
  mutate(case = 'Diabetes') %>% 
  head(3)

hr_model_diab_steps




```

```{r}
activeDiab_steps <- activeDiab_steps %>% 
  mutate(activity_groups = as.character(activity_groups)
        , new_groups = as.factor(ifelse(activity_groups_num == 0, 'No Active', activity_groups))) %>% 
  mutate(new_groups = fct_relevel(new_groups, c('No Active', 'Mild', 'Moderate', 'High')))


ch_diab <- survival::survfit(survival::Surv(Duration, newOccurance) ~ new_groups,
                               data = activeDiab_steps, id = guid_tz)

plt_2 = create_figure(the_model = ch_diab, accuracy = 0.1, high_limit = 0.012, title = "B. Diabetes Mellitus", legend.position = "none")
  
```

## CVA

```{r}
activeStroke_steps <- prepareData(Data_steps
                            , history = stroke_history
                            , occuranceDate = stroke_date
                            , titleName = 'CVA')
activeStroke_steps %>% select(activity_groups_num, newOccurance) %>% table()


model_CVA_steps <- coxph(Surv(Duration, newOccurance) ~ activity_groups_num ,data= activeStroke_steps)
give_HR(model_CVA_steps)
 
 model_CVA_steps <- coxph(Surv(Duration, newOccurance) ~ activity_groups_num +
                       age_at_index_date +
                       I(age_at_index_date * age_at_index_date) +
                       bmi_value +
                       Gender +
                       log(imputed_pn_gp_visits_5yr) +
                       as.factor(SES) +
                       historic_DT2_diag  +
                       historicCVD +
                       HTN_diag +
                       AF_diag,
                      data = activeStroke_steps)
 
S_R(model_CVA_steps)

hr_model_CVA_steps <-  give_HR(model_CVA_steps) %>% 
    mutate(case = 'CVA') %>% 
    head(3)
hr_model_CVA_steps
```

```{r}
activeStroke_steps <- activeStroke_steps %>% 
  mutate(activity_groups = as.character(activity_groups)
    , new_groups = as.factor(ifelse(activity_groups_num == 0, 'Non Active users', activity_groups))) %>% 
  mutate(new_groups = fct_relevel(new_groups, c('Non Active users', 'Mild', 'Moderate', 'High')))


ch_cva <- survival::survfit(survival::Surv(Duration, newOccurance) ~ new_groups,
                               data = activeStroke_steps, id = guid_tz)
 
plt_3 = create_figure(the_model = ch_cva, accuracy = 0.1, high_limit = 0.012, title = "A. Stroke", legend.position = "none")
```

## Any diagnosis

```{r}
active_any_steps <- prepareData(Data_steps
                          , history = any_diag_history
                          , occuranceDate = any_diag_date
                          , titleName = 'Any Diagnosis')

active_any_steps %>% select(activity_groups_num, newOccurance) %>% table()


model_any_steps <- coxph(Surv(Duration, newOccurance) ~ activity_groups_num , data = active_any_steps)
S_R(model_any_steps)
give_HR(model_any_steps)

model_any_steps <- coxph(Surv(Duration, newOccurance) ~ activity_groups_num +
                       age_at_index_date +
                       I(age_at_index_date * age_at_index_date) +
                       bmi_value +
                       Gender +
                       log(imputed_pn_gp_visits_5yr) +
                       as.factor(SES) +
                       HTN_diag +
                       hyperlipidemia_diag,
                       data = active_any_steps)
 
S_R(model_any_steps)

hr_model_any_steps <- give_HR(model_any_steps)%>% 
  mutate(case = 'Any Diagnosis') %>% 
  head(3)
hr_model_any_steps
```

```{r}
active_any_steps <- active_any_steps %>% 
  mutate(activity_groups = as.character(activity_groups)
    , new_groups = as.factor(ifelse(activity_groups_num == 0, 'Non Active users', activity_groups))) %>% 
  mutate(new_groups = fct_relevel(new_groups, c('Non Active users', 'Mild', 'Moderate', 'High')))

ch_any <- survival::survfit(survival::Surv(Duration, newOccurance) ~ new_groups,
                               data = active_any_steps, id = guid_tz)
 
plt_4 = create_figure(the_model = ch_any, accuracy = 1, high_limit = 0.05, title = "D. Any Diagnosis", legend.position = "right")
```

```{r}
( plt_3 + plt_2  + plt_1 + plt_4)
getwd()
ggsave(plot = last_plot(), filename = 'acivity_groups.jpeg', dpi = 300, width = 40, height= 30, units = 'cm')

```

# By gender

## CVD by gender

### Males

```{r}
activeCVD_steps_M <- prepareData(Data_steps[Data_steps$Gender == 'M', ]
                         , history = historicCVD
                         , occuranceDate = CVDoccuranceDate
                         , titleName = 'CVD')

activeCVD_steps_M %>% select(activity_groups_num, newOccurance) %>% table()
 model_cvd_M_steps <- coxph(Surv(Duration, newOccurance) ~ activity_groups_num , data = activeCVD_steps_M)
S_R(model_any_steps)
give_HR(model_cvd_M_steps)

 model_cvd_steps_M <- coxph(Surv(Duration, newOccurance) ~ activity_groups_num +
                       age_at_index_date +
                       I(age_at_index_date * age_at_index_date) +
                       bmi_value +
                       log(imputed_pn_gp_visits_5yr) +
                       as.factor(SES) +
                       historic_DT2_diag  +
                       HTN_diag +
                       hyperlipidemia_diag ,
                   data = activeCVD_steps_M)
 
S_R(model_cvd_steps_M)

hr_model_cvd_steps_M <- give_HR(model_cvd_steps_M) %>% 
   mutate(case = 'CVD'
          , Gender = 'Male') %>% 
  head(3)
```

### Feales

```{r}
activeCVD_steps_F <- prepareData(Data_steps[Data_steps$Gender == 'F', ]
                         , history = historicCVD
                         , occuranceDate = CVDoccuranceDate
                         , titleName = 'CVD')

activeCVD_steps_F %>% select(activity_groups_num, newOccurance) %>% table()
 model_cvd_F_steps <- coxph(Surv(Duration, newOccurance) ~ activity_groups_num , data = activeCVD_steps_F)
S_R(model_any_steps)
give_HR(model_cvd_F_steps)

 model_cvd_steps_F <- coxph(Surv(Duration, newOccurance) ~ activity_groups_num +
                       age_at_index_date +
                       I(age_at_index_date * age_at_index_date) +
                       bmi_value +
                       log(imputed_pn_gp_visits_5yr) +
                       as.factor(SES) +
                       historic_DT2_diag  +
                       HTN_diag +
                       hyperlipidemia_diag ,
                   data = activeCVD_steps_F)
 
S_R(model_cvd_steps_F)

hr_model_cvd_steps_F <- give_HR(model_cvd_steps_F) %>% 
   mutate(case = 'CVD'
          , Gender = 'Female') %>% 
  head(3)
hr_model_cvd_steps_F
```

## Diabetes

### Males

```{r}
activeDiabetes_steps_M <- prepareData(Data_steps[Data_steps$Gender == 'M', ]
                         , history = historic_DT2_diag
                         , occuranceDate = diab_occurance_date
                         , titleName = 'diabetes')
activeDiabetes_steps_M %>% select(activity_groups_num, newOccurance) %>% table()
 model_daib_M_steps <- coxph(Surv(Duration, newOccurance) ~ activity_groups_num , data = activeDiabetes_steps_M)
S_R(model_any_steps)
give_HR(model_daib_M_steps)

 model_Diabetes_steps_M <- coxph(Surv(Duration, newOccurance) ~ activity_groups_num +
                       age_at_index_date +
                       I(age_at_index_date * age_at_index_date) +
                       bmi_value +
                       log(imputed_pn_gp_visits_5yr) +
                       as.factor(SES) +
                       historicCVD  +
                       HTN_diag +
                       hyperlipidemia_diag ,
                   data = activeDiabetes_steps_M)
 
S_R(model_Diabetes_steps_M)

hr_model_Diabetes_steps_M <- give_HR(model_Diabetes_steps_M) %>% 
   mutate(case = 'Diabetes'
          , Gender = 'Male') %>% 
  head(3)
hr_model_Diabetes_steps_M
```

### Feales

```{r}
activeDiabetes_steps_F <- prepareData(Data_steps[Data_steps$Gender == 'F', ]
                         , history = historic_DT2_diag
                         , occuranceDate = diab_occurance_date
                         , titleName = 'diabetes')
activeDiabetes_steps_F %>% select(activity_groups_num, newOccurance) %>% table()

 model_daib_F_steps <- coxph(Surv(Duration, newOccurance) ~ activity_groups_num , data = activeDiabetes_steps_F)
S_R(model_any_steps)
give_HR(model_daib_F_steps)
 model_Diabetes_steps_F <- coxph(Surv(Duration, newOccurance) ~ activity_groups_num +
                       age_at_index_date +
                       I(age_at_index_date * age_at_index_date) +
                       bmi_value +
                       log(imputed_pn_gp_visits_5yr) +
                       as.factor(SES) +
                       historicCVD  +
                       HTN_diag +
                       hyperlipidemia_diag ,
                   data = activeDiabetes_steps_F)
 
S_R(model_Diabetes_steps_F)

hr_model_Diabetes_steps_F <- give_HR(model_Diabetes_steps_F) %>% 
   mutate(case = 'Diabetes'
          , Gender = 'Female') %>% 
  head(3)
hr_model_Diabetes_steps_F
```

## CVA by gender

### Males

```{r}
activeCVA_steps_M <- prepareData(Data_steps[Data_steps$Gender == 'M', ]
                         , history = stroke_history
                         , occuranceDate = stroke_date
                         , titleName = 'CVA')
activeCVA_steps_M %>% select(activity_groups_num, newOccurance) %>% table()

 model_CVA_M_steps <- coxph(Surv(Duration, newOccurance) ~ activity_groups_num , data = activeCVA_steps_M)
S_R(model_any_steps)
give_HR(model_CVA_M_steps)

 model_CVA_steps_M <- coxph(Surv(Duration, newOccurance) ~ activity_groups_num +
                       age_at_index_date +
                       I(age_at_index_date * age_at_index_date) +
                       bmi_value +
                       log(imputed_pn_gp_visits_5yr) +
                       as.factor(SES) +
                       historic_DT2_diag  +
                       HTN_diag +
                       AF_diag ,
                   data = activeCVA_steps_M)
 
S_R(model_CVA_steps_M)

hr_model_CVA_steps_M <- give_HR(model_CVA_steps_M) %>% 
   mutate(case = 'Stroke'
          , Gender = 'Male') %>% 
  head(3)
hr_model_CVA_steps_M
```

### Feales

```{r}
activeCVA_steps_F <- prepareData(Data_steps[Data_steps$Gender == 'F', ]
                         , history = stroke_history
                         , occuranceDate = stroke_date
                         , titleName = 'CVA')
activeCVA_steps_F %>% select(activity_groups_num, newOccurance) %>% table()
 model_CVA_F_steps <- coxph(Surv(Duration, newOccurance) ~ activity_groups_num , data = activeCVA_steps_F)
S_R(model_any_steps)
give_HR(model_CVA_F_steps)

 model_CVA_steps_F <- coxph(Surv(Duration, newOccurance) ~ activity_groups_num +
                       age_at_index_date +
                       I(age_at_index_date * age_at_index_date) +
                       bmi_value +
                       log(imputed_pn_gp_visits_5yr) +
                       as.factor(SES) +
                       historic_DT2_diag  +
                       HTN_diag +
                       AF_diag ,
                   data = activeCVA_steps_F)
 
S_R(model_CVA_steps_F)

hr_model_CVA_steps_F <- give_HR(model_CVA_steps_F) %>% 
   mutate(case = 'Stroke'
          , Gender = 'Female') %>% 
  head(3)
hr_model_CVA_steps_F
```

## Any diagnosis

### Male

```{r}
active_any_steps_M <- prepareData(Data_steps[Data_steps$Gender == 'M', ]
                          , history = any_diag_history
                          , occuranceDate = any_diag_date
                          , titleName = 'Any Diagnosis')
active_any_steps_M %>% select(activity_groups_num, newOccurance) %>% table()
 model_any_M_steps <- coxph(Surv(Duration, newOccurance) ~ activity_groups_num , data = active_any_steps_M)
S_R(model_any_steps)
give_HR(model_any_M_steps)


 model_any_steps_M <- coxph(Surv(Duration, newOccurance) ~ activity_groups_num +
                       age_at_index_date +
                       I(age_at_index_date * age_at_index_date) +
                       bmi_value +
                       log(imputed_pn_gp_visits_5yr) +
                       as.factor(SES) +
                       
                       HTN_diag +
                       hyperlipidemia_diag ,
                   data = active_any_steps_M)
 
S_R(model_any_steps_M)

hr_model_any_steps_M <- give_HR(model_any_steps_M)%>% 
  mutate(case = 'Any Diagnosis'
          , Gender = 'Male') %>% 
  head(3)
hr_model_any_steps_M
```

### Femlae

```{r}
active_any_steps_F <- prepareData(Data_steps[Data_steps$Gender == 'F', ]
                          , history = any_diag_history
                          , occuranceDate = any_diag_date
                          , titleName = 'Any Diagnosis')

active_any_steps_F %>% select(activity_groups_num, newOccurance) %>% table()
 model_any_F_steps <- coxph(Surv(Duration, newOccurance) ~ activity_groups_num , data = active_any_steps_F)
S_R(model_any_steps)
give_HR(model_any_F_steps)

 model_any_steps_F <- coxph(Surv(Duration, newOccurance) ~ activity_groups_num +
                       age_at_index_date +
                       I(age_at_index_date * age_at_index_date) +
                       bmi_value +
                       log(imputed_pn_gp_visits_5yr) +
                       as.factor(SES) +
                       
                       HTN_diag +
                       hyperlipidemia_diag ,
                   data = active_any_steps_F)
 
S_R(model_any_steps_F)

hr_model_any_steps_F <- give_HR(model_any_steps_F)%>% 
  mutate(case = 'Any Diagnosis'
          , Gender = 'Female') %>% 
  head(3)
hr_model_any_steps_F
```

```{r}
M <- rbind(hr_model_any_steps_M, hr_model_any_steps_F, hr_model_cvd_steps_M, hr_model_cvd_steps_F , hr_model_Diabetes_steps_M, hr_model_Diabetes_steps_F, hr_model_any_steps_M, hr_model_any_steps_F, hr_model_CVA_steps_M, hr_model_CVA_steps_F) %>% 
  as.data.frame() %>% 
  mutate(case = as.factor(case)
         , term = as.factor(if_else(term == 'activity_groups_num1', 'Mild activity',
                                               if_else(term == 'activity_groups_num2', 'Moderate activity', 'High activity')))
         ) %>% 
  
  mutate(term = fct_relevel(term, c('High activity' , 'Moderate activity', 'Mild activity'))
         , case = fct_relevel(case, c('Diabetes', 'Stroke', 'CVD', 'Any Diagnosis')))  %>% 
  ggplot(aes(x = term, y = estimate, ymin = conf.low, ymax = conf.high, color =  Gender)) +
    geom_pointrange( fatten = 8, position = position_dodge(width = 0.5))+
    geom_hline(yintercept = 1, linetype = 2) +
    xlab('') +
  ylab("Hazard Ratio (95% Confidence Interval)") +
    facet_wrap(vars(case), strip.position = "left", nrow = 12, scales = "free_y") +
    coord_flip() +
    theme_classic() +
        theme(plot.title   = element_text(size = 20, face = "bold"),
          axis.text.y  = element_text(size = 18, face = "bold"),
          axis.ticks.y = element_blank(),
          axis.text.x  = element_text(size = 18, face = "bold"),
          axis.title   = element_text(size = 20, face = "bold"),
          strip.text.y = element_text(size = 18, hjust = 0.5, vjust = 1, angle = 180, face = "bold"),
          strip.text.x = element_text(size = 18, hjust = 0, vjust = 1, angle = 180, face = "bold"),
          legend.text   = element_text(size = 20) 
          , legend.title   = element_text(size = 20, face = "bold")
          , panel.border = element_rect(fill = NA)
        ) +
scale_color_manual(values=c("Black", "gray")) 

ggsave(plot = M, filename = 'gender.jpeg', dpi = 300, width = 40, height= 30, units = 'cm')

getwd()
```
