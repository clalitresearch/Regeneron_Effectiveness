---
title: "results_loss_to_FU_20231001"
author: "Tom Mushkat"
format: 
  html:
    code-fold: true
    code-summary: "Show the code"
editor: visual
---

```{r, warning=FALSE, message=FALSE}
pacman::p_load(tidyverse, cri.utils, skimr, gtsummary, ggfortify, ggplot2, gmodels, survival, grDevices, patchwork, smoothHR)
cri.utils::activate_keytab("tommuspim") # innovationprod
con <- cri.utils::cri_connection()
Data <- connect_table(con, "db.[scheme]", "table_name") %>%
  collect() %>%
      mutate(imputed_pn_gp_visits_5yr = ifelse(is.na(pn_gp_visits_5yr), round(median(pn_gp_visits_5yr, na.rm = T)), pn_gp_visits_5yr)
             , historic_DT2_diag = if_else(historic_DT2_diag == 2, 1, 0))
```

```{r, warning=FALSE, message=FALSE}


prepareData <- function(data
                        , history
                        , occuranceDate
                        , titleName){ 
      # 
    data <- data %>%
            filter({{history}} == 0)
      
     if(titleName %in% c('CVA')){ # 'CVD',

      data <- data %>%
            group_by(subclass) %>%
            mutate(lead_cvd = lead({{history}})
             , lag_cvd = lag({{history}})) %>%
             mutate(summing = if_else(({{history}} + lead_cvd == 1) | ({{history}} + lag_cvd == 1), 'remove', 'ok')) %>%
            filter(is.na(summing)) %>%
            ungroup()

    }

    
  data <- data %>%
       mutate(newOccurance = if_else(is.na({{occuranceDate}}), 0, 1)) %>% # 
       mutate(Duration     = if_else(newOccurance == 0, as.numeric(loss_to_FU - index_date), as.numeric({{occuranceDate}} - index_date))) %>%    
       mutate(newOccurance = if_else(Duration > 1000, 0, newOccurance)
              , Duration     = if_else(Duration > 1000, 1000, Duration))
   
    return(data)
}
 





S_R <- function(data) {

        test.ph = cox.zph(data)
      zph_df = as.data.frame(test.ph$table)
      ph_vars_list = zph_df %>% 
        mutate( var = dimnames(zph_df)[[1]]) %>% 
        filter(p < 0.05) %>% 
        filter(var != 'GLOBAL') %>% 
        pull  (var) 
# ph_vars_list = c('exposure') # TODO only the exposure is of intereset to us
for (v in ph_vars_list) {
          p = survminer::ggcoxzph(test.ph, var = c(v),
                       point.size = 0.5, point.alpha = 0.2, font.main = 12, caption = "Schoinfeld residulas")
          print(v)
          print(p)
          
  }
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
  
  
  
 # palette = c("gray70", "black"),  legend.labs = c("Non-app users", "App users"), 

  plot_ch_cvd <- survminer::ggsurvplot(the_model, fun = 'cumhaz', xlab = "", ylab = c(""), legend.title = "",legend.labs = c("Non-app users", "App users"), conf.int = TRUE, legend = "right", size = 0.05, censor.size = 0.01, xlim = c(0, 1000), font.main = c(20, "bold"), surv.scale = c("percent") ,  risk.table = TRUE) + guides(color = guide_legend(title = NULL))

plot_ch_cvd_data <- plot_ch_cvd$plot 

plt_1 <- plot_ch_cvd_data$data %>% 
  ggplot(aes(x = time, y = surv, group = strata, colour = strata)) + 
  theme_classic() +
  scale_color_manual(values = colors) +
  geom_errorbar(aes(ymin = lower, ymax = upper), color = "grey80") +
  geom_line() + 
  ggtitle(title) +
  theme(axis.text.x    = element_text(size = 20)
        , axis.text.y   = element_text(size = 20)
        , axis.title.y  = element_text(size = 20)
        , plot.title    = element_text(size = 20, face = 'bold')
        , legend.text   = element_text(size = 20)) + 
  xlab('') + ylab('') +
  theme(legend.position = legend.position) +
  scale_y_continuous(labels = scales::percent_format(accuracy = accuracy), limits = c(0, high_limit)) + 
  scale_x_continuous(breaks = c(0,  200, 400, 600, 800, 1000)) +
  guides(color = guide_legend(title = NULL))
  
  return(plt_1)
  
}

```

# General population

## CVD

```{r, warning=FALSE, message=FALSE}
activeCVD <- prepareData(Data
                         , history = historicCVD
                         , occuranceDate = CVDoccuranceDate
                         , titleName = 'CVD') %>% 
  mutate(Condition = as.factor(Condition)) %>% 
  mutate(Condition = fct_relevel(Condition, c('not_active', 'Active')))


activeCVD %>% select(Condition, newOccurance) %>% table()
an_cvd <- coxph(Surv(Duration, newOccurance) ~ Condition
      , data = activeCVD)
give_HR(an_cvd)


model_cvd <- coxph(Surv(Duration, newOccurance) ~ Condition +
                       age_at_index_date +
                       I(age_at_index_date * age_at_index_date) +
                       bmi_value +
                       log(imputed_pn_gp_visits_5yr) +
                       historic_DT2_diag +
                       HTN_diag +
                       hyperlipidemia_diag,
                   data = activeCVD)

S_R(model_cvd)
give_HR(model_cvd)

ch_cvd <- survival::survfit(survival::Surv(Duration, newOccurance) ~ Condition,
                               data = activeCVD, id = guid_tz)
plt_1 = create_figure(the_model = ch_cvd, accuracy = 1, high_limit = 0.05, title = "C. Cardiovascular Disease", colors = c("black", "gray"), legend.position = "none")

# https://stackoverflow.com/questions/59574997/how-to-remove-automated-strata-text-in-ggsurvplot-legend
```

### Diabetes

```{r, warning=FALSE, message=FALSE}

activeDiab <- prepareData(Data
                          , history = historic_DT2_diag
                          , occuranceDate = diab_occurance_date
                          , 'diabetes') %>% 
  mutate(Condition = fct_relevel(Condition, c('not_active', 'Active')))


activeDiab %>% select(Condition, newOccurance) %>% table()
an_diab <- coxph(Surv(Duration, newOccurance) ~ Condition
      , data = activeDiab)
give_HR(an_diab)

model_Diab <- coxph(Surv(Duration, newOccurance) ~ Condition +
                       age_at_index_date +
                       I(age_at_index_date * age_at_index_date) +
                       bmi_value +
                       log(imputed_pn_gp_visits_5yr) +
                       historicCVD +
                       HTN_diag +
                       hyperlipidemia_diag,
                   data = activeDiab)

S_R(model_Diab)
give_HR(model_Diab)


ch_surc_diab <- survival::survfit(survival::Surv(Duration, newOccurance) ~ Condition,
                               data = activeDiab, id = guid_tz)
plt_2 = create_figure(the_model = ch_surc_diab, accuracy = 0.1, high_limit = 0.015, title = "B. Diabetes Mellitus", colors = c("black", "gray"), legend.position = "none")
```

## Stroke

```{r, warning=FALSE, message=FALSE}

activeStroke <- prepareData(Data
                            , history = stroke_history
                            , occuranceDate = stroke_date
                            , titleName = 'CVA') %>% 
  mutate(Condition = fct_relevel(Condition, c('not_active', 'Active')))

activeStroke %>% select(Condition, newOccurance) %>% table()
an_stroke <- coxph(Surv(Duration, newOccurance) ~ Condition
      , data = activeStroke)
give_HR(an_stroke)

model_STKOKE <- coxph(Surv(Duration, newOccurance) ~ Condition +
                       age_at_index_date +
                       I(age_at_index_date * age_at_index_date) +
                       bmi_value +
                       log(imputed_pn_gp_visits_5yr) +
                       historic_DT2_diag +
                       HTN_diag +
                       AF_diag,
                   data = activeStroke)

S_R(model_STKOKE)
give_HR(model_STKOKE)


ch_stroke <- survival::survfit(survival::Surv(Duration, newOccurance) ~ Condition,
                               data = activeStroke, id = guid_tz)
plt_3 = create_figure(the_model = ch_stroke, accuracy = 0.1, high_limit = 0.015, title = "A. Stroke", colors = c("black", "gray"), legend.position = "none")

```

# Any diagnosis

```{r}

active_any <- prepareData(Data #[!c(Data$historic_DT2_diag == 1 & Data$historicCVD == 1), ]
                          , history = any_diag_history
                          , occuranceDate = any_diag_date
                          , titleName = 'Any Diagnosis') %>% 
  mutate(Condition = fct_relevel(Condition, c('not_active', 'Active')))

an_any <- coxph(Surv(Duration, newOccurance) ~ Condition
      , data = active_any)
give_HR(an_any)



model_any <- coxph(Surv(Duration, newOccurance) ~ Condition +
                     
                       age_at_index_date +
                       I(age_at_index_date * age_at_index_date) +
                       bmi_value +
                       log(imputed_pn_gp_visits_5yr) +
                       HTN_diag +
                       hyperlipidemia_diag,
                   data = active_any)

S_R(model_any)
give_HR(model_any)




 
ch_any <- survival::survfit(survival::Surv(Duration, newOccurance) ~ Condition,
                               data = active_any, id = guid_tz)
plt_4 = create_figure(the_model = ch_any, accuracy = 1, high_limit = 0.05, title = "D. Any Diagnosis", colors = c("black", "gray"), legend.position = "right")

```

```{r}
library(patchwork)
a = (plt_3 + plt_2  + plt_1 + plt_4)
a  
getwd()
ggsave(plot = last_plot(), filename = 'treatment_control.jpeg', dpi = 300, width = 40, height= 30, units = 'cm')
```
