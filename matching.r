---
title: "Matching"
author: "Tom Mushkat"
format: html
editor: visual
---

```{r}
pacman::p_load(tidyverse, cri.utils, MatchIt, skimr, optmatch, gtsummary, cobalt, parallel)

activate_keytab("username")
system2("kinit" ,args=("-k -t "), stdout = TRUE, stderr = TRUE)
con <- cri_connection()
raw_data_con <- connect_table(con, "db.[scheme]", "table_name")
raw_data <- raw_data_con %>% dplyr::collect() %>%
  mutate(Condition = if_else(Condition == 'Active', 1, 0))

# getwd()
write_rds(raw_data, 'file')
raw_data <- read_rds('file')

noNA_data <- raw_data %>% 
  # drop_na() %>% 
  mutate(Sector = as.character(Sector),
         corona_COM_CAT = as.character(corona_COM_CAT),
         kod_machoz = as.character(kod_machoz),
         SUG_YISHUV = as.character(SUG_YISHUV),
         monthYear = as.character(monthYear)
         )

```

# Matching functions- goal: to prepar a new data set for the next matching iteration.

## Those two functions take the data frame than been used to the matching process and remove IDs that were already matched. For patients

### returnData prepare the data set and been used inside returnLeftIDs which return the IDs

```{r}

returnData <- function(matchingOutcome, dataFor_matchingOutcome){
  
    # parameters fo test
      # matchingOutcome <- Round_4
      # dataFor_matchingOutcome <- newData_3
  
  # The function take the result of the matching and the data that was created it and return as a data frame only subclass that were matched

  matchingData <- match.data(matchingOutcome, "all") %>% as.data.frame()
  set.seed(123)
  selectedSubclass <- matchingData %>% 
    filter(Condition == 0) %>%
    group_by(teudat_zehut) %>% 
    mutate(RN = row_number()) %>%
    filter(RN == 1) %>% 
    pull(subclass) 
    
  selectedSubclass <- unique(selectedSubclass)
  
  selectedIDs <- matchingData %>% 
    filter(subclass %in% selectedSubclass)
  
  return(selectedIDs)
  
}


returnLeftIDs <- function(matchingOutcome, dataFor_matchingOutcome){
  
  # parameters fo test
    # matchingOutcome <- Round_4
    # dataFor_matchingOutcome <- newData_3
  
  # This function return the data without the patients who already were matched using the returnData function

  selectedIDs <- returnData(matchingOutcome)
  
  data_1 <- 
    dataFor_matchingOutcome %>% 
    filter(!teudat_zehut %in% selectedIDs$teudat_zehut) 
  
  return(data_1)
  
}




```


### Matching process using 5 rounds

```{r echo=FALSE}

FORMULA_ALL <- Condition ~ Gender + Sector + corona_COM_CAT + monthYear + BMI_CAT + SUG_YISHUV + kod_machoz + HbA1c_Base_HbA1c_CAT + historic_diab_complication + historicCVD + historic_DT2_diag + age_at_index_date

FORMULA_EXACT <- Condition ~ Gender + Sector + corona_COM_CAT + monthYear + BMI_CAT + SUG_YISHUV + kod_machoz + HbA1c_Base_HbA1c_CAT + historic_diab_complication + historicCVD + historic_DT2_diag

Limits <- c(age_at_index_date = 5)

set.seed(42)
Round_1  = matchit(formula = FORMULA_ALL,
                   data    = noNA_data,
                   exact   =  FORMULA_EXACT,
                   caliper = Limits,
                   std.caliper = FALSE,
                   # na.action = "na.pass",
                   method  = "nearest", distance = "mahalanobis", verbose = F, replace = F, k2k = T)

saveRDS(Round_1, '/workspace//Clalit_Active_General_pop_nov22//Clalit_Active_General_pop_nov22//Round_1.rds')
newData_1 <- returnLeftIDs(matchingOutcome = Round_1, dataFor_matchingOutcome = noNA_data)
saveRDS(newData_1, '/workspace//Clalit_Active_General_pop_nov22//Clalit_Active_General_pop_nov22//newData_1.rds')
noNA_data <- NULL

gc()
set.seed(345)
Round_2  = matchit(formula = FORMULA_ALL,
                   data    = newData_1,
                   exact   =  FORMULA_EXACT,
                   caliper = Limits,
                   std.caliper = FALSE,
                   # na.action = "na.pass",
                   method = "nearest", distance = "mahalanobis", verbose = F, replace = F, k2k = T)

saveRDS(Round_2, '/workspace//Clalit_Active_General_pop_nov22//Clalit_Active_General_pop_nov22//Round_2.rds')
newData_2 <- returnLeftIDs(matchingOutcome = Round_2, dataFor_matchingOutcome = newData_1)
saveRDS(newData_2, '/workspace//Clalit_Active_General_pop_nov22//Clalit_Active_General_pop_nov22//newData_2.rds')
newData_1 <- NULL
Round_2 <- NULL


```

```{r}
getwd()
gc()
Round_1 <- read_rds('/workspace/Clalit_Active_General_pop_nov22/Clalit_Active_General_pop_nov22/Round_1.rds')
Round_2 <- read_rds('/workspace/Clalit_Active_General_pop_nov22/Clalit_Active_General_pop_nov22/Round_2.rds')
newData_1 <- read_rds('/workspace//Clalit_Active_General_pop_nov22//Clalit_Active_General_pop_nov22//newData_1.rds')

plot_df_1 <- returnData(Round_1, noNA_data) %>%
  mutate(Condition = ifelse(Condition == 0, "not_active", "Active")) %>%
  mutate(subclass = paste0(subclass, '_1'))

gc()
plot_df_2 <- returnData(Round_2, newData_1) %>%
  mutate(Condition = ifelse(Condition == 0, "not_active", "Active")) %>%
  mutate(subclass = paste0(subclass, '_2'))

full_data <- rbind(plot_df_1, plot_df_2)
saveRDS(full_data, '/workspace//Clalit_Active_General_pop_nov22//Clalit_Active_General_pop_nov22//full_data.rds')
write_sql(con, table = full_data, name = 'Clalit_Active_General_pop_nov22', schema = NULL, overwrite = T, append = F)
```


### Covariate balance after matching

```{r love_plot,  echo=F, warning=F, message=F}
LOVE_1 <-
love.plot(Round_1,
          drop.distance = TRUE,
          var.order = "alphabetical",
          abs = T,
          line = F,
          thresholds = c(m = .1))

# LOVE_1
saveRDS(LOVE_1, 'LOVE_1.rds')
```




