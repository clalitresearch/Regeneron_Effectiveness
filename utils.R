#############################################################################################
show_worst_matches = function(matched_df, var, simplified_var, quantile = 0.85, nrows = 5) {
#############################################################################################
  diffs = matched_df %>%
    group_by(subclass) %>%
    mutate(diff = max( {{ var }}) - min( {{ var }})) %>% 
    mutate(diff_simp =  max( {{ simplified_var }}) - min( {{ simplified_var }})) %>% 
    ungroup() %>%
    arrange(desc(diff_simp),subclass) %>%
    select(subclass,{{var}}, {{simplified_var}}, diff, diff_simp)  %>%
    head(nrows*2)
  
  colname            = diffs %>% select({{var}}) %>%  colnames()
  simplified_colname = diffs %>% select({{simplified_var}}) %>%  colnames()
  
  print(glue("{100*quantile}% tile of {colname} diffs = {quantile(diffs$diff, quantile, na.rm = T)}")) 
  print(glue("{100*quantile}% tile of {simplified_colname}  diffs = {quantile(diffs$diff_simp, quantile)}"))
  return(diffs)
}

#############################################################################################
plot_var = function(matched_df, var){
#############################################################################################
  var_name = matched_df %>% select({{var}}) %>%  colnames()
  ggplot(data =  matched_df, mapping = aes(x = {{var}}, y = ..density.., color = arm)) +
    scale_color_manual(values = c("red" , "green")) +
    theme_bw() +
    ggtitle(glue("Distribution of {var_name} on the matched data")) +
    geom_freqpoly(bins =10 ,position = position_jitter(width  = 0.035)) 
}



#############################################################################################
print_total = function(.data, msg) {
#############################################################################################
  total = .data %>%
    dplyr::summarise(total = n()) %>% 
    pull(total)
  treated = .data %>%
    dplyr::summarise(treated = sum(if_else(is.na(regeneron_date), 0, 1))) %>% 
    pull(treated)
 untreated = .data %>%
   dplyr::summarise(untreated = sum(if_else(is.na(regeneron_date), 1, 0))) %>% 
    pull(untreated)
  message(msg,": total:",total," treated:", treated, " untreated:", untreated);
  .data
}



####################################################
match_it = function(df, formula, exact_formula, num_controls, max_controls = NULL, min_controls = NULL, distance_method = "mahalanobis") {
####################################################
  stopifnot(distance_method %in% c("mahalanobis", "ps"))
  if (distance_method == "mahalanobis")
    return (match_it_mahalanobis(df, formula, exact_formula, num_controls, max_controls, min_controls))
  if (distance_method == "ps")
    return (match_it_ps(df, formula, exact_formula, num_controls, max_controls, min_controls))
}


####################################################
match_it_ps = function(df, formula, exact_formula, num_controls, max_controls = NULL, min_controls = NULL) {
####################################################
  print(glue("propensity score matching "))
  print(glue("num_controls={num_controls} , max_controls = {max_controls}, min_controls = {min_controls}"))
  mm  = matchit(formula,
                data = df,
                method = "optimal",
                distance = "glm",
                link = "logit",
                verbose = T,
                max.controls = max_controls ,
                ratio = num_controls,
                min.controls = min_controls,
                k2k = T)
  print(glue("will call match.data on mm "))
  final_df = match.data(mm,data = df) 
  print ("created final_df from mm")
  final_df %>% 
    mutate(arm = ifelse(exposure_cat==0, "control", "intervention"))
}


####################################################
match_it_mahalanobis = function(df, formula, exact_formula, num_controls, max_controls = NULL, min_controls = NULL) {
####################################################
  print("matching with mahalanobis distance")
  mm  = matchit(formula,
              data = df,
              exact = exact_formula,
              method = "optimal",
              distance = "mahalanobis",
              verbose = T,
              max.controls = max_controls ,
              ratio = num_controls,
              min.controls = min_controls,
              k2k = T)

  final_df = match.data(mm) %>% 
    mutate(arm = ifelse(exposure_cat==0, "control", "intervention"))
  return(final_df)
}



####################################################
get_pairs_by_immortal_time_bias_status = function(df, itb_status = T) {
####################################################

  treated = df %>% 
    mutate (immortal_time   = as.numeric(difftime(regeneron_date, index_date,   units = "days"))) %>%
    select(subclass, guid_tz, arm, index_date, immortal_time) %>% 
    filter(arm == "intervention")
  
  controls = df %>%
    filter(arm == "control") %>% 
    mutate(days_to_earliest_outcome_control = pmin(days_to_hosp, days_to_severe, days_to_death, days_to_end_fwp, na.rm = T)) %>% 
    select(subclass, guid_tz, arm, index_date, days_to_earliest_outcome_control)
     
  joined_df = controls %>% 
    inner_join(treated, by=c('subclass')) %>% 
    mutate(immortal_bias = immortal_time  > days_to_earliest_outcome_control) %>% 
    filter (immortal_bias   == itb_status) %>% 
    rename(tz_control = guid_tz.x) %>% 
    rename(tz_treated = guid_tz.y) %>% 
    select(tz_control, tz_treated,immortal_time, days_to_earliest_outcome_control)
   
  return (joined_df)
}

####################################################
match_multiple_controls_without_immortal_time_bias = function(df, formula, exact_formula, num_controls=1, distance_method = "mahalanobis") {
####################################################
# this function does repeated iterations of  optimal matching but verifying that matched pairs are not in a state
# of immortal time for the treated twin, if it does , the matched control is discarded from being matched in the future

## create an empty data frame with the same columns as df
  correct_pairs = df[FALSE,]
  subclass_dedup_coeffecient = 1
  num_treated = df %>%  filter(exposure == 1) %>%  nrow()
  num_pairs = (num_treated * num_controls) + 1
  print(glue("starting the matching loop with formula {Reduce(paste, deparse(formula))}"))
  while (T) {
    print("---- start of iteration ----")
    # we need to manually make sure that all the subclasses created by matchit are unique
    # hence in iteration K we will add to each subclass id , K * num_controls times the total number of possible pairs 
    subclass_dedup_coeffecient = subclass_dedup_coeffecient + num_pairs
    
    print (glue("Matching on a df with {nrow(df)} rows ..."))
    matched_df = match_it(df, formula, exact_formula, num_controls = num_controls-0.000001, min_controls = num_controls-1,
                          max_controls = num_controls, distance_method = distance_method)
    print (glue("Finished Matching. returned df has {nrow(matched_df)} rows ..."))
    matched_df = matched_df%>% 
      mutate(subclass = factor(subclass_dedup_coeffecient + as.numeric(subclass)))
    print (glue("Finished Matching. returned df has {nrow(matched_df)} rows ..."))
    
    valid_pairs = get_pairs_by_immortal_time_bias_status(matched_df, FALSE)
    invalid_pairs = get_pairs_by_immortal_time_bias_status(matched_df, TRUE)
    print (glue("Threre are {nrow(valid_pairs)} pairs without immortal time bias"))
    print (glue("Threre are {nrow(invalid_pairs)} pairs with immortal time bias"))
    
    treated_correctly_matched = matched_df %>% 
      filter(arm == "intervention") %>% 
      semi_join(valid_pairs, by=c("guid_tz" = "tz_treated")) #keep all treated obs. in matched_df that have a match in valid_pairs
    
    controls_correctly_matched = matched_df %>% 
      filter(arm == "control") %>% 
      semi_join(valid_pairs, by=c("guid_tz" = "tz_control"))  #keep all control obs. in matched_df that have a match in valid_pairs
    
    correct_pairs = bind_rows(correct_pairs, treated_correctly_matched, controls_correctly_matched)
    
    n_distinct_treated_with_a_match = correct_pairs %>% 
      filter(arm == "intervention") %>% 
      distinct(guid_tz) %>% 
      nrow()
    
    print(glue("n_distinct_treated_with_a_match = {n_distinct_treated_with_a_match}, num_treated = {num_treated}"))
    # if every treated subject has at least one valid control, we are ok, can exit the loop.
    if (num_treated - n_distinct_treated_with_a_match == 0 ) 
      return (correct_pairs)

    updated_controls =  df %>% 
      filter(exposure == 0) %>% 
      anti_join(correct_pairs, by="guid_tz" )  %>%               #drop all  valid controls from previous iterations
      anti_join(invalid_pairs, by=c("guid_tz" = "tz_control") )  #drop all invalid controls from this iteration
    
    needs_rematch = df %>% 
      filter(exposure == 1) %>% 
      semi_join(invalid_pairs, by=c("guid_tz" = "tz_treated")) 
     
    df = bind_rows(needs_rematch, updated_controls) 
    
    num_treated_unmatched = nrow(needs_rematch) 
    print(glue("preparing for another iteration on {num_treated_unmatched} treated yet unmatched... and {nrow(updated_controls)} controls "))
     
  }
  return (correct_pairs)
  
}
