arcsine_two_way <- function(df, expression) {
  main_analysis_acc_aggregate <-
    df %>% 
    group_by(participant_id, is_congruent, is_previous_congruent) %>% 
    summarise(acc_conditional_mean = mean(is_correct, na.rm = T)) %>% 
    ungroup() %>% 
    mutate(acc_conditional_mean = !! expression,
           is_congruent = as_factor(is_congruent),
           is_previous_congruent = as_factor(is_previous_congruent),
           participant_id = as_factor(participant_id))
  
  anova_acc <-
    ezANOVA(data = main_analysis_acc_aggregate, dv = acc_conditional_mean, wid = participant_id, within = .(is_congruent, is_previous_congruent), type = 3)
  
  anova_acc <-
    as_tibble(anova_acc$ANOVA)
  
  anova_acc_f_value <-
    anova_acc %>% 
    filter(Effect %in% c("is_congruent",
                         "is_previous_congruent",
                         "is_congruent:is_previous_congruent")) %>% 
    transmute(term = case_when(Effect == "is_congruent" ~ "main_effect_congruent",
                               Effect == "is_previous_congruent" ~ "main_effect_previous_congruent",
                               Effect == "is_congruent:is_previous_congruent" ~ "interaction_effect"),
              f_value = F,
              df = DFd)
  
  main_analysis_acc_participant_level_raw_effect <-
    main_analysis_acc_aggregate %>% 
    ungroup() %>% 
    mutate(condition = case_when(is_previous_congruent == 0L & is_congruent == 0L ~ "ii",
                                 is_previous_congruent == 0L & is_congruent == 1L ~ "ic",
                                 is_previous_congruent == 1L & is_congruent == 0L ~ "ci",
                                 is_previous_congruent == 1L & is_congruent == 1L ~ "cc",
                                 TRUE ~ NA_character_)) %>% 
    select(-is_previous_congruent, -is_congruent) %>% 
    spread(key = condition, value = acc_conditional_mean) %>% 
    mutate(congruency_effect =  ((cc + ic) / 2) - ((ci + ii) / 2),
           previous_congruency_effect = ((ci + cc) / 2) - ((ii + ic) / 2),
           interaction_effect = (cc - ci) - (ic - ii))
  
  main_analysis_acc_raw_effect <-
    main_analysis_acc_participant_level_raw_effect %>% 
    summarise(main_effect_congruent = mean(congruency_effect, na.rm = T),
              main_effect_previous_congruent = mean(previous_congruency_effect, na.rm = T),
              interaction_effect = mean(interaction_effect, na.rm = T)) %>%  
    gather(key = "term", value = "raw_effect") %>%
    inner_join(., anova_acc_f_value, by = "term") %>% 
    mutate(se = abs(raw_effect / sqrt(f_value)))
  
  interaction_se <- 
    main_analysis_acc_raw_effect %>% 
    filter(term == "interaction_effect") %>% 
    pull(se)
  
  interaction_effect <- 
    main_analysis_acc_raw_effect %>% 
    filter(term == "interaction_effect") %>% 
    pull(raw_effect)
  
  interaction_df <- 
    main_analysis_acc_raw_effect %>% 
    filter(term == "interaction_effect") %>% 
    pull(df)
  
  Bf(sd = interaction_se,
     obtained = interaction_effect,
     dfdata = interaction_df,
     meanoftheory = 0,
     sdtheory = .03,
     dftheory = 10^10,
     tail = 1)
}

arcsine_three_way <- function(df, expression) {
  main_analysis_acc_aggregate <-
    df %>% 
    group_by(participant_id, reversal_group, is_congruent, is_previous_congruent) %>% 
    summarise(acc_conditional_mean = mean(is_correct, na.rm = T)) %>% 
    ungroup() %>% 
    mutate(acc_conditional_mean = !! expression,
           is_congruent = as_factor(is_congruent),
           is_previous_congruent = as_factor(is_previous_congruent),
           participant_id = as_factor(participant_id),
           reversal_group = as_factor(reversal_group))
  
  anova_acc <-
    ezANOVA(data = main_analysis_acc_aggregate, dv = acc_conditional_mean, wid = participant_id, within = .(is_congruent, is_previous_congruent), between = reversal_group, type = 3)
  
  anova_acc <-
    as_tibble(anova_acc$ANOVA)
  
  anova_acc_f_value <-
    anova_acc %>% 
    filter(Effect %in% c("is_congruent",
                         "is_previous_congruent",
                         "reversal_group:is_congruent:is_previous_congruent")) %>% 
    transmute(term = case_when(Effect == "is_congruent" ~ "main_effect_congruent",
                               Effect == "is_previous_congruent" ~ "main_effect_previous_congruent",
                               Effect == "reversal_group:is_congruent:is_previous_congruent" ~ "interaction_effect"),
              f_value = F,
              df = DFd)
  
  main_analysis_acc_participant_level_raw_effect <-
    main_analysis_acc_aggregate %>% 
    ungroup() %>% 
    mutate(condition = case_when(is_previous_congruent == 0L & is_congruent == 0L ~ "ii",
                                 is_previous_congruent == 0L & is_congruent == 1L ~ "ic",
                                 is_previous_congruent == 1L & is_congruent == 0L ~ "ci",
                                 is_previous_congruent == 1L & is_congruent == 1L ~ "cc",
                                 TRUE ~ NA_character_)) %>% 
    select(-is_previous_congruent, -is_congruent) %>% 
    spread(key = condition, value = acc_conditional_mean) %>% 
    mutate(congruency_effect = ((cc + ic) / 2) - ((ci + ii) / 2),
           previous_congruency_effect = ((ci + cc) / 2) - ((ii + ic) / 2),
           interaction_effect = (cc - ci) - (ic - ii),
           current_incongruent = ii - ci)
  
  main_analysis_acc_raw_effect <-
    main_analysis_acc_participant_level_raw_effect %>% 
    group_by(reversal_group) %>% 
    summarise(main_effect_congruent = mean(congruency_effect, na.rm = T),
              main_effect_previous_congruent = mean(previous_congruency_effect, na.rm = T),
              interaction_effect = mean(interaction_effect, na.rm = T)) %>%  
    gather(key = "term", value = "raw_effect", -reversal_group) %>%
    spread(key = "reversal_group", value = "raw_effect") %>% 
    ungroup() %>% 
    mutate(raw_effect = high - low) %>% 
    select(-high,
           -low) %>% 
    inner_join(., anova_acc_f_value, by = "term") %>% 
    mutate(se = abs(raw_effect / sqrt(f_value)))
  
  acc_sd <-
    main_analysis_acc_raw_effect %>% 
    filter(term == "interaction_effect") %>% 
    pull(se)
  
  acc_obtained <-
    main_analysis_acc_raw_effect %>% 
    filter(term == "interaction_effect") %>% 
    pull(raw_effect)
  
  acc_df <-
    main_analysis_acc_raw_effect %>% 
    filter(term == "interaction_effect") %>% 
    pull(df)
  
  Bf(sd = acc_sd,
     obtained = acc_obtained,
     dfdata = acc_df,
     meanoftheory = 0,
     sdtheory = .03,
     dftheory = 10^10,
     tail = 1)
}

