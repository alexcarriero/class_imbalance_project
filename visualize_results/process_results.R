# ------------------------ set up-----------------------------------------------

algorithm   <- c("logistic_regression", 
                  "support_vector_machine", 
                  "random_forest", 
                  "xgboost", 
                  "rusboost", 
                  "easy_ensemble")

correction  <- c("control", 
                  "rus", 
                  "ros", 
                  "smo", 
                  "sen")

pairs <- 
  expand_grid(correction, algorithm) %>% 
  mutate(pair_id = c(1:30)) %>% 
  relocate(pair_id)

# ------------------ get summary -----------------------------------------------

get_summary <- function(scN, recalibrated = F){
  
  if(recalibrated == F){
    df <- 
      list.files(path = paste0("./simulation_results/raw_results/", scN,"/per_iter_results/"), pattern = ".rds") %>%
      paste0("./simulation_results/raw_results/", scN,"/per_iter_results/", .) %>%
      map(readRDS) %>% 
      bind_rows() %>% 
      rename("a_warn" = "a_warning")
  }
  
  if(recalibrated == T){
    df <- 
      list.files(path = paste0("./simulation_results/recalibrated_results/", scN,"/per_iter_results/"), pattern = ".rds") %>%
      paste0("./simulation_results/recalibrated_results/", scN,"/per_iter_results/", .) %>%
      map(readRDS) %>% 
      bind_rows()
  }
  
  tbl <- 
    as.data.frame(table(df$pair_id)) %>% 
    mutate(pair_id = as.numeric(as.character(Var1))) %>%
    dplyr::select(pair_id, Freq) %>%
    arrange(pair_id, decreasing = F)
  
  # problem_iterations <- 
  #   df %>% 
  #   mutate(
  #     iter = as.factor(iter), 
  #     pair_id = as.factor(pair_id), 
  #     scenario = as.numeric(scenario), 
  #     a_problem = as.numeric(a_problem),
  #     obs_ef = as.numeric(obs_ef), 
  #     new_ef = as.numeric(obs_ef)
  #   ) %>% 
  #   filter (a_problem == 1 | var == 0) 
  
  any_na_row <- df[!complete.cases(df), ]
  
  summary_stats <- 
    df %>% 
    group_by(pair_id) %>% 
    mutate(# scenario = as.numeric(scenario),
           int = as.numeric(int), 
           slp = as.numeric(slp)) %>% 
    summarise(
      # scenario = mean(scenario),
      auc_med  = median(auc, na.rm = T), 
      auc_sd   = sd  (auc, na.rm = T), 
      bri_med  = median(bri, na.rm = T), 
      bri_sd   = sd  (bri, na.rm = T), 
      int_med  = median(int, na.rm = T), 
      int_sd   = sd  (int, na.rm = T),
      slp_med  = median(slp, na.rm = T), 
      slp_sd   = sd  (slp, na.rm = T)
    ) 
  
  return(list("overall_dataframe" = df, "rows_per_pair" = tbl,  "summary"= summary_stats))
} 

# ------------------ re-calibration --------------------------------------------

recalibrate_everything <- function(i, scN){
  
  #  import probabilities from simulation
  df <- 
    i %>% 
    paste0("./simulation_results/raw_results/", scN, "/predictions/", scN, "_", . , ".rds") %>% 
    map(readRDS) 
  
  # save re-calibrated probabilities
  r_df <- 
    df %>%
    recalibrate_one_iteration()
  
  saveRDS(r_df, paste0("./simulation_results/recalibrated_results/", scN, "/predictions/", scN, "_", i, ".rds"))
  
  # save re-calibrated per-iteration statistics
  r_stats <- 
    r_df %>% 
    get_stats_recalibrated() 
  
  saveRDS(r_stats, paste0("./simulation_results/recalibrated_results/", scN, "/per_iter_results/", scN, "_", i, ".rds"))
  
  # clean environment
  rm(df, r_df, r_stats)
}

recalibrate_one_iteration <- function(dataframe){
  
  # group data frame by pair_id and re-calibrate
  dataframe %>% 
    as.data.frame() %>% 
    mutate(
      pair_id = as.factor(pair_id), 
      pred = as.numeric(pred), 
      class = as.factor(class)
    ) %>%
    group_by(pair_id) %>%
    group_modify(~recalibrate(.)) 

}

recalibrate <- function(dataframe){
  
  # retrieve re-calibrated probabilities
  re_calibrate <- 
    recalibrated_probabilities(dataframe$pred, dataframe$class)
  
  # save warning of from glm
  r_warn <<- ifelse(is.null(re_calibrate$warning) == TRUE, "0", as.character(re_calibrate$warning))
  
  # update probabilities to re-calibrated probabilities
  dataframe <- 
    dataframe %>% 
    mutate(
      pred = re_calibrate$pred,
      recalibrated = 1
    ) %>% 
    mutate(class = as.numeric(class) -1)
  
  return(dataframe)
}

recalibrated_probabilities <- function(probs, class){
  
  # computing re-calibrated probabilities
  if(any(is.na(probs)) == TRUE){
    return(list("pred" = NA, "warning" = as.character("there were NAs")))
  }
  if(any(is.na(probs)) == FALSE & any(probs < 0) == TRUE){
    probs[probs < 0 ] <- 1e-10
  }
  if(any(is.na(probs)) == FALSE & any(probs > 1) == TRUE){
    probs[probs > 1 ] <- 1 - 1e-10
  }
  if(sum(probs == 1) != 0){
    probs[probs == 1] <- 1 - 1e-10
  }
  if(sum(probs == 0) != 0){
    probs[probs == 0] <- 1e-10
  }
  
  # ensure correct factor levels 
  class  <- factor(class, levels=c(0,1), ordered = TRUE)
  
  # model 
  a <- tryCatch.W.E(
    expr    = { 
      mod  <- glm(class ~ 1, offset = log(probs/(1-probs)), family = "binomial")
      pred <- predict(mod, type = "response")
    }
  )
  
  return(list("pred" = a$value, "warning" = a$warning))
}

# ------------------- calibration plots ----------------------------------------

reorder_and_rename <- function(dataframe){
  dataframe %>% 
  mutate(
      algorithm  = factor(algorithm,
                          levels = c(
                            'logistic_regression',
                            'support_vector_machine',
                            'random_forest',
                            'xgboost',
                            'rusboost',
                            'easy_ensemble')),
      correction = factor(correction,
                          levels = c(
                            'control', 'rus', 'ros', 'smo', 'sen'))) %>%
    mutate(
      algorithm = recode(algorithm,
                         "logistic_regression"    = "Logistic Regression",
                         "support_vector_machine" = "Support Vector Machine",
                         "random_forest"          = "Random Forest",
                         "xgboost"                = "XGBoost" ,
                         "rusboost"               = "RUSBoost",
                         "easy_ensemble"          = "EasyEnsemble"),
      correction = recode(correction,
                          "control" = "Control",
                          "rus"    = "RUS",
                          "ros"    = "ROS",
                          "smo"    = "SMOTE",
                          "sen"    = "SMOTE-ENN"))
}

ready2plot <- function(dataframe){
  dataframe %>% 
  # filter(problem != 1) %>% 
  drop_na() %>% 
  mutate(
      pred  = as.numeric(pred), 
      class = as.numeric(class)) %>%
  reorder_and_rename()
}

loess_regression<- function(dataframe){
  
  x <- c(seq(0, 1, length.out = 200))
  iter <- as.numeric(dataframe$iter[1])
  
  if(any(is.na(dataframe$pred)) == TRUE){
    y <- rep(NA, 200)
  } else{
    
    a <- 
    tryCatch.W.E(
      
      expr  = {
        mod  <- dataframe %>% loess(formula = class ~ pred)
        y    <- predict(mod, newdata = x)}
  
    )
    
    if(is.vector(a$value) == FALSE){
      y <- rep(NA, 200)
    }
  }
  
  out  <- tibble(iter = rep(iter, 200), x  = x, y  = y)
  return(out)
}

save_plot_coords <- function(i, scN, recalibrated = F){
  
  if(recalibrated == T){
    df <- 
      i %>% 
      paste0("./simulation_results/recalibrated_results/", scN, "/predictions/", scN, "_", . , ".rds") %>% 
      map(readRDS) %>% 
      as.data.frame() %>%
      ready2plot() %>% 
      group_by(pair_id) %>% 
      group_modify(~loess_regression(.)) %>% 
      merge(pairs, by = "pair_id") %>% 
      reorder_and_rename() %>% 
      mutate(iter = as.factor(iter))
    
      saveRDS(df, paste0("./processed_results/plot_coords/recalibrated_results/", scN, "/", scN, "_",  i, "_plot_coords.rds"))
  }
  
  if(recalibrated == F){
    df <- 
      i %>% 
      paste0("./simulation_results/raw_results/", scN, "/predictions/", scN, "_", . , ".rds") %>% 
      map(readRDS) %>% 
      as.data.frame() %>%
      ready2plot() %>% 
      group_by(pair_id) %>% 
      group_modify(~loess_regression(.)) %>% 
      merge(pairs, by = "pair_id") %>% 
      reorder_and_rename() %>% 
      mutate(iter = as.factor(iter))
    
    saveRDS(df, paste0("./processed_results/plot_coords/raw_results/", scN, "/", scN, "_",  i, "_plot_coords.rds"))
  }
    
}

pre_plot <- function(df){
  
  # filter for appropriate values 
  df <- df %>% filter(y < 1 & y > 0) 
  
  # determine which groups don't have enough obs to make a plot
  df <- df %>% mutate(iter = as.numeric(as.character(iter)))
  
  check   <- df %>% count(iter, pair_id)
  problem <- which(check$n < 4)
  
  if(length(problem) > 0){
    for (i in 1:length(problem)){
      
      key   <- check[problem[i], ]
      index <- which(df$iter == key$iter & df$pair_id == key$pair_id)
      df    <- df[-index, ]
      
    }
  }
  
  return(df)
}

plot_from_coords <- function(df){
  p <-
    ggplot() +
    geom_abline(slope = 1, intercept = 0, linewidth = 1) +
    geom_smooth(data = df, aes(x = x, y = y, group = iter), 
                method = "loess",
                formula = y ~ x, 
                se = F,
                linewidth = 0.01, 
                alpha = 0.01,
                color = "#00008B")+ 
    theme_minimal() +
    xlab(" ") + 
    ylab(" ") +
    xlim(0,1) + 
    ylim(0,1) + 
    facet_grid(algorithm ~ correction) +
    theme(legend.position = "none") 
  
  return(p) 
}


# ----------------performance metric plots -------------------------------------

tidy_for_pm_plts <- function(dataframe){
  dataframe %>% 
    mutate(
      algorithm  = factor(algorithm, 
                          levels = c(
                            'easy_ensemble',
                            'rusboost',
                            'xgboost',
                            'random_forest',
                            'support_vector_machine',
                            'logistic_regression')),
      correction = factor(correction, 
                          levels = c(
                            'control', 'rus', 'ros', 'smo', 'sen'))) %>% 
    mutate(
      algorithm = recode(algorithm, 
                         "logistic_regression"    = "Logistic Regression",
                         "support_vector_machine" = "Support Vector Machine",
                         "random_forest"          = "Random Forest",
                         "xgboost"                = "XGBoost" ,
                         "rusboost"               = "RUSBoost",
                         "easy_ensemble"          = "EasyEnsemble"), 
      correction = recode(correction,
                          "control" = "Control", 
                          "rus"    = "RUS", 
                          "ros"    = "ROS", 
                          "smo"    = "SMOTE", 
                          "sen"    = "SMOTE-ENN")) %>% 
    mutate(int = as.numeric(int), 
           slp = as.numeric(slp))
}

pm_plot <- function(dataframe, method, xmin, xmax){
  
  if(method == "auc"){
    title <- "Concordance Statistic"
    ref   <- 0.85
  }
  if(method == "slp"){
    title <- "Calibration Slope"
    ref   <- 1
  }
  if(method == "int"){
    title <- "Calibration Intercept"
    ref <- 0
  }
  if(method == "bri"){
    title <- "Scaled Brier Score"
    ref <- 0
  }
  
  dataframe %>% 
    tidy_for_pm_plts() %>%
    ggplot() + 
    geom_hline(yintercept = ref, color = "black", linewidth = 1) + 
    geom_violin(aes_string(x = "algorithm", y = method),
                fill  = "#b3cd41", 
                color = "grey",
                draw_quantiles = 0.5,
                linewidth = 0.25,
                alpha = 0.5) + 
    # scale_fill_viridis(discrete = T, direction = -1, option = "D") + 
    theme_minimal() +
    xlab(" ") + 
    ylab(" ") +
    ylim(xmin, xmax) +
    coord_flip() + 
    facet_wrap(~correction, nrow = 1) + 
    ggtitle(title)
}