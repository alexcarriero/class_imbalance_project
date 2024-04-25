# ------------------------ set up-----------------------------------------------

# Create a data frame which serves as a key for each prediction model 
# This data frame has 30 rows and 3 columns 
#
# Each row represents a unique prediction model 
# Columns are: 
#
# algorithm - character specifying the algorithm name 
# correction - character specifying the correction name 
# pair_id - an integer between 1 and 30 that is unique for each prediction model


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

get_brier_score <- function(dataframe){
  
  # input : a dataframe
  # it is necessary for the dataframe to have the columns "pred" and "class"
  #
  # pred  : predicted risks for observations in one data set
  # class : true observed class for observations in one data set
  #
  # it returns the brier score, as a single integer saved in a tibble. 
  
  tibble(
    brier_score = brier_score(dataframe$pred, dataframe$class)
  )
}

save_brier_scores <- function(i){
  
  # input: a single integer (i) 
  # this integer represents the simulation scenario 
  # 
  # for each iteration in scenario i 
  # the function imports the raw predicted risks from the appropriate folder, 
  # computes the Brier score based on the raw_predicted risks and 
  # saves the value of the Brier score for all prediction models in the iteration
  # 
  # Results are then merged, for all models in all iterations into one dataframe 
  # and saved in a new file called scX_brier_scores.rds. 
  
  scN <- paste0("sc", i)
  
  out = tibble()
  
  for(i in 1:2000){
    
    df <- 
      i %>% 
      paste0("./simulation_results/raw_results/", scN, "/predictions/", scN, "_", . , ".rds") %>% 
      map(readRDS) %>% 
      as.data.frame()
    
    temp <- 
      df %>% 
      mutate(
        pair_id = as.factor(pair_id), 
        pred = as.numeric(pred), 
        class = as.factor(class)
      ) %>%
      group_by(pair_id) %>%
      group_modify(~get_brier_score(.)) %>% 
      mutate(iter = i) %>% 
      bind_rows(out)
    
    out <- temp
  }
  
  saveRDS(out, paste0("./simulation_results/raw_results/", scN,"/", scN, "_brier_scores.rds"))
}

get_summary <- function(i, recalibrated = F){
  
  # input: 
  # i - integer representing the simulation scenario 
  # recalibrated - boolean representing specification of raw results (F) or 
  #                recalibrated results (T)
  # 
  # 
  # For simulation scenario i
  # this function returns a list with 3 objects: 
  #
  # overall_dataframe: 
  # 60,000 rows, one for every iteration (2000) x prediction model (30)
  # this is equivalent to binding the rows of all files in  per_iter_results for scenario i 
  #
  # rows_per_pair: 
  # 30 rows, one for each prediction model
  # a column indicating how many observations there are for each model 
  # this is a check to ensure there are 2000 
  #
  # summary: 
  # 30 rows, one for each prediction model 
  # summary statistics as columns representing the median value of each 
  # empirical performance metrics over the 2000 simulation iterations 
  # and the corresponding MCMC errors. 
  #
  # auc = concordance statistic 
  # bri = Brier score 
  # int = calibration intercept 
  # slp = calibration slope
  
  
  scN <- paste0("sc", i)
  print(scN)
  if(recalibrated == F){
    
    brier_scores <- readRDS(paste0("./simulation_results/raw_results/", scN,"/", scN, "_brier_scores.rds"))
    
    df <- 
      list.files(path = paste0("./simulation_results/raw_results/", scN,"/per_iter_results/"), pattern = ".rds") %>%
      paste0("./simulation_results/raw_results/", scN,"/per_iter_results/", .) %>%
      map(readRDS) %>% 
      bind_rows() %>%
      rename("a_warn" = "a_warning",
             "s_bri" =  "bri") %>%
      merge(brier_scores, by = c("pair_id", "iter")) %>%
      rename("bri" = "brier_score")
  }
  
  if(recalibrated == T){
    df <- 
      list.files(path = paste0("./simulation_results/recalibrated_results/", scN,"/per_iter_results/"), pattern = ".rds") %>%
      paste0("./simulation_results/recalibrated_results/", scN,"/per_iter_results/", .) %>%
      map(readRDS) %>% 
      bind_rows() %>% 
      merge(pairs, by = "pair_id")
  }
  
  tbl <- 
    as.data.frame(table(df$pair_id)) %>% 
    mutate(pair_id = as.numeric(as.character(Var1))) %>%
    dplyr::select(pair_id, Freq) %>%
    arrange(pair_id, decreasing = F)

  
  summary_stats <- 
    df %>% 
    group_by(pair_id) %>% 
    mutate(
      int = as.numeric(int), 
      slp = as.numeric(slp)) %>% 
    summarise(
      auc_med  = median(auc, na.rm = T), 
      auc_sd   = sd  (auc, na.rm = T), 
      bri_med  = median(bri, na.rm = T), 
      bri_sd   = sd  (bri, na.rm = T), 
      int_med  = median(int, na.rm = T), 
      int_sd   = sd  (int, na.rm = T),
      slp_med  = median(slp, na.rm = T), 
      slp_sd   = sd  (slp, na.rm = T)
    ) 
  
  out <- list("overall_dataframe" = df, "rows_per_pair" = tbl,  "summary"= summary_stats)
  
  if (recalibrated == F){
    saveRDS(out, paste0("./../visualize_results/results/results_", scN ,".RData"))
  }
  if (recalibrated == T){
    saveRDS(out, paste0( "./../visualize_results/results/results_", scN ,"_recalibrated.RData"))
  }
} 

# ------------------ re-calibration --------------------------------------------

recalibrate_everything <- function(i, scN){
  
  # input: 
  # i - an integer between 1 and 2000 
  # this integer represents the simulation iteration
  #
  # scN - a character specifying the simulation scenario 
  # here N can take on values between 1 and 18 (e.g., sc3)
  # 
  # for a given iteration i in scenario N
  # the function imports the raw predicted risks from the appropriate folder, 
  # implements logistic re-calibration for each prediction model 
  # and 
  # 1. stores the re-calibrated predicted risks in a new file 
  # 2. calculates and stores all empirical performance metrics 
  #    based on the re-calibrated predicted risks in a new file
  
  
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
  
  # input: data frame 
  #
  # it is necessary for this data frame to have columns: 
  #
  # pair_id - a unique indicator of the prediction model (integer between 1-30)
  # pred  - the predicted risks for observations from one iteration
  # class - the true class values for observations from one iteration
  #
  # The function turns pair_id into a factor 
  # such that a re-calibration procedure can occur for each prediction model.
  #
  # It returns a dataframe identical to the dataframe that was input
  # except in the output dataframe, pred houses the recalibrated predicted risks
  # instead of the raw predicted risks
  
  
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
  
  # retrieve re-calibrated probabilities from one data set
  r_preds <- recalibrate_one_dataset(dataframe$pred, dataframe$class)
  
  # save warning from glm during recalibration
  r_warn <<- ifelse(is.null(r_preds$warning) == TRUE, "0", as.character(r_preds$warning))
  
  # update probabilities to re-calibrated probabilities
  dataframe <- 
    dataframe %>% 
    mutate(
      pred = r_preds$pred,
      recalibrated = 1
    ) %>% 
    mutate(class = as.numeric(class) -1)
  
  return(dataframe)
}

recalibrate_one_dataset <- function(probs, class){
  
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
  
  # this function takes a dataframe as input and re-orders and renames
  # the factors correction and algorithm 
  # such that the names are suitable for calibration plots
  
  
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
  
  # this function takes a data frame as input and ensures that 
  # all variables are in the correct form to generate calibration plots, 
  # any NAs are dropped, as specified in Supplementary Material section B
  
  dataframe %>% 
  drop_na() %>% 
  mutate(
      pred  = as.numeric(pred), 
      class = as.numeric(class)) %>%
  reorder_and_rename()
}

loess_regression<- function(dataframe){
  
  # input: dataframe 
  # 
  # it is necessary for this data frame to have columns:
  # pred  - predicted risks for one data set 
  # class - true observed class for for one data set
  #
  # this function conducts loess regression with class ~ pred for one data set
  # and saves 200 coordinates of the resulting loess curve in a new file
  
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
  
  # input: 
  #
  # i - integer between 1 and 2000 indicating the iteration 
  # scN - a indicator of the simulation scenario (e.g., "sc3")
  # recalibrated - a boolean indicating if processing raw or recalibrated results
  #
  # for iteration i in scenario N  this function 
  # imports the relevant raw or recalibrated predictions
  # computes loess regression / saves plot coordinates for each prediction model 
  # stores resulting coordinates for the iteration in a new file. 
  
  
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
  
  # input: a data frame 
  # it is necessary for the data frame to have columns 
  #
  # iter - indicating the iteration 
  # y - indicating the y coordinate for a given flexible calibration curve
  # 
  # this function prevents errors in the calibration_plot() function. 
  # if there are too few valid coordinates for a given flexible calibration curve
  # it cannot be represented in a ggplot 
  #
  # this function counts the number of coordinates for each flexible calibration curve
  # and removes coordinates for curves with fewer than 4 coordinates l
  
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

plot_from_coords <- function(df, hex = "#101011", restricted_range = F){
  
  # input: 
  # 
  # df - a dataframe
  # it is necessary that this data frame be generated by merging all files
  # in the plot_coords folder for a given scenario 
  #
  # hex - a hex code indicating the color of the calibration curves 
  #
  # restricted_range - a boolean specifying if the calibration plot 
  # should have a restricted range on the x-axis. 
  #
  # restricted_range = F: the x-axis ranges from 0 to 1 
  # restricted_range = T: the x-axis ranges from 0 to 0.25
  #
  # 
  # This function generates the plot object of a calibration plot for a given scenario 
  
  
  if(restricted_range == F){
    p <-
      ggplot() +
      geom_abline(slope = 1, intercept = 0, linewidth = 1, color = "darkgray") +
      # geom_smooth(data = df, aes(x = x, y = y, group = iter, color = chance_imbalance), 
      geom_smooth(data = df, aes(x = x, y = y, group = iter), 
                  method = "loess",
                  formula = y ~ x, 
                  se = F,
                  linewidth = 0.01, 
                  color = hex,
                  alpha = 0.01) + 
      theme_minimal() +
      xlim(0,1) + 
      ylim(0,1) + 
      facet_grid(algorithm ~ correction) +
      theme(legend.position = "none") + 
      xlab("Predicted Risk") + 
      ylab("Observed Proportion") 
  }
  
  if(restricted_range == T){
    p <-
      ggplot() +
      geom_abline(slope = 1, intercept = 0, linewidth = 1, color = "darkgray") +
      # geom_smooth(data = df, aes(x = x, y = y, group = iter, color = chance_imbalance), 
      geom_smooth(data = df, aes(x = x, y = y, group = iter), 
                  method = "loess",
                  formula = y ~ x, 
                  se = F,
                  linewidth = 0.01, 
                  color = hex,
                  alpha = 0.01) + 
      theme_minimal() +
      facet_grid(algorithm ~ correction) +
      theme(legend.position = "none") + 
      xlab("Predicted Risk") + 
      ylab("Observed Proportion") + 
      coord_cartesian(xlim = c(0,0.25), ylim = c(0,0.25))
  }
  
  return(p) 
}

calibration_plot <- function(i, recalibrated = F){
  
  # input 
  # 
  # i - a single integer representing the simulation scenario 
  # recalibrated - boolean representing if processing raw or re-calibrated results
  #
  # this function generates and saves a calibration plot, 
  # with flexible calibration curves for all prediction models 
  # and all simulation iterations for a given simulation scenario. 
  
  
  scN <- paste0("sc", i)
  
  if (recalibrated == F){
    df <- 
      list.files(path = paste0("./processed_results/plot_coords/raw_results/", scN, "/"), pattern = ".rds") %>%
      paste0("./processed_results/plot_coords/raw_results/", scN, "/", .) %>%
      map(readRDS) %>% 
      bind_rows() %>% 
      as.data.frame() %>% 
      drop_na() %>%
      pre_plot() 
  
    p <- plot_from_coords(df)
  
    ggsave(paste0("./../visualize_results/results/calibration_plot_", scN, ".png"), 
           plot = p,  width = 14, height = 14)
  
  }
  
  if(recalibrated == T){
    df <- 
      list.files(path = paste0("./processed_results/plot_coords/recalibrated_results/", scN, "/"), pattern = ".rds") %>%
      paste0("./processed_results/plot_coords/recalibrated_results/", scN, "/", .) %>%
      map(readRDS) %>% 
      bind_rows() %>% 
      as.data.frame() %>%
      drop_na() %>%
      pre_plot() 
    
    p <- plot_from_coords(df)
    
    ggsave(paste0("./../visualize_results/results/calibration_plot_", scN, "_recalibrated.png"), 
           plot = p,  width = 14, height = 14)
  }
}


# ----------------performance metric plots -------------------------------------

tidy_for_pm_plts <- function(dataframe){
  
  # this function takes a dataframe as input and tidies it such that the 
  # the variables are in the appropriate form to generate violin plots
  
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
                            'sen', 'smo', 'ros', 'rus', 'control'))) %>%
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
                          "sen"    = "SENN")) %>%
    mutate(int = as.numeric(int), 
           slp = as.numeric(slp))
}

pm_plot <- function(dataframe, method, xmin, xmax){
  
  # input: 
  # 
  # dataframe - a data frame 
  # method - the empirical performance metric for visualization 
  # xmin - the minimum value for the range on the x-axis 
  # xmax - the maximum value for the range on the x-axis 
  #
  # this function generates violin plots visualizing the results across 
  # 2000 simulation iterations for a given empirical performance metric 
  # for all 30 prediction models in a given scenario
  
  if(method == "auc"){
    title <- "Concordance Statistic:"
    ref   <- 0.85
  }
  if(method == "slp"){
    title <- "Calibration Slope:"
    ref   <- 1
  }
  if(method == "int"){
    title <- "Calibration Intercept:"
    ref <- 0
  }
  if(method == "bri"){
    title <- "Brier Score:"
    ref <- 0
  }
  
  dataframe %>% 
    tidy_for_pm_plts() %>%
    ggplot() + 
    geom_hline(yintercept = ref, color = "darkgray", linewidth = 1) + 
    geom_violin(aes_string(x = "correction", y = method),
                fill  = "#b3cd41", 
                color = "grey",
                draw_quantiles = 0.5,
                linewidth = 0.25,
                alpha = 0.4) + 
    theme_minimal() +
    xlab(" ") + 
    ylab(" ") +
    ylim(xmin, xmax) +
    coord_flip() + 
    facet_wrap(~algorithm, nrow = 1) + 
    ggtitle(title) + 
    theme(text=element_text(family="Times"))
}

# ------------------- scenario plots -------------------------------------------

reorder_and_rename2 <- function(dataframe){
  
  # this function takes a data frame as input 
  # and tidies it such that the variables algorithm and correction 
  # have the appropriate names for the plots displayed in Appendix A
  
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
                         "logistic_regression"    = "LR",
                         "support_vector_machine" = "SVM",
                         "random_forest"          = "RF",
                         "xgboost"                = "XG" ,
                         "rusboost"               = "RB",
                         "easy_ensemble"          = "EE"),
      correction = recode(correction,
                          "control" = "Control",
                          "rus"    = "RUS",
                          "ros"    = "ROS",
                          "smo"    = "SMOTE",
                          "sen"    = "SENN"))
}

all_summaries <- function(recalibrated = F){
  
  # this function imports all the results (summaries) for all simulation scenarios 
  # and merges them into one data frame 
  #
  # the result is a data frame which stores the median value of each empirical 
  # performance metric, for all prediction models, in all simulation scenarios
  
  
  tutto <- data.frame()
  
  if(recalibrated == F){
    for(i in 1:18){
    
      df <- paste0("./../visualize_results/results/results_sc", i, ".RData") %>% readRDS() 
      df <- df$summary %>% mutate(scenario = i) %>% relocate(scenario) 
      df <- df %>% mutate(pair_id = pair_id %>% as.numeric()) %>% arrange(pair_id)
    
      tutto <- bind_rows(tutto, df) 
    }
  }
  
  if(recalibrated == T){
    for(i in 1:18){
      
      df <- paste0("./../visualize_results/results/results_sc", i, "_recalibrated.RData") %>% readRDS() 
      df <- df$summary %>% mutate(scenario = i) %>% relocate(scenario) 
      df <- df %>% mutate(pair_id = pair_id %>% as.numeric()) %>% arrange(pair_id)
      
      tutto <- bind_rows(tutto, df) 
    }
  }
  
  return(tutto)
}

scenario_graph <- function(dataframe, method){
  
  # this function produces the figures displayed in Appendix A 
  # 
  # input: 
  # dataframe - must be the ouput from the function all_summaries()
  # method - must specify the performance metric desired for visualization
  #
  # output: 
  # figures displayed in Appendix A
  
  if(method == "auc_med"){
    title <- "Concordance Statistic"
    ref   <- 0.85
  }
  if(method == "slp_med"){
    title <- "Calibration Slope"
    ref   <- 1
  }
  if(method == "int_med"){
    title <- "Calibration Intercept"
    ref <- 0
  }
  if(method == "bri_med"){
    title <- "Brier Score"
    ref <- 0
  }
  
  sim_sets  <- sim_sets %>% mutate(scenario = sc)
  
  df <- 
    dataframe %>% 
    merge(pairs, by = "pair_id") %>%
    merge(sim_sets, by = "scenario") %>%
    mutate(ef = as.factor(ef)) %>%
    mutate(
      ef = recode(ef, 
                  "0.5"  = "Event Fraction = 0.5",
                  "0.2"  = "Event Fraction = 0.2",
                  "0.02" = "Event Fraction = 0.02"
      ))%>% 
    mutate(npred = as.factor(npred)) %>% 
    mutate(
      npred = recode(npred, 
                     "8"  = "8 Predictors", 
                     "16" = "16 Predictors"
      ))
  
  ggplot(data = df, aes_string(x = "pair_id", y = method , color = "n")) + 
    geom_point(size = 0.6) + 
    geom_line(linewidth = 0.2) + 
    geom_hline(yintercept = ref, color = "darkgray", linewidth = 1) + 
    coord_flip() + 
    facet_nested(npred ~ ef) + 
    scale_color_manual("Sample Size", values = c("#7CB4B8", "#7F95D1", "#344055"))+
    theme_minimal() + 
    theme(panel.background = element_rect(
      color="#FEFEFA", fill="#FFFFFF", size=1.5, linetype="solid"
    ),
    strip.background = element_rect(
      color="#FEFEFA", fill="#F5F5DC", size=1.5, linetype="solid"
    ), 
    panel.spacing.y = unit(3, "lines"),
    panel.spacing.x = unit(3, "lines"),
    legend.key = element_rect( color="#FEFEFA", fill = "#FFFFFF"),
    text=element_text(family="Times", size = 20),
    axis.text.y = element_text(hjust=0.95),
    strip.text = element_text(size = rel(1.25))) + 
    xlab("Prediction Model") + 
    ylab(title) + 
    scale_x_discrete(limits = c(
      "Control-LR", "Control-SVM", "Control-RF", "Control-XG", "Control-RB", "Control-EE", 
      "RUS-LR", "RUS-SVM", "RUS-RF", "RUS-XG", "RUS-RB", "RUS-EE", 
      "ROS-LR", "ROS-SVM", "ROS-RF", "ROS-XG", "ROS-RB", "ROS-EE", 
      "SMOTE-LR", "SMOTE-SVM", "SMOTE-RF", "SMOTE-XG", "SMOTE-RB", "SMOTE-EE", 
      "SENN-LR", "SENN-SVM", "SENN-RF", "SENN-XG", "SENN-RB", "SENN-EE"))
}

# please note: the figures displayed in Appendix A have a slightly different 
# y-axis label appearance.  The graph is identical, but y-axis labels were edited
# in power point to improve readability. 


#-------------------- error handling -------------------------------------------

# create a data frame which functions as a key for the simulation scenario 
#
# rows: 18, one for each simulation scenario 
# columns: 
#
# npred: number of predictors 
# ef   : event fraction 
# n    : the setting for sample size 
# sc   : an indicator of the simulation scenario

npred <- c(8,16)
ef    <- c(0.5,0.2,0.02)
n     <- c("0.5N", "N", "2N")
sets  <- expand.grid("ef" = ef, "n" = n, "npred" = npred)

sim_sets <- sets %>% mutate(sc = c(1:18))
sim_sets <- sim_sets[1:18,]

#-------------------------------------------------------------------------------
# error handling continued

non_zeros <- function(x){
  
  # a function which counts all non-zero entries in a column
  # and returns the total sum in a dataframe with a single value
  
  count <- sum(x!="0")
  return(as.data.frame(count))
  
} 

count_correction_errors <- function(dataframe){
  dataframe %>%
    select(-seed, -sc_ef, -auc, -s_bri, -bri, -int, -slp, -min, -max) %>% 
    mutate(correction = as.factor(correction)) %>% 
    group_by(correction) %>% 
    group_modify(~non_zeros(.x$c_err)) %>% 
    mutate(count = count / 6) 
}

count_alg_errors <- function(dataframe){
  dataframe %>% 
    select(-seed, -sc_ef, -auc, -s_bri, -bri, -int, -slp, -min, -max) %>% 
    mutate(pair_id = as.factor(pair_id)) %>% 
    group_by(pair_id) %>% 
    group_modify(~non_zeros(.x$a_err))
}

count_rb_invalid <- function(dataframe){
  dataframe %>% 
    select(-seed, -sc_ef, -auc, -s_bri, -bri, -int, -slp, -min, -max) %>% 
    mutate(pair_id = as.factor(pair_id)) %>% 
    mutate(invalid = invalid_0 + invalid_1) %>% 
    group_by(pair_id) %>% 
    group_modify(~non_zeros(.x$invalid))
}

count_lg_separation <- function(dataframe){
  dataframe %>% 
    select(-seed, -sc_ef, -auc, -s_bri, -bri, -int, -slp, -min, -max) %>% 
    mutate(pair_id = as.factor(pair_id)) %>% 
    group_by(pair_id) %>% 
    group_modify(~non_zeros(.x$a_warn))
}

save_errors_corr <- function(i){
  
  scN <- paste0("sc", i)
  df  <- readRDS(paste0("./../visualize_results/results/results_", scN, ".RData"))$overall_dataframe
  
  err <- count_correction_errors(df) %>% pull(count)
  return(err)
  
}

save_errors_alg <- function(i, name){
  
  scN <- paste0("sc", i)
  df  <- readRDS(paste0("./../visualize_results/results/results_", scN, ".RData"))$overall_dataframe
  
  err <- 
    count_alg_errors(df) %>% 
    merge(pairs, by = "pair_id") %>% 
    filter(algorithm == name) %>%
    mutate(pair_id = pair_id %>% as.character() %>% as.numeric()) %>% 
    arrange(pair_id, decreasing = F) %>% 
    pull(count)
  
  return(err)
}

save_errors_rb <- function(i){
  
  scN <- paste0("sc", i)
  df  <- readRDS(paste0("./../visualize_results/results/results_", scN, ".RData"))$overall_dataframe
  
  err <- count_rb_invalid(df) %>% 
    merge(pairs, by = "pair_id") %>% 
    mutate(pair_id = pair_id %>% as.character() %>% as.numeric()) %>% 
    arrange(pair_id) %>% 
    filter(algorithm == "rusboost") %>% 
    pull(count)
  
  return(err)
  
}

save_lg_separation <- function(i){
  
  scN <- paste0("sc", i)
  df  <- readRDS(paste0("./../visualize_results/results/results_", scN, ".RData"))$overall_dataframe
  
  err <- 
    count_lg_separation(df) %>% 
    merge(pairs, by = "pair_id") %>% 
    mutate(pair_id = pair_id %>% as.character() %>% as.numeric()) %>% 
    arrange(pair_id) %>% 
    filter(algorithm == "logistic_regression") %>% 
    pull(count)
  
  return(err)
  
}

make_table_corrections <- function(){
  
  corr_err <- data.frame()

  for (i in 1:18){
    row      <- save_errors_corr(i) %>% t() %>% as.data.frame() 
    corr_err <- bind_rows(corr_err, row)
  }

  colnames(corr_err) <- c("control", "rus", "ros", "smo", "sen")

  table <- 
    corr_err %>% 
    mutate(sc = c(1:18)) %>% 
    merge(sim_sets, by = "sc") %>% 
    relocate(sc, npred, n, ef)

  colnames(table) <- c("Scenario", "No. Predictors", "Sample Size", "Event Fraction", 
                     "Control", "RUS", "ROS", "SMOTE", "SENN")

  saveRDS(table, "./../manuscript/error_table_corrections.RDS")
}

make_table_algorithm <- function(algorithm){
  
  alg_err <- data.frame()
  
  for (i in 1:18){
    row      <- save_errors_alg(i, algorithm) %>% t() %>% as.data.frame() 
    alg_err  <- bind_rows(alg_err, row)
  }
  
  colnames(alg_err) <- c("control", "rus", "ros", "smo", "sen")
  
  table <- 
    alg_err %>% 
    mutate(sc = c(1:18)) %>% 
    merge(sim_sets, by = "sc") %>% 
    relocate(sc, npred, n, ef)
  
  colnames(table) <- c("Scenario", "No. Predictors", "Sample Size", "Event Fraction", 
                       "Control", "RUS", "ROS", "SMOTE", "SENN")
  
  print(table)
  saveRDS(table, paste0("./../manuscript/error_table_", algorithm, ".RDS"))
}

make_table_rb_invalid <- function(){
  
  invalid_probs <- data.frame()
  
  for (i in 1:18){
    row            <- save_errors_rb(i) %>% t() %>% as.data.frame() 
    invalid_probs  <- bind_rows(invalid_probs, row)
  }
  
  colnames(invalid_probs) <- c("control", "rus", "ros", "smo", "sen")
  
  table <- 
    invalid_probs %>% 
    mutate(sc = c(1:18)) %>% 
    merge(sim_sets, by = "sc") %>% 
    relocate(sc, npred, n, ef)
  
  colnames(table) <- c("Scenario", "No. Predictors", "Sample Size", "Event Fraction", 
                       "Control", "RUS", "ROS", "SMOTE", "SENN")
  
  print(table)
  saveRDS(table, paste0("./../manuscript/rb_invalid_probs.RDS"))
}

make_table_lg_separation <- function(){
  
  invalid_probs <- data.frame()
  
  for (i in 1:18){
    row            <- save_lg_separation(i) %>% t() %>% as.data.frame() 
    invalid_probs  <- bind_rows(invalid_probs, row)
  }
  
  colnames(invalid_probs) <- c("control", "rus", "ros", "smo", "sen")
  
  table <- 
    invalid_probs %>% 
    mutate(sc = c(1:18)) %>% 
    merge(sim_sets, by = "sc") %>% 
    relocate(sc, npred, n, ef)
  
  colnames(table) <- c("Scenario", "No. Predictors", "Sample Size", "Event Fraction", 
                       "Control", "RUS", "ROS", "SMOTE", "SENN")
  
  print(table)
  saveRDS(table, paste0("./../manuscript/lg_separation.RDS"))
}