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

# ------------------ re-calibration --------------------------------------------

recalibrate_everything <- function(i, scN){
  
  #  import probabilities from simulation
  df <- 
    i %>% 
    paste0("./sim_results/raw_results/predictions/", scN, "_", . , ".rds") %>% 
    map(readRDS) 
  
  # save re-calibrated probabilities
  r_df <- 
    df %>%
    recalibrate_one_iteration()
  
  saveRDS(r_df, paste0("./sim_results/recalibrated_results/predictions/", scN, "_", i, ".rds"))
  
  # save re-calibrated per-iteration statistics
  r_stats <- 
    r_df %>% 
    get_stats_recalibrated() 
  
  saveRDS(r_stats, paste0("./sim_results/recalibrated_results/per_iter_results/", scN, "_", i, ".rds"))
  
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
      paste0("./sim_results/recalibrated_results/predictions/", scN, "_", . , ".rds") %>% 
      map(readRDS) %>% 
      as.data.frame() %>%
      ready2plot() %>% 
      group_by(pair_id) %>% 
      group_modify(~loess_regression(.)) %>% 
      merge(pairs, by = "pair_id") %>% 
      reorder_and_rename() %>% 
      mutate(iter = as.factor(iter))
    
      saveRDS(df, paste0("./plot_coords/recalibrated_results/", scN, "_",  i, "_plot_coords.rds"))
  }
  
  if(recalibrated == F){
    df <- 
      i %>% 
      paste0("./sim_results/raw_results/predictions/", scN, "_", . , ".rds") %>% 
      map(readRDS) %>% 
      as.data.frame() %>%
      ready2plot() %>% 
      group_by(pair_id) %>% 
      group_modify(~loess_regression(.)) %>% 
      merge(pairs, by = "pair_id") %>% 
      reorder_and_rename() %>% 
      mutate(iter = as.factor(iter))
    
    saveRDS(df, paste0("./plot_coords/raw_results/", scN, "_",  i, "_plot_coords.rds"))
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


