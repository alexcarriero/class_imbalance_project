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

get_brier_score <- function(dataframe){
  tibble(
    brier_score = brier_score(dataframe$pred, dataframe$class)
  )
}

save_brier_scores <- function(i){
  
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
    
    saveRDS(df, paste0("./../masters_thesis/visualize_results/results/results_", scN ,"_new.RData"))
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
    saveRDS(out, paste0("./../masters_thesis/visualize_results/results/results_", scN ,".RData"))
  }
  if (recalibrated == T){
    saveRDS(out, paste0( "./../masters_thesis/visualize_results/results/results_", scN ,"_recalibrated.RData"))
  }
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

plot_from_coords <- function(df, hex = "#101011", restricted_range = F){
  
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
  
    ggsave(paste0("./../masters_thesis/visualize_results/results/calibration_plot_", scN, ".png"), 
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
    
    ggsave(paste0("./../masters_thesis/visualize_results/results/calibration_plot_", scN, "_recalibrated.png"), 
           plot = p,  width = 14, height = 14)
  }
}


# ----------------performance metric plots -------------------------------------

tidy_for_pm_plts <- function(dataframe){
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
  
  tutto <- data.frame()
  
  if(recalibrated == F){
    for(i in 1:18){
    
      df <- paste0("./../masters_thesis/visualize_results/results/results_sc", i, ".RData") %>% readRDS() 
      df <- df$summary %>% mutate(scenario = i) %>% relocate(scenario) 
      df <- df %>% mutate(pair_id = pair_id %>% as.numeric()) %>% arrange(pair_id)
    
      tutto <- bind_rows(tutto, df) 
    }
  }
  
  if(recalibrated == T){
    for(i in 1:18){
      
      df <- paste0("./../masters_thesis/visualize_results/results/results_sc", i, "_recalibrated.RData") %>% readRDS() 
      df <- df$summary %>% mutate(scenario = i) %>% relocate(scenario) 
      df <- df %>% mutate(pair_id = pair_id %>% as.numeric()) %>% arrange(pair_id)
      
      tutto <- bind_rows(tutto, df) 
    }
  }
  
  return(tutto)
}

scenario_graph <- function(dataframe, method){
  
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

#-------------------- error handling -------------------------------------------

npred <- c(8,16)
ef    <- c(0.5,0.2,0.02)
n     <- c("0.5N", "N", "2N")
sets  <- expand.grid("ef" = ef, "n" = n, "npred" = npred)

sim_sets <- sets %>% mutate(sc = c(1:18))
sim_sets <- sim_sets[1:18,]


non_zeros <- function(x){
  
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
  df  <- readRDS(paste0("./../masters_thesis/visualize_results/results/results_", scN, ".RData"))$overall_dataframe
  
  err <- count_correction_errors(df) %>% pull(count)
  return(err)
  
}

save_errors_alg <- function(i, name){
  
  scN <- paste0("sc", i)
  df  <- readRDS(paste0("./../masters_thesis/visualize_results/results/results_", scN, ".RData"))$overall_dataframe
  
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
  df  <- readRDS(paste0("./../masters_thesis/visualize_results/results/results_", scN, ".RData"))$overall_dataframe
  
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
  df  <- readRDS(paste0("./../masters_thesis/visualize_results/results/results_", scN, ".RData"))$overall_dataframe
  
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

  saveRDS(table, "./../masters_thesis/manuscript/error_table_corrections.RDS")
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
  saveRDS(table, paste0("./../masters_thesis/manuscript/error_table_", algorithm, ".RDS"))
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
  saveRDS(table, paste0("./../masters_thesis/manuscript/rb_invalid_probs.RDS"))
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
  saveRDS(table, paste0("./../masters_thesis/manuscript/lg_separation.RDS"))
}