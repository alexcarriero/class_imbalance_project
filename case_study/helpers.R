
# imbalance corrections --------------------------------------------------------

do_rus <- function(dataframe){
  
  message <- 0
  count   <- 0
  warn    <- 0
  formula <- y~.
  data    <- dataframe
  
  df.rus <- 
    tryCatch.W.E(
      expr = {
        df.rus <- ovun.sample(formula, data, method = "under")$data
      }
    )
  
  if(is.data.frame(df.rus$value) == TRUE){
    c_data  <- df.rus$value
  } else{
    message <- as.character(df.rus$value)
    count   <- 1
    c_data  <- dataframe
  }
  
  warn <- as.character(df.rus$warning)
  out  <- list("c_dataframe" = c_data, "correction_err" = message, "count_err" = count, "correction_warning" = warn)
  return(out)
}

do_ros <- function(dataframe){
  
  message <- 0
  count   <- 0
  warn    <- 0
  formula <- outcome~.
  data    <- dataframe
  
  df.ros <- 
    tryCatch.W.E(
      expr = {
        df.ros <-  ovun.sample(formula, data, method = "over")$data
      }
    )
  if(is.data.frame(df.ros$value) == TRUE){
    c_data  <- df.ros$value
  } else{
    message <- as.character(df.ros$value)
    count   <- 1
    c_data  <- dataframe
  }
  
  warn <- as.character(df.ros$warning)
  out  <- list("c_dataframe" = c_data, "correction_err" = message, "count_err" = count, "correction_warning" = warn)
  return(out)
}

do_smo <- function(df, percOver, k = 5){
  
  message <- 0
  count   <- 0
  warn    <- 0
  
  df.smote <- 
    tryCatch.W.E(
      expr = {
        df.smote <- SMOTE(x = df[-1], y = df[1], percOver, k)
      }
    )
  if(is.data.frame(df.smote$value) == TRUE){
    c_data  <- df.smote$value %>% relocate("outcome")
  } else{
    message <- as.character(df.smote$value)
    count   <- 1
    c_data  <- df
  }
  
  warn <- as.character(df.smote$warning)
  out  <- list("c_dataframe" = c_data, "correction_err" = message, "count_err" = count, "correction_warning" = warn)
  return(out)
}

do_sen <- function(df, percOver, k1 = 5, k2 = 3){
  
  message <- 0
  count   <- 0
  warn    <- 0
  
  df.smote.enn <- 
    tryCatch.W.E(
      expr = {
        df.smote.enn <- SmoteENN(x = df[-1], y = df[1], percOver, k1, k2)
      }
    )
  if(is.data.frame(df.smote.enn$value) == TRUE){
    c_data  <- df.smote.enn$value %>% relocate("outcome")
  } else{
    message <- as.character(df.smote.enn$value)
    count   <- 1
    c_data  <- df
  }
  
  warn <- as.character(df.smote.enn$warning)
  out  <- list("c_dataframe" = c_data, "correction_err" = message, "count_err" = count, "correction_warning" = warn)
  return(out)
}



# model development ------------------------------------------------------------

deviance <- function(data, lev = NULL, model = NULL) {
  obs  <- as.numeric(data$obs) - 1
  pred <- data$one
  
  pred[pred == 0] <- 1e-16
  pred[pred == 1] <- 1-1e-16
  
  dev <- -2*sum(obs*log(pred) + (1-obs)*log(1-pred))
  
  c(Deviance = dev)
}


do_lrg <- function(df_train, df_test){
  
  df_train$y  <- factor(df_train$y, levels=c(0,1), ordered = TRUE)
  
  # model 
  mod  <- glm( y~., family = "binomial", data = df_train)
  pred <- predict(mod,  newdata = df_test, type = "response")
  
  
  # output 
  output   <- cbind(pred, "class" = as.numeric(df_test$y)-1) %>% as.data.frame() %>%  `rownames<-`( NULL )
  return(output)
}


do_svm <- function(df_train, df_test){
  
  df_train$y     = factor(df_train$y, levels=c(0,1), ordered = TRUE)
  
  # model
  levels(df_train$y) <- c("zero", "one")
  
  train_ctrl <- 
    trainControl(
      method = "cv", 
      number = 5, 
      summaryFunction = deviance, 
      classProbs = T,
      allowParallel = F
    )
  
  mod <- 
    train(
      y ~., 
      data = df_train, 
      method = "svmRadial", 
      metric = "Deviance",
      maximize = FALSE, 
      trControl = train_ctrl,
      preProcess = c("center", "scale"), 
      tuneLength = 3
    ) 
  
  pred <- predict(mod, newdata = df_test, type = "prob")["one"]
  colnames(pred) <- "pred"
  pred <- pred
  
  # output 
  output   <- cbind(pred, "class" = as.numeric(df_test$y)-1)
  return(output)
}

do_rnf <- function(df_train, df_test){
  
  df_train$y     = factor(df_train$y, levels=c(0,1), ordered = TRUE)
  
  # model
  npred = dim(df_train)[2] -1 
  levels(df_train$y) <- c("zero", "one")
  
  train_ctrl <- 
    trainControl(
      method = "cv", 
      number = 5, 
      summaryFunction = deviance, 
      classProbs = T,
      allowParallel = F
    )
  
  mod <- 
    train(
      y ~., 
      data = df_train, 
      method = "ranger", 
      metric = "Deviance",
      maximize = FALSE, 
      trControl = train_ctrl,  
      tuneGrid = 
        expand.grid(
          mtry = 1:npred,
          min.node.size = 1:10,
          splitrule = "gini"
        )
    )    
  
  pred <- predict(mod, newdata = df_test, type = "prob")["one"]
  colnames(pred) <- "pred"
  pred <- pred
  
  
  # output 
  output   <- cbind(pred, "class" = as.numeric(df_test$y)-1)
  return(output)
}

do_xgb <- function(df_train, df_test){
  
  df_train$y   = factor(df_train$y, levels=c(0,1), ordered = TRUE)
  
  # model

      levels(df_train$y) <- c("zero", "one")
      
      train_ctrl <- 
        trainControl(
          method = "cv", 
          number = 5, 
          summaryFunction = deviance, 
          classProbs = T,
          allowParallel = F
        )
      
      mod <- 
        train(
          y ~., 
          data = df_train, 
          method = "xgbTree", 
          metric = "Deviance",
          maximize = FALSE, 
          trControl = train_ctrl, 
          verbosity = 0, 
          verbose = FALSE, 
          tuneLength = 3
        )
      
      pred <- predict(mod, newdata = df_test, type = "prob")["one"]
      colnames(pred) <- "pred"
      pred <- pred
  
  # output 
  output   <- cbind(pred, "class" = as.numeric(df_test$y)-1)
  return(output)
}


do_rub <- function(df_train, df_test){
  
  # model 
      mod  <- rus(y ~., size = 10, alg = "svm", data = df_train, svm.ker = "radial")
      pred <- predict(mod, newdata = df_test)
  
  # output 
  output <- cbind(pred, "class" = as.numeric(df_test$y)-1)
  return(output)
}


do_eas <- function(df_train, df_test){

      train_x <- df_train %>% dplyr::select(-y)
      train_y <- df_train %>% dplyr::pull(y)
      test_x  <- df_test  %>% dplyr::select(-y)
      
      mod  <- EasyEnsemble(train_x, train_y)
      pred <- predict(mod, test_x, type = "probability")$X1
  
  # output 
  output   <- cbind(pred, "class" = as.numeric(df_test$y)-1)
  return(output)
  
}



# performance measures ---------------------------------------------------------


auroc <- function (probs, outcome){
  if(any(is.na(probs)) == TRUE){
    return(NA)
  }
  
  if(is.numeric(outcome) == FALSE) {outcome <- as.numeric(outcome) - 1}
  
  auc  <- as.numeric(pROC::roc(outcome~probs, quiet = TRUE)$auc)
  
  return(auc)
}

brier_score <- function(probs, outcome){
  if(any(is.na(probs)) == TRUE){
    return(NA)
  }
  
  if(is.numeric(outcome) == FALSE) {outcome <- as.numeric(outcome) - 1}
  
  obs = outcome
  brier_score  = mean((probs - obs)^2)
  return(brier_score)
}

calibration_intercept <- function(probs, outcome){
  if(any(is.na(probs)) == TRUE){
    return(NA)
  }
  if(sum(probs == 1) != 0){
    probs[probs == 1] <- 1 - 1e-10
  }
  if(sum(probs == 0) != 0){
    probs[probs == 0] <- 1e-10
  }
  mod <- glm(outcome ~ 1, offset = log(probs/(1-probs)), family = "binomial")
  int <- coef(mod)[1]
  return(int)
}

calibration_slope <- function(probs, outcome){
  if(any(is.na(probs)) == TRUE){
    return(NA)
  }
  if(sum(probs == 1) != 0){
    probs[probs == 1] <- 1 - 1e-10
  }
  if(sum(probs == 0) != 0){
    probs[probs == 0] <- 1e-10
  }
  mod <- glm(outcome ~ log(probs/(1-probs)), family = "binomial")
  slp <- coef(mod)[2]
  return(slp)
}

performance_statistics <- function(pred, class){
  tibble(
    min = min(pred),
    max = max(pred),
    var = var(pred),
    auc = auroc(pred, class),
    bri = brier_score(pred, class),
    int = calibration_intercept(pred, class),
    slp = calibration_slope(pred, class)
  )
}

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
                          "sen"    = "SENN"))
}

tidy_for_pm_plts <- function(dataframe){
  
  # this function takes a dataframe as input and tidies it such that the 
  # the variables are in the appropriate form to generate violin plots
  
  dataframe %>% 
    mutate(
      algorithms  = factor(algorithms,
                          levels = c(
                            'logistic_regression',
                            'support_vector_machine',
                            'random_forest',
                            'xgboost',
                            'rusboost',
                            'easy_ensemble')),
      corrections = factor(corrections,
                          levels = c(
                            'sen', 'smo', 'ros', 'rus', 'control'))) %>%
    mutate(
      algorithms = recode(algorithms,
                         "logistic_regression"    = "Logistic Regression",
                         "support_vector_machine" = "Support Vector Machine",
                         "random_forest"          = "Random Forest",
                         "xgboost"                = "XGBoost" ,
                         "rusboost"               = "RUSBoost",
                         "easy_ensemble"          = "EasyEnsemble"),
      corrections = recode(corrections,
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
    ref   <- 1
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
    geom_hline(yintercept = ref, color = "grey", linewidth = 1) + 
    geom_point(aes_string(x = "corrections", y = method),
                color = "#48494B") + 
    theme_minimal() +
    xlab(" ") + 
    ylab(" ") +
    ylim(xmin, xmax) +
    coord_flip() + 
    facet_wrap(~algorithms, nrow = 1) + 
    ggtitle(title) + 
    theme(text=element_text(family="Times"))
}
