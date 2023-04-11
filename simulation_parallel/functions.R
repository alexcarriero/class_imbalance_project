
# data generation --------------------------------------------------------------

generate_data <- function(inn, test = FALSE){
  
  npred   = inn[[1]]
  ev      = inn[[2]]
  n.level = inn[[3]]
  n       = inn[[4]]
  mu0     = inn[[5]]
  mu1     = inn[[6]]
  sigma0  = inn[[7]]
  sigma1  = inn[[8]]
  
  # test set
  if (test == TRUE){n  = n*10}
  
  # positive class
  n1      <- rbinom(1, n, ev)
  class_1 <- mvrnorm(n1, mu1, sigma1)
  
  # negative class
  n0      <- n - n1
  class_0 <- mvrnorm(n0, mu0, sigma0)
  
  outcome <- c(rep(1, n1), rep(0, n0))
  
  # format data frame
  df <- cbind(outcome, rbind(class_1, class_0)) %>% 
    as.data.frame() %>% 
    mutate(outcome = as.factor(outcome))
  
  # names
  colnames(df) <- c("outcome", paste0("x", 1:npred))
  
  # standardize 
  # df <- df %>% mutate(across(where(is.numeric), scale))
  
  return(df)
}


# class imbalance corrections --------------------------------------------------

do_rus <- function(dataframe){
  k = 0
  
  formula <- outcome~.
  data    <- dataframe
  
  
  while(k < 10){
    
    k = k + 1
    
    df.rus <- 
      tryCatch(
        
        expr = {
          df.rus <-  ovun.sample(formula, data, method = "under")$data
          return(df.rus)
        },
        
        error = function(e){ 
          df.rus <- "problematic"
          return(df.rus)
        }
      )
    
    if (df.rus != "problematic")
      break
    
  }
  return(df.rus)
}


do_ros <- function(dataframe){
  k = 0
  
  formula <- outcome~.
  data    <- dataframe
  
  
  while(k < 10){
    
    k = k + 1
    
    df.ros <- 
      tryCatch(
        
        expr = {
          df.ros <-  ovun.sample(formula, data, method = "over")$data
          return(df.ros)
        },
        
        error = function(e){ 
          df.ros <- "problematic"
          return(df.ros)
        }
      )
    
    if (df.ros != "problematic")
      break
    
  }
  return(df.ros)
}


do_smo <- function(df, percOver, k = 5){
  df.smote <- SMOTE(x = df[-1], y = df[1], percOver, k)
  df.smote <-
    df.smote %>%
    relocate("outcome")
  return(df.smote)
}


do_sen <- function(df, percOver, k1 = 5, k2 = 3){
  df.smote.enn <- SmoteENN(x = df[-1], y = df[1], percOver, k1, k2)
  df.smote.enn <- 
    df.smote.enn %>% 
    relocate("outcome")
  return(df.smote.enn)
}


# model implementation ---------------------------------------------------------

do_lrg <- function(df, test, iter, s){
  
  error       <- 0
  warning     <- 0
  problematic <- 0
  
  # ensure correct factor levels 
  df$outcome  <- factor(df$outcome, levels=c(0,1), ordered = TRUE)
  
  # model 
  a <- tryCatch.W.E(
    expr    = { 
      set.seed(s)
      mod  <- glm( outcome~., family = "binomial", data = df)
      pred <- predict(mod,  newdata = test, type = "response")
    }
  )
  
  # sort predictions and errors
  if (is.vector(a$value) == FALSE){
    pred  <- NA
    error <- a$value
    problematic <- 1
  } else{
    pred <- a$value
  }
  
  # sort warnings
  if(is.null(a$warning) == FALSE){
    warning <- a$warning
    problematic <- 1
  }
  
  # output 
  pred   <- cbind(pred, "iter" = iter, "class" = as.numeric(test$outcome) - 1, "problematic" = problematic)
  output <- list("pred" = pred, "warning" = warning, "error" = error, "problem" = problematic)
  return(output)
}


deviance <- function(data, lev = NULL, model = NULL) {
  obs  <- as.numeric(data$obs) - 1
  pred <- data$one
  
  pred[pred == 0] <- 1e-16
  pred[pred == 1] <- 1-1e-16
  
  dev <- -2*sum(obs*log(pred) + (1-obs)*log(1-pred))
  
  c(Deviance = dev)
}


do_svm <- function(df, test, iter, s){
  
  error       <- 0
  warning     <- 0
  problematic <- 0
  
  # ensure correct factor levels, and character labels 
  df$outcome     = factor(df$outcome, levels=c(0,1), ordered = TRUE)
  
  # model
  a <- tryCatch.W.E(
    expr    = { 
      set.seed(s)
      
      levels(df$outcome) <- c("zero", "one")
      
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
          outcome ~., 
          data = df, 
          method = "svmRadial", 
          metric = "Deviance",
          maximize = FALSE, 
          trControl = train_ctrl,
          preProcess = c("center", "scale"), 
          tuneLength = 3
      ) 
      
      pred <- predict(mod, newdata = test, type = "prob")["one"]
      colnames(pred) <- "pred"
      pred <- pred
    }
  )
  
  # sort predictions and errors
  if (is.data.frame(a$value) == FALSE){
    pred  <- NA
    error <- a$value
    problematic <- 1
  } else{
    pred <- a$value
  }
  
  # sort warnings
  if(is.null(a$warning) == FALSE){
    warning <- a$warning
    problematic <- 1
  }
  
  # output 
  pred   <- cbind(pred, "iter" = iter, "class" = as.numeric(test$outcome) - 1, "problematic" = problematic)
  output <- list("pred" = pred, "warning" = warning, "error" = error, "problem" = problematic)
  return(output)
}


do_rnf <- function(df, test, iter){
  
  error       <- 0
  warning     <- 0
  problematic <- 0
  
  # ensure correct factor levels, and character labels 
  df$outcome     = factor(df$outcome, levels=c(0,1), ordered = TRUE)
  
  # model
  a <- tryCatch.W.E(
    expr    = { 
      
      npred = dim(df)[2] -1 
      levels(df$outcome) <- c("zero", "one")
      
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
          outcome ~., 
          data = df, 
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
      
      pred <- predict(mod, newdata = test, type = "prob")["one"]
      colnames(pred) <- "pred"
      pred <- pred
    }
  )
  
  # sort predictions and errors
  if (is.data.frame(a$value) == FALSE){
    pred  <- NA
    error <- a$value
    problematic <- 1
  } else{
    pred <- a$value
  }
  
  # sort warnings
  if(is.null(a$warning) == FALSE){
    warning <- a$warning
    problematic <- 1
  }
  
  # output 
  pred   <- cbind(pred, "iter" = iter, "class" = as.numeric(test$outcome) - 1, "problematic" = problematic)
  output <- list("pred" = pred, "warning" = warning, "error" = error, "problem" = problematic)
  return(output)
}


do_xgb <- function(df, test, iter){
  
  error       <- 0
  warning     <- 0
  problematic <- 0
  
  # ensure correct factor levels, and character labels 
  df$outcome     = factor(df$outcome, levels=c(0,1), ordered = TRUE)
  
  # model
  a <- tryCatch.W.E(
    expr    = { 
      
      levels(df$outcome) <- c("zero", "one")
      
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
          outcome ~., 
          data = df, 
          method = "xgbTree", 
          metric = "Deviance",
          maximize = FALSE, 
          trControl = train_ctrl, 
          verbosity = 0, 
          verbose = FALSE, 
          tuneLength = 3
      )
      
      pred <- predict(mod, newdata = test, type = "prob")["one"]
      colnames(pred) <- "pred"
      pred <- pred
    }
  )
  
  # sort predictions and errors
  if (is.data.frame(a$value) == FALSE){
    pred  <- NA
    error <- a$value
    problematic <- 1
  } else{
    pred <- a$value
  }
  
  # sort warnings
  if(is.null(a$warning) == FALSE){
    warning <- a$warning
    problematic <- 1
  }
  
  # output 
  pred   <- cbind(pred, "iter" = iter, "class" = as.numeric(test$outcome) - 1, "problematic" = problematic)
  output <- list("pred" = pred, "warning" = warning, "error" = error, "problem" = problematic)
  return(output)
}


do_rub <- function(df, test, iter){
  
  error       <- 0
  warning     <- 0
  problematic <- 0
  
  # ensure correct levels of the factors
  df$outcome = factor(df$outcome, levels=c(0,1), ordered = TRUE)
  
  # model 
  a <- tryCatch.W.E(
    expr    = { 
      mod  <- rus(outcome ~., size = 10, alg = "svm", data = df, svm.ker = "radial")
      pred <- predict(mod, newdata = test)
    }
  )
  
  # sort predictions and errors
  if (is.vector(a$value) == FALSE){
    pred  <- NA
    error <- a$value
    problematic <- 1
  } else{
    pred <- a$value
  }
  
  # sort warnings
  if(is.null(a$warning) == FALSE){
    warning <- a$warning
    problematic <- 1
  }
  
  # output 
  pred   <- cbind(pred, "iter" = iter, "class" = as.numeric(test$outcome) - 1, "problematic" = problematic)
  output <- list("pred" = pred, "warning" = warning, "error" = error, "problem" = problematic)
  return(output)
}


do_eas <- function(df, test, iter){
  
  error       <- 0
  warning     <- 0
  problematic <- 0
  
  # model 
  a <- tryCatch.W.E(
    expr    = { 
      train_x <- df %>% dplyr::select(-outcome)
      train_y <- df %>% dplyr::pull(outcome)
      test_x  <- test %>% dplyr::select(-outcome)
      
      mod  <- EasyEnsemble(train_x, train_y)
      pred <- predict(mod, test_x, type = "probability")$X1
    }
  )
  
  # sort predictions and errors
  if (is.vector(a$value) == FALSE){
    pred  <- NA
    error <- a$value
    problematic <- 1
  } else{
    pred <- a$value
  }
  
  # sort warnings
  if(is.null(a$warning) == FALSE){
    warning <- a$warning
    problematic <- 1
  }
  
  # output 
  pred   <- cbind(pred, "iter" = iter, "class" = as.numeric(test$outcome) - 1, "problematic" = problematic)
  output <- list("pred" = pred, "warning" = warning, "error" = error, "problem" = problematic)
  return(output)

}

# performance metrics ----------------------------------------------------------

auroc <- function (probs, outcome){
  if(any(is.na(probs)) == TRUE){
    return(NA)
  }
  
  if(is.numeric(outcome) == FALSE) {outcome <- as.numeric(outcome) - 1}
  # auc  <- as.numeric(pROC::auc(outcome~probs))
  auc  <- as.numeric(pROC::roc(outcome~probs, quiet = TRUE)$auc)
  
  return(auc)
}

scaled_brier_score <- function(probs, outcome){
  if(any(is.na(probs)) == TRUE){
    return(NA)
  }
  
  obs       = as.numeric(outcome)-1
  unscaled  = mean((probs - obs)^2)
  brier_max = mean(obs)*(1 - mean(obs))
  scaled    = 1 - (unscaled / brier_max)
  
  return(scaled)
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
  if(any(is.na(pred)) == FALSE & any(pred < 0) == TRUE){
    pred[pred < 0 ] <- 1e-10
  }
  if(any(is.na(pred)) == FALSE & any(pred > 1) == TRUE){
    pred[pred < 1 ] <- 1- 1e-10
  }
  
  a = tryCatch.W.E(expr = {calibration_intercept(pred, class)})
  b = tryCatch.W.E(expr = {calibration_slope(pred,class)})
  
  return(
    tibble(
      min = min(pred),
      max = max(pred),
      var = var(pred),
      auc = auroc(pred, class),
      bri = scaled_brier_score(pred, class), 
      int = as.character(a$value), 
      slp = as.character(b$value),
      int_warn = ifelse(is.null(a$warning), "0", as.character(a$warning)),
      slp_warn = ifelse(is.null(b$warning), "0", as.character(b$warning))
    )
  )
}

get_stats <- function(dataframe){
    pm <- 
      dataframe %>%
      mutate(iter  = as.factor(iter),
             class = as.factor(class)) %>%
      group_by(iter) %>%
      group_modify(~performance_statistics(pred = .x$pred, class = .x$class))
    
    return(pm)
}

# --------------- simulation set up --------------------------------------------
  
algorithms   <- c("logistic_regression", "support_vector_machine", "random_forest", "xgboost", "rusboost", "easy_ensemble")
corrections  <- c("control", "rus", "ros", "smo", "sen")

pairs <- 
  expand_grid(corrections, algorithms) %>% 
  mutate(pair_id = c(1:30)) %>% 
  relocate(pair_id)

# --------------- simulation code ----------------------------------------------

sim_in_parallel <- function(scenario, pair, start, stop){
  
  # set up 
  sc    = scenario
  p     = pair
  out   = c()
  info  = c()
  
  # save scenario input and correction/algorithm pair 
  input      <- set[[sc]]
  pair_id    <- as.character(pairs[p,1])
  correction <- as.character(pairs[p,2])
  algorithm  <- as.character(pairs[p,3])
  
  # simulation
  for (i in start:stop){
    
    # set seed
    seed <- pull(seed.farm[sc])[i]
    set.seed(seed)
    
    # train and test data
    df <- generate_data(input)
    tf <- generate_data(input, test = TRUE)
    ef <- mean(as.numeric(df$outcome)-1)

    # calculate percent over sample needed to balance
    if (ef < 0.5){
      percO <- (1/ef - 1/0.5)*100
    } else{
      temp_ef <- 1-ef
      percO <- (1/temp_ef - 1/0.5)*100
    }

    # pre-process step 
    if (ef < 0.485 | ef > 0.515){
      if(correction == "rus"){ df <- do_rus(df) }
      if(correction == "ros"){ df <- do_ros(df) }
      if(correction == "smo"){ df <- do_smo(df, percO) }
      if(correction == "sen"){ df <- do_sen(df, percO) }
    }
    
    # implement algorithm 
    if(algorithm == "logistic_regression")   { o <- do_lrg(df, tf, i, seed) }
    if(algorithm == "support_vector_machine"){ o <- do_svm(df, tf, i, seed) }
    if(algorithm == "random_forest")         { o <- do_rnf(df, tf, i) }
    if(algorithm == "xgboost")               { o <- do_xgb(df, tf, i) }
    if(algorithm == "rusboost")              { o <- do_rub(df, tf, i) }
    if(algorithm == "easy_ensemble")         { o <- do_eas(df, tf, i) }
    
    # save iteration information
    iter_info <- 
      cbind("scenario" = sc,
            "pair_id"  = p,
            "iter"     = i,
            "seed"     = seed, 
            "problem"  = o$problem,
            "warning"  = as.character(o$warning), 
            "err"      = as.character(o$error), 
            "sc_ef"    = input$event_frac,
            "obs_ef"   = round(ef, 3),
            "new_ef"   = round(mean(as.numeric(df$outcome)-1), 3)
      )
    
    # store results
    out  <- rbind(out, o$pred)      # predictions
    info <- rbind(info, iter_info)  # iteration info
  }
  
  # save predicted probabilities
  out  <- cbind("scenario" = sc, "pair_id" = p, out)
  out  <- as.data.frame(out)  %>% `rownames<-`( NULL )
  # saveRDS(out, file = paste0("sim_results/predictions/sc", sc, "_p", pair_id, "_", start, "_", stop,".rds"))
  
  # save iteration info 
  info <- as.data.frame(info) %>% `rownames<-`( NULL )
  # saveRDS(info, file = paste0("sim_results/iteration_info/sc", sc, "_p", pair_id, "_", start, "_", stop, ".rds"))
  
  # save per iteration results
  results <- get_stats(out) 
  results <- merge(info, results, by = "iter")
  # saveRDS(results, paste0("sim_results/per_iter_results/sc", sc, "_p", pair_id, "_", start, "_", stop,".rds"))
  
  return(list(out, info, results))
}

# ------------------------ view results ----------------------------------------

calibration_plot <- function(dataframe){
  
  dataframe %>%
    filter(problematic != 1) %>% 
    drop_na() %>% 
    ggplot(aes(x = pred, y = class, group = as.factor(iter))) +
    geom_line(stat = "smooth",
              method = stats::loess, 
              formula = y ~ x,
              se = F, 
              linewidth = 0.05, 
              alpha = 0.1,
              color = "#0d0887") +
    geom_abline(slope = 1, intercept = 0, size = 1) +
    theme_minimal() +
    xlab(" ") + 
    ylab(" ") +
    xlim(0,1) + 
    ylim(0,1)  
}