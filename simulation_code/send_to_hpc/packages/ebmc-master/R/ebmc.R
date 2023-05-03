##### Set new class #####
setClass("modelBag", representation = "list")
setClass("modelBst", representation = "list")

##### Funtion creation #####
# Random under-sampling
.ru <- function(target, data, ir = 1)    # ir = Imbalance Ratio. (how many times majority instances are over minority instances)
{
  p <- data[which(data[ ,target] == "1"), ]
  n <- data[which(data[ ,target] == "0"), ]
  n <- n[sample(nrow(n), nrow(p) * ir, replace = TRUE), ]
  result <- rbind(p, n)
  return(result)
}


# Weight update/ pseudo-loss calculation for AdaBoost.M2
.wt.update <- function(probability, prediction, actual, wt, smooth)
{
  fp <- which(ifelse(prediction == "1" & actual == "0", TRUE, FALSE) == TRUE)
  fn <- which(ifelse(prediction == "0" & actual == "1", TRUE, FALSE) == TRUE)
  p_loss <- 0.5 * sum( wt[fp] * (1 - probability[fp, ][ ,"0"] + probability[fp, ][ ,"1"]),  # pseudo-loss
                       wt[fn] * (1 - probability[fn, ][ ,"1"] + probability[fn, ][ ,"0"]) )
  a <- (p_loss + smooth) / (1 - p_loss + smooth) # weight updater with prediction smoothing, dealing with a == 0
  wt[c(fp, fn)] <- rep(1/(length(fp) + length(fn)), (length(fp) + length(fn)))
  wt[fn] <- wt[fn] * a^(0.5 * (1 + probability[fn, ][ ,"1"] - probability[fn, ][ ,"0"]))
  wt[fp] <- wt[fp] * a^(0.5 * (1 + probability[fp, ][ ,"0"] - probability[fp, ][ ,"1"]))
  wt <- wt / sum(wt)
  result <- list()
  result[[1]] <- wt
  result[[2]] <- a
  return(result)
}


# AdaBoost.M2
adam2 <- function(formula, data, size, alg, rf.ntree = 50, svm.ker = "radial")
{
  target <- gsub(" ", "", unlist(strsplit(format(formula), split = "~"))[1])
  data[ ,target] <- as.factor(data[ ,target])

  list_model <- list()
  a <- 0
  w <- rep(1/nrow(data), nrow(data))
  label <- data[ ,target]
  for(i in 1:size)
  {
    if(i == 1) { samp <- sample(nrow(data), nrow(data), replace = FALSE) } # no replacement in first iteration
    else if(i != 1) { samp <- sample(nrow(data), nrow(data), replace = TRUE, prob = w) }
    train <- data[samp, ]
    if(alg == "svm") {
      list_model[[i]] <- e1071::svm(formula, data = train, kernel = svm.ker, probability = TRUE)
      prob <- as.data.frame(attr(predict(list_model[[i]], data, probability = TRUE), "prob"))
    }
    else if(alg == "cart") {
      list_model[[i]] <- rpart::rpart(formula, data = train)
      prob <- as.data.frame(predict(list_model[[i]], data, type = "prob"))
    }
    else if(alg == "c50") {
      list_model[[i]] <- C50::C5.0(formula, data = train)
      prob <- as.data.frame(predict(list_model[[i]], data, type = "prob"))
    }
    else if(alg == "nb") {
      list_model[[i]] <- e1071::naiveBayes(formula, data = train)
      prob <- as.data.frame(predict(list_model[[i]], data, type = "raw"))
    }
    else if(alg == "rf") {
      list_model[[i]] <- randomForest::randomForest(formula, data = train, ntree = rf.ntree)
      prob <- as.data.frame(predict(list_model[[i]], data, type = "prob"))
    }
    pred <- as.factor(ifelse(prob[ ,"1"] >= 0.5, 1, 0))
    new <- .wt.update(probability = prob, prediction = pred, actual = label, wt = w, smooth = 1/nrow(data))
    w <- new[[1]]
    a[i] <- new[[2]]
  }
  result <- list(weakLearners = list_model, errorEstimation = a)
  attr(result, "class") <- "modelBst"
  return(result)
}


# RUSBoost
rus <- function(formula, data, size, alg, ir = 1, rf.ntree = 50, svm.ker = "radial")
{
  target <- gsub(" ", "", unlist(strsplit(format(formula), split = "~"))[1])
  data[ ,target] <- as.factor(data[ ,target])

  list_model <- list()
  a <- 0
  data$w <- rep(1/nrow(data), nrow(data))
  label <- data[ ,target]
  for(i in 1:size)
  {
    train <- .ru(target, data, ir)
    train$w <- train$w / sum(train$w) # normalize sample weights
    train <- train[sample(nrow(train), nrow(train), replace = TRUE, prob = train$w), ] # equivalent to pass w' to learner
    train$w <- NULL # remove weight otherwise it will be used as a variable in when training
    if(alg == "svm") {
      list_model[[i]] <- e1071::svm(formula, data = train, kernel = svm.ker, probability = TRUE)
      prob <- as.data.frame(attr(predict(list_model[[i]], data, probability = TRUE), "prob"))
    }
    else if(alg == "cart") {
      list_model[[i]] <- rpart::rpart(formula, data = train)
      prob <- as.data.frame(predict(list_model[[i]], data, type = "prob"))
    }
    else if(alg == "c50") {
      list_model[[i]] <- C50::C5.0(formula, data = train)
      prob <- as.data.frame(predict(list_model[[i]], data, type = "prob"))
    }
    else if(alg == "nb") {
      list_model[[i]] <- e1071::naiveBayes(formula, data = train)
      prob <- as.data.frame(predict(list_model[[i]], data, type = "raw"))
    }
    else if(alg == "rf") {
      list_model[[i]] <- randomForest::randomForest(formula, data = train, ntree = rf.ntree)
      prob <- as.data.frame(predict(list_model[[i]], data, type = "prob"))
    }
    pred <- as.factor(ifelse(prob[ ,"1"] >= 0.5, 1, 0))
    new <- .wt.update(probability = prob, prediction = pred, actual = label, wt = data$w, smooth = 1/nrow(train))
    data$w <- new[[1]]
    a[i] <- new[[2]]
  }
  result <- list(weakLearners = list_model, errorEstimation = a)
  attr(result, "class") <- "modelBst"
  return(result)
}


# SMOTEBoost
sbo <- function(formula, data, size, alg, over = 100, smote.k = 5, rf.ntree = 50, svm.ker = "radial")
{
  target <- gsub(" ", "", unlist(strsplit(format(formula), split = "~"))[1])
  data[ ,target] <- as.factor(data[ ,target])

  list_model <- list()
  a <- 0
  n <- data[which(data[ ,target] == "0"), ]
  p <- data[which(data[ ,target] == "1"), ]
  data$w <- rep(1/nrow(data), nrow(data))
  label <- data[ ,target]
  for(i in 1:size)
  {
    n <- data[which(data[ ,target] == "0"), ]
    smote <- smotefamily::SMOTE(X = data[ ,colnames(data) != target], target = data[ ,target], K = smote.k, dup_size = round(over/100))
    smote <- rbind(smote$orig_P, smote$syn_dat)
    colnames(smote)[which(colnames(smote) == "class")] <- target # correct column name 'class' bug in package 'smotefamily'

    train <- rbind(n, smote)
    train$w <- train$w / sum(train$w) # normalize sample weights
    train <- train[sample(nrow(train), nrow(train), replace = TRUE, prob = train$w), ] # equivalent to pass w' to learner
    train$w <- NULL # remove weight otherwise it will be used as a variable in when training
    if(alg == "svm") {
      list_model[[i]] <- e1071::svm(formula, data = train, kernel = svm.ker, probability = TRUE)
      prob <- as.data.frame(attr(predict(list_model[[i]], data, probability = TRUE), "prob"))
    }
    else if(alg == "cart") {
      list_model[[i]] <- rpart::rpart(formula, data = train)
      prob <- as.data.frame(predict(list_model[[i]], data, type = "prob"))
    }
    else if(alg == "c50") {
      list_model[[i]] <- C50::C5.0(formula, data = train)
      prob <- as.data.frame(predict(list_model[[i]], data, type = "prob"))
    }
    else if(alg == "nb") {
      list_model[[i]] <- e1071::naiveBayes(formula, data = train)
      prob <- as.data.frame(predict(list_model[[i]], data, type = "raw"))
    }
    else if(alg == "rf") {
      list_model[[i]] <- randomForest::randomForest(formula, data = train, ntree = rf.ntree)
      prob <- as.data.frame(predict(list_model[[i]], data, type = "prob"))
    }
    pred <- as.factor(ifelse(prob[ ,"1"] >= 0.5, 1, 0))
    new <- .wt.update(probability = prob, prediction = pred, actual = label, wt = data$w, smooth = 1/nrow(data))
    data$w <- new[[1]]
    a[i] <- new[[2]]
  }
  result <- list(weakLearners = list_model, errorEstimation = a)
  attr(result, "class") <- "modelBst"
  return(result)
}


# SMOTEBagging
sbag <- function(formula, data, size, alg, smote.k = 5, rf.ntree = 50, svm.ker = "radial")
{
  target <- gsub(" ", "", unlist(strsplit(format(formula), split = "~"))[1])
  data[ ,target] <- as.factor(data[ ,target])

  list_train <- list()
  list_model <- list()
  b <- rep(seq(10, 100, 10), (size %/% 10 + 1))
  b <- b[1:size]  # b% detemining proportions of random over-sampling and SMOTE
  n <- data[which(data[ ,target] == "0"), ]
  p <- data[which(data[ ,target] == "1"), ]
  smote <- smotefamily::SMOTE(X = data[ ,colnames(data) != target], target = data[ ,target], K = smote.k, dup_size = round(nrow(n)/nrow(p)*965/100))
  smote <- smote$syn_data
  colnames(smote)[which(colnames(smote) == "class")] <- target # correct column name 'class' bug in package 'smotefamily'

  for(i in 1:size) {
    resamp_rate <- nrow(n)/nrow(p) * (b[i] / 100)
    resamp <- p[sample(nrow(p), round(nrow(p) * resamp_rate), replace = TRUE), ] # random over-sampling
    smt <- smote[sample(nrow(smote), (nrow(n) - nrow(resamp)), replace = TRUE), ] # SMOTE
    list_train[[i]] <- rbind(n, resamp, smt) # create training set
  }
  for(i in 1:size) {
    if(alg == "svm") {
      list_model[[i]] <- e1071::svm(formula, data = list_train[[i]], kernel = svm.ker, probability = TRUE)
    }
    else if(alg == "cart") {
      list_model[[i]] <- rpart::rpart(formula, data = list_train[[i]])
    }
    else if(alg == "c50") {
      list_model[[i]] <- C50::C5.0(formula, data = list_train[[i]])
    }
    else if(alg == "nb") {
      list_model[[i]] <- e1071::naiveBayes(formula, data = list_train[[i]])
    }
    else if(alg == "rf") {
      list_model[[i]] <- randomForest::randomForest(formula, data = list_train[[i]], ntree = rf.ntree)
    }
  }
  attr(list_model, "class") <- "modelBag"
  return(list_model)
}


# UnderBagging
ub <- function(formula, data, size, alg, ir = 1, rf.ntree = 50, svm.ker = "radial")
{
  target <- gsub(" ", "", unlist(strsplit(format(formula), split = "~"))[1])
  data[ ,target] <- as.factor(data[ ,target])

  list_train <- list()
  list_model <- list()
  for(i in 1:size) {
    list_train[[i]] <- .ru(target, data, ir)
  }
  for(i in 1:size) {
    if(alg == "svm") {
      list_model[[i]] <- e1071::svm(formula, data = list_train[[i]], kernel = svm.ker, probability = TRUE)
    }
    else if(alg == "cart") {
      list_model[[i]] <- rpart::rpart(formula, data = list_train[[i]])
    }
    else if(alg == "c50") {
      list_model[[i]] <- C50::C5.0(formula, data = list_train[[i]])
    }
    else if(alg == "nb") {
      list_model[[i]] <- e1071::naiveBayes(formula, data = list_train[[i]])
    }
    else if(alg == "rf") {
      list_model[[i]] <- randomForest::randomForest(formula, data = list_train[[i]], ntree = rf.ntree)
    }
  }
  attr(list_model, "class") <- "modelBag"
  return(list_model)
}


# Prediction for Boosting-based method
predict.modelBst <- function(object, newdata, type = "prob", ...)
{
  list_model <- object[[1]]
  a <- object[[2]]
  a <- log(1/a, base = exp(1)) / sum(log(1/a, base = exp(1))) # normalize alpha values into percentage
  if(attr(list_model[[1]], "class")[2] %in% "svm") {
    prob <- lapply(lapply(list_model, predict, newdata, probability = TRUE), attr, which = "probabilities")
    prob <- lapply(prob, subset, select = "1")
  }
  else if(attr(list_model[[1]], "class")[1] == "rpart") {
    prob <- lapply(lapply(list_model, predict, newdata, type = "prob"), subset, select = "1")
  }
  else if(attr(list_model[[1]], "class")[1] == "C5.0") {
    prob <- lapply(lapply(list_model, predict, newdata, type = "prob"), subset, select = "1")
  }
  else if(attr(list_model[[1]], "class")[1] == "naiveBayes") {
    prob <- lapply(lapply(list_model, predict, newdata, type = "raw"), subset, select = "1")
  }
  else if(attr(list_model[[1]], "class")[2] == "randomForest") {
    prob <- lapply(lapply(list_model, predict, newdata, type = "prob"), subset, select = "1")
  }
  prob <- rowSums(mapply("*", prob, a))
  if(type == "class") {
    pred <- as.factor(ifelse(prob > 0.5, 1, 0))
    return(pred)
  }
  else if(type == "prob") { return(prob) }
}


# Prediction for Bagging-based method
predict.modelBag <- function(object, newdata, type = "prob", ...)
{
  a <- rep(1/length(object), length(object)) # voting weight
  if(attr(object[[1]], "class")[2] %in% "svm") {
    prob <- lapply(lapply(object, predict, newdata, probability = TRUE), attr, which = "probabilities")
    prob <- lapply(prob, subset, select = "1")
  }
  else if(attr(object[[1]], "class")[1] == "rpart") {
    prob <- lapply(lapply(object, predict, newdata, type = "prob"), subset, select = "1")
  }
  else if(attr(object[[1]], "class")[1] == "C5.0") {
    prob <- lapply(lapply(object, predict, newdata, type = "prob"), subset, select = "1")
  }
  else if(attr(object[[1]], "class")[1] == "naiveBayes") {
    prob <- lapply(lapply(object, predict, newdata, type = "raw"), subset, select = "1")
  }
  else if(attr(object[[1]], "class")[2] == "randomForest") {
    prob <- lapply(lapply(object, predict, newdata, type = "prob"), subset, select = "1")
  }
  prob <- rowSums(mapply("*", prob, a))
  if(type == "class") {
    pred <- as.factor(ifelse(prob > 0.5, 1, 0))
    return(pred)
  }
  else if(type == "prob") {
    return(prob)
  }
}


# Performance measurement calculation
measure <- function(label, probability, metric, threshold = 0.5)
{
  if(metric == "auc") {
    result <- pROC::auc(pROC::roc(response = label, predictor = probability))[1]
  }
  else if(metric %in% c("acc", "tpr", "tnr", "gmean", "f")) {
    pred <- as.factor(ifelse(probability >= threshold, 1, 0))
    pred <- factor(pred, levels = c("0", "1"))
    label <- factor(label, levels = c("0", "1"))
    cm <- table(label, pred)
    if(metric == "acc") {
      result <- mean(label == pred)
    }
    else if(metric == "tpr"){
      result <- cm["1", "1"] / sum(cm["1", ])
    }
    else if(metric == ("tnr")) {
      result <- cm["0", "0"] / sum(cm["0", ])
    }
    else if(metric == "gmean") {
      tpr <- cm["1", "1"] / sum(cm["1", ])
      tnr <- cm["0", "0"] / sum(cm["0", ])
      result <- (tpr * tnr)^0.5
    }
    else if(metric == "f") {
      precision <- cm["1", "1"] / sum(cm[ ,"1"])
      recall <- cm["1", "1"] / sum(cm["1", ])
      result <- (2 * precision * recall) / (precision + recall)
    }
  }
  return(result)
}
