
seeds <- c(940, 615, 963, 604, 335, 505, 607, 85, 965, 529, 937, 454, 519, 226,
           591, 860, 768, 120, 911, 231, 368, 788, 896, 778, 552, 419, 103)

generate_all_the_data <- function(sc, iter = 2000){
  
  set <- readRDS("./set.RData")
  
  scenario      <- set[[sc]]

  storage_data  <- vector(mode = "list", length = iter)
  storage_test  <- vector(mode = "list", length = iter)
  
  for (i in 1:iter){
    storage_data[[i]] <- generate_data(scenario)
    storage_test[[i]] <- generate_data(scenario, test = TRUE)
  }
  
  out <- list("data" = storage_data, "test" = storage_test)
  saveRDS(out, file = paste0("./generated_data/scenario_", sc,".rds"))
}


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
  
  return(df)
}
