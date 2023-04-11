# sample seeds 

df <- c()
set.seed(1)

for(i in 1:27){
  sd <- sample(1e7, 2000)
  df <- cbind(df, sd)
}

colnames(df) <- c(paste0("sc", 1:27))
seed.farm    <- as.data.frame(df)  %>% `rownames<-`( NULL )
saveRDS(seed.farm, "seedfarm.RData")


# check if all unique
sum(duplicated(seed.farm)) # all unique seeds 