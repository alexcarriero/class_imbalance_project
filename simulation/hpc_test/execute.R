
library(tidyverse)
library(MASS)

source("./functions.R")

# scenario number 
sc <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

# seed
sd <- seeds[sc]

# generate data
set.seed(sd)
generate_all_the_data(sc, iter = 2000)
