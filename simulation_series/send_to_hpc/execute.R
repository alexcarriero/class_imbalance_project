# Libraries --------------------------------------------------------------------
library(tidyverse)
library(MASS)
library(caret)
library(pROC)

# Manual Override Package Installation
# ROSE 
source("packages/ROSE-master/R/data_balancing_funcs.R")
source("packages/ROSE-master/R/estimation_funcs.R")
source("packages/ROSE-master/R/ROSE_eval.R")

# IRIC 
library(RANN)
library(foreach)
source("packages/IRIC-master/R/Ensemble-based level/EasyEnsemble.R")
source("packages/IRIC-master/R/Data level/SmoteENN.R")
source("packages/IRIC-master/R/Data level/SMOTE.R")

# ebmc 
source("packages/ebmc-master/R/ebmc.R")

# simsalapar
source("packages/simsalapar-master/R/tryCatchWE.R")


# Helpers ----------------------------------------------------------------------

# functions
source("functions2.R")

# seeds 
seed.farm <- readRDS("seeds/seedfarm.RData")

# simulation scenarios 
set <- readRDS("data-generating-mechanism/set.RData")

# Simulation -------------------------------------------------------------------

slurm_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

start_iter <- c(seq(1911, 1991, by = 10))
stop_iter  <- c(seq(1920, 2000, by = 10))

which_iters <- cbind("start_iter" = start_iter, "stop_iter" = stop_iter)

scenario <- 6
start    <- which_iters[slurm_id, 1]
stop     <- which_iters[slurm_id, 2]

sim_in_series(scenario, start, stop)
