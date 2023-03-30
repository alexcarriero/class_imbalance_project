# Libraries --------------------------------------------------------------------
library(tidyverse)
library(MASS)
library(caret)
library(pROC)

setwd("/Users/alexcarriero/Desktop/thesis github/masters_thesis/simulation_series/")

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

# slurm_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

scenario <- 3
start    <- 1200 # which_iters[slurm_id,1]
stop     <- 1200 # which_iters[slurm_id,2]

sim_in_series(scenario, start, stop)
