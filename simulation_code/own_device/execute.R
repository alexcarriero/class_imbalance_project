# Libraries --------------------------------------------------------------------
library(tidyverse)
library(MASS)
library(caret)
library(pROC)

# setwd("PLEASE SET YOUR WORKING DIRECTORY")

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
source("functions.R")

# seeds 
seed.farm <- readRDS("seeds/seedfarm.RData")

# simulation scenarios 
set <- readRDS("data-generating-mechanism/set.RData")

# Simulation -------------------------------------------------------------------

# select scenario
scenario <- 1

# select iterations to run 
start    <- 1
stop     <- 2000

sim_in_series(scenario, start, stop)
