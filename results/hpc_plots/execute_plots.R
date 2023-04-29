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

source("functions.R")
source("process_results.R")

# Generate Plots ---------------------------------------------------------------

i   <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

scN <- "sc18"

# save_plot_coords(i, scN, recalibrated = F)

save_plot_coords(i, scN, recalibrated = T)

