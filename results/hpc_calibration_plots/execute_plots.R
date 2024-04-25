# Libraries --------------------------------------------------------------------
library(tidyverse)
library(MASS)
library(caret)
library(pROC)

# simsalapar
source("packages/simsalapar-master/R/tryCatchWE.R")

# Helpers ----------------------------------------------------------------------
source("functions.R")

# Generate Plots ---------------------------------------------------------------

i   <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

scN <- "sc6"

save_plot_coords(i, scN, recalibrated = F)

# save_plot_coords(i, scN, recalibrated = T)