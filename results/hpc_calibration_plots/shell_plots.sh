#!/bin/bash
#
#SBATCH --job=6plots
#SBATCH --output="./hello_from_hpc_6plots.txt"
#
#SBATCH --ntasks=1
#SBATCH --time=15:00:00
#SBATCH --mem-per-cpu=10G
#SBATCH --chdir="./"
#SBATCH --mail-type=ALL
#SBATCH --mail-user=a.j.carriero@umcutrecht.nl
#
#SBATCH --array=1000-2000

Rscript "./execute_plots.R" $SLURM_ARRAY_TASK_ID

