#!/bin/bash
#
#SBATCH --job=plots_r18
#SBATCH --output="./hello_from_hpc_plots_r18.txt"
#
#SBATCH --ntasks=1
#SBATCH --time=01:00:00
#SBATCH --mem-per-cpu=5G
#SBATCH --chdir="./"
#SBATCH --mail-type=ALL
#SBATCH --mail-user=a.j.carriero@umcutrecht.nl
#
#SBATCH --array=1801-2000

Rscript "./execute_plots.R" $SLURM_ARRAY_TASK_ID

