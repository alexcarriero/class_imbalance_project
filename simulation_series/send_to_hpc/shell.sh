#!/bin/bash
#
#SBATCH --job=sc6
#SBATCH --output="./hello_from_hpc6.txt"
#
#SBATCH --ntasks=1
#SBATCH --time=35:00:00
#SBATCH --mem-per-cpu=10G
#SBATCH --chdir="./"
#SBATCH --mail-type=ALL
#SBATCH --mail-user=a.j.carriero@umcutrecht.nl
#
#SBATCH --array=1-9

Rscript "./execute.R" $SLURM_ARRAY_TASK_ID

