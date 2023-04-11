#!/bin/bash
#
#SBATCH --job=TESTA
#SBATCH --output="./hello_from_hpc.txt"
#
#SBATCH --ntasks=1
#SBATCH --time=3:00:00
#SBATCH --mem-per-cpu=20G
#SBATCH --chdir="./"
#SBATCH --mail-type=ALL
#SBATCH --mail-user=a.j.carriero@umcutrecht.nl
#
#SBATCH --array=1-30

Rscript "./execute.R" $SLURM_ARRAY_TASK_ID

