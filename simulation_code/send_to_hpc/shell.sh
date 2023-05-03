#!/bin/bash
#
#SBATCH --job=sc18
#SBATCH --output="./hello_from_hpc_sc18.txt"
#
#SBATCH --ntasks=1
#SBATCH --time=25:00:00
#SBATCH --mem-per-cpu=10G
#SBATCH --chdir="./"
#SBATCH --mail-type=ALL
#SBATCH --mail-user=a.j.carriero@umcutrecht.nl
#
#SBATCH --array=1701-2000

Rscript "./execute.R" $SLURM_ARRAY_TASK_ID

