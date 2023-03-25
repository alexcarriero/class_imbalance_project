#!/bin/bash
#
#SBATCH --job-name=sc3(1401-2k)
#SBATCH --output="./hello_from_hpc.txt"
#
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=20G
#SBATCH --chdir="./"
#SBATCH --mail-type=ALL
#SBATCH --mail-user=a.j.carriero@umcutrecht.nl
#
#SBATCH --array=1-4,6-10,12-16,18-22,24-28,30

Rscript "./execute.R" $SLURM_ARRAY_TASK_ID

