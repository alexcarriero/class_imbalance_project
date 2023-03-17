#!/bin/bash
#
#SBATCH --job-name=test-generate-data
#SBATCH --output="./generated_data/hello_from_hpc.txt"
#
#SBATCH --ntasks=1
#SBATCH --time=02:00:00
#SBATCH --mem-per-cpu=20G
#SBATCH --chdir="./"
#SBATCH --mail-type=ALL
#SBATCH --mail-user=a.j.carriero@umcutrecht.nl
#
#SBATCH --array=1-6

Rscript "./execute.R" $SLURM_ARRAY_TASK_ID

