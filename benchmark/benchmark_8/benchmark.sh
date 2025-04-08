#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --time=00:10:00
#SBATCH --tasks=1
#SBATCH --mem=8G
#SBATCH --error=benchmark.err
#SBATCH --output=benchmark.txt
module purge
module load R/4.4.1-gfbf-2023b
Rscript ../CINsim_HPC_benchmark.R
