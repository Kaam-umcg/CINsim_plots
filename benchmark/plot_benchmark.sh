#!/bin/bash

#SBATCH --job-name=Benchmark_plot
#SBATCH --time=00:10:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4GB
#SBATCH --error=errors/plot_benchmark.err
#SBATCH --output=outputs/plot_benchmark.log

######
# CHANGE THE NAME OF THE JOB TO THE SCRIPT NAME - OTHERWISE A LOT OF THINGS BREAK
######
START_TIME=$(date -d @$SLURM_JOB_START_TIME)
END_TIME=$(date -d @$SLURM_JOB_END_TIME)
echo $SLURM_JOB_NAME starting at $START_TIME with prospective end time of $END_TIME

# copy the script to the node
cp "plot_benchmark.R" $TMPDIR 

# modules for running scripts
module purge
module load R/4.4.1-gfbf-2023b

# run your script
Rscript $TMPDIR/plot_benchmark.R

# copies all relevant results to SCRATCH
mkdir -p $SCRATCH/CINsim/benchmark
cp $TMPDIR $SCRATCH/CINsim/benchmark -r