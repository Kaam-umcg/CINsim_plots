#!/bin/bash

#SBATCH --job-name=no_evolution_low_misseg
#SBATCH --time=00:10:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4GB
#SBATCH --error=errors/no_evolution_low_misseg.err
#SBATCH --output=outputs/no_evolution_low_misseg.log

######
# CHANGE THE NAME OF THE JOB TO THE SCRIPT NAME - OTHERWISE A LOT OF THINGS BREAK
######
START_TIME=$(date -d @$SLURM_JOB_START_TIME)
END_TIME=$(date -d @$SLURM_JOB_END_TIME)
echo $SLURM_JOB_NAME starting at $START_TIME with prospective end time of $END_TIME

# copy the script to the node
cp "../scripts/no_evolution_low_misseg.R" $TMPDIR 

# modules for running scripts
module purge
module load R/4.4.1-gfbf-2023b

# the directory we'll be loading the data from
DATA_SOURCE="/scratch/p319788/CINsim/Figure_2/tmp/results/T_ALL_params"
export DATA_SOURCE

# run your script
Rscript $TMPDIR/no_evolution_low_misseg.R

# copies all relevant results to SCRATCH
mkdir -p $SCRATCH/CINsim/no_evolution_low_misseg
cp $TMPDIR $SCRATCH/CINsim/no_evolution_low_misseg -r