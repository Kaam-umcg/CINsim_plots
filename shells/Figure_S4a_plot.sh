#!/bin/bash

#SBATCH --job-name=Fig_S4a_plot_10
#SBATCH --time=00:10:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=2GB
#SBATCH --error=errors/Figure_S4a_plot_10.err
#SBATCH --output=outputs/Figure_S4a_plot_10.log

######
# CHANGE THE NAME OF THE JOB TO THE SCRIPT NAME - OTHERWISE A LOT OF THINGS BREAK
######
START_TIME=$(date -d @$SLURM_JOB_START_TIME)
END_TIME=$(date -d @$SLURM_JOB_END_TIME)
echo $SLURM_JOB_NAME starting at $START_TIME with prospective end time of $END_TIME

# copy the script to the node
cp "../scripts/" $TMPDIR -r
cp "../data/" $TMPDIR -r

# modules for running scripts
module purge
module load R/4.4.1-gfbf-2023b

export "BASE_PATH"="/scratch/p319788/CINsim"
export "DIPLOID_PATH"="Figure_S4a_diploid/tmp/results"
export "WGD_PATH"="Figure_S4a_WGD/tmp/results"
export "OUTPUT_PATH"="/scratch/p319788/CINsim/Figure_S4a_10_cells/"

# run your script
Rscript $TMPDIR/scripts/Figure_S4a_plot.R

# copies all relevant results to SCRATCH
mkdir -p $OUTPUT_PATH
cp $TMPDIR $OUTPUT_PATH -r