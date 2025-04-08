##########################################################################
# Fill in your variables here
MIN_CORES=1     # minimum amount of cores to reserve
MAX_CORES=16    # max amount of cores to reserve
STEP_CORES=1    # how much difference between benchmarking points
# PARTITION=genoa  # which HPC partition you want to run on, just the one for now

CURRENT_CORES=$MIN_CORES
# CURRENT_CORES=$(( $CORES_STEP + 16 ))
#########################################################################

# Checks if the required .tpr exists for running the mdrun
#if [ /input_data/benchmark.tpr ]; then
#    
#

while [ $CURRENT_CORES -le $MAX_CORES ]
do
    mkdir -p benchmark_$CURRENT_CORES
    FILENAME=benchmark_$CURRENT_CORES/benchmark.sh

    # we need to remove the old shell file (if it exists). Otherwise we keep appending to it
    rm $FILENAME -f
    touch $FILENAME

    # Makes the shell scripts that do the mdrun
    echo "#!/bin/bash"                                                       >> $FILENAME
    echo "#SBATCH --nodes=1"                                                 >> $FILENAME
    echo "#SBATCH --cpus-per-task=$CURRENT_CORES"                            >> $FILENAME
    echo "#SBATCH --time=00:10:00"                                           >> $FILENAME
    echo "#SBATCH --tasks=1"                                                 >> $FILENAME
    echo "#SBATCH --mem=8G"                                 >> $FILENAME
#SBATCH --error=errors/DE_MHC_new.err
#SBATCH --output=outputs/DE_MHC_new.txt

    echo "#SBATCH --error=benchmark.err"                    >> $FILENAME
    echo "#SBATCH --output=benchmark.txt"                   >> $FILENAME
    echo "module purge"                                     >> $FILENAME
    echo "module load R/4.4.1-gfbf-2023b"                   >> $FILENAME
    echo "Rscript ../CINsim_HPC_benchmark.R"                     >> $FILENAME

    echo "echo "$FILENAME is done!""

    cd benchmark_$CURRENT_CORES
    echo "Would run the sbatch for $CURRENT_CORES now..."
    # sbatch benchmark.sh
    cd ..
    CURRENT_CORES=$(( $CURRENT_CORES + $STEP_CORES ))

done
