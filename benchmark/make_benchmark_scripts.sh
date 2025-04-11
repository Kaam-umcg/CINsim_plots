##########################################################################
# Fill in your variables here
MIN_CORES=2     # minimum amount of cores to reserve
MAX_CORES=16    # max amount of cores to reserve
STEP_CORES=2    # how much difference between benchmarking points

CURRENT_CORES=$MIN_CORES
# CURRENT_CORES=$(( $CORES_STEP + 16 ))
#########################################################################

mkdir -p results/max_g
mkdir -p results/max_cells
mkdir -p results/iterations

while [ $CURRENT_CORES -le $MAX_CORES ]
do
    mkdir -p benchmark_$CURRENT_CORES
    FILENAME=benchmark_$CURRENT_CORES/benchmark.sh

    # We need to remove the old shell file (if it exists). 
    # Otherwise we keep appending to it
    rm $FILENAME -f
    touch $FILENAME

    # Makes the shell scripts that run the benchmark
    echo "#!/bin/bash"                                      >> $FILENAME
    echo "#SBATCH --nodes=1"                                >> $FILENAME
    echo "#SBATCH --cpus-per-task=$CURRENT_CORES"           >> $FILENAME
    echo "#SBATCH --time=04:00:00"                          >> $FILENAME
    echo "#SBATCH --tasks=1"                                >> $FILENAME
    echo "#SBATCH --mem=4G"                                 >> $FILENAME
    echo "#SBATCH --error=benchmark.err"                    >> $FILENAME
    echo "#SBATCH --output=benchmark.txt"                   >> $FILENAME

    echo "module purge"                                     >> $FILENAME
    echo "module load R/4.4.1-gfbf-2023b"                   >> $FILENAME
    echo "Rscript ../CINsim_HPC_benchmark.R"                >> $FILENAME

    echo "$FILENAME is done!"

    cd benchmark_$CURRENT_CORES
    echo "Starting benchmark for $CURRENT_CORES now..."
    sbatch benchmark.sh
    cd ..
    CURRENT_CORES=$(( $CURRENT_CORES + $STEP_CORES ))
done
