#!/bin/bash

# This script reruns all the shells present in this directory
# I manually update this, as I prefer giving a last check before eating a ton
# of computation time.

####
# MAIN FIGURE SCRIPTS
####
sbatch Figure_1.sh
sbatch Figure_2.sh

# split this into two for performance and parallelisation reasons
sbatch Figure_3A.sh
sbatch Figure_3C.sh

# 3E needs to have results from a few other simulations, so you have to wait with
# starting that script till later in the worflow
# sbatch Figure_3E.sh
sbatch Figure_4.sh

sbatch Figure_5.sh
sbatch Figure_S1.sh
sbatch Figure_S2.sh

# no need for scripting for S3, as it is from the scWGS data
# Figure S4a needs a lot of hyperparam tuning, so we split it into two:
# diploid and WGD (tetraploid) versions of the same script
sbatch Figure_S4a_diploid.sh
sbatch Figure_S4a_WGD.sh
sbatch Figure S4c.sh

# fig 5 needs a lot of hyperparam searching that is done for
# supplemental figure S4, so we do it at the end
sbatch Figure_5.sh

# since 3E needs results from all the other scripts before, 
# we can't realistically run it parallel with the other scripts
sbatch Figure_3E.sh 