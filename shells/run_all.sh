# This script reruns all the shells present in this directory
# I manually update this, as I prefer giving a last check before eating a ton
# of computation time.

# sets all the scripts for the main figures
sbatch Figure_1.sh
sbatch Figure_2.sh
sbatch Figure_3A.sh
sbatch Figure_3C.sh

# since 3E needs results from all the other scripts before, 
# we can't realistically run it parallel with the other scripts
# sbatch Figure_3E.sh 

sbatch Figure_4.sh

# fig 5 needs a lot of hyperparam searching, so we sbatch it as such
sbatch Figure_5_14T.sh
sbatch Figure_5_16T.sh
sbatch Figure_5_9T.sh
sbatch Figure_5_24TB.sh


# and all the scripts for the supplemental figures
