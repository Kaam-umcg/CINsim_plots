# TODO - update!

## Code required for reproducing CINsim plots
This repository contains all required code for reproducing every plot as found in the [CINsim paper, BioRxiv](https://www.biorxiv.org/content/10.1101/2023.02.14.528596v1). 
The scripts are designed to run on HPC systems using SLURM - running everything locally would require a very strong desktop (>10 fast cores) and a lot of time (>2 weeks). 
Assuming that is not the case, the instructions in shells/run_all.sh can be followed to run every shell scripts that invokes the associated .R script. 
There are some hardcoded paths that require special consideration to ensure directories exist - I've opened an issue in this repository if I ever get the time to fix those issues.

### Directory information
`scripts` contains all .R files used for running simulations and creating plots when relevant. 

`shells` contains all .sh files used to `sbatch` jobs into a SLURM queue. The relevant SLURM options _should_ be enough to complete every script with minimal inefficiency, but your mileage may vary depending on CPU speed
of your specific HPC cluster and any updates that may have happened to CINsim. Since some simulations fail and SNOW maps 1 core to a simulation, some idle CPU time is unavoiable in the current setup.

`plots` contains all the plots as output from the scripts.

`data` contains data needed to run scripts correctly - this is mostly copy number frequency objects that come from scWGS data that would be outside the scope of this repository to include in the analysis.

`benchmark` is a directory specific to the benchmarking done to gage the computational performance of CINsim. It is excluded from the standard structure of the rest of the repository to ensure some code features
such as `list.dirs()` work correctly. It also contains its own README for further information on how to reproduce the plots related to benchmarking.
