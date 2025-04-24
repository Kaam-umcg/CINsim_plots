## Computational Benchmarking for CINsim
This folder contains all required code for benchmarking the computational performance of CINsim. This directory is split into a few main components:

`CINsim_HPC_benchmark.R` is an R script that performs some benchmarking for CINsim

`make_benchmark_scripts.sh` is the main meat of the benchmarking - it creates and invokes .sh files that run `CINsim_HPC_benchmark.R` with a variety of cores. 
It does this through creation of directories with the pattern `benchmark_#cores`, where `#cores` is the number of cores for that benchmarking run. This script can be dangerous, as 
it invokes `sbatch` within an `sbatch` - this is something that would signal deep red flags for any sysadmin. For my HPC cluster this is allowed - yours might be more
strict with recursive sbatch calls. 

#### Make sure you properly understand all the options in the `make_benchmark_scripts.sh` file before running it!

The other scripts `plot_benchmark.R/.sh` are simple scripts that make the figures related to benchmarking for the [CINsim paper, BioRxiv](https://www.biorxiv.org/content/10.1101/2023.02.14.528596v1).
