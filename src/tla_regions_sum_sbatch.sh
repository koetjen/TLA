#!/bin/bash -l

#SBATCH --job-name=TLA_setup_sum
#SBATCH -N 1                         # number of compute nodes
#SBATCH -n 1                         # number of concurrently running jobs
#SBATCH --mem=10G                    # requests memory on the node
#SBATCH -p general
#SBATCH -q public
#SBATCH --time=0-0:05:00             # upper bound time: d-hh:mm:ss
#SBATCH --export=NONE                # purge the job-submitting shell environment
#SBATCH --output=log/slurm-%A.%a.out # stdout file
#SBATCH --error=log/slurm-%A.%a.err  # stderr file

# Load the python interpreter
module load mamba/latest
source activate tlaenv

# summarize study statistics 
python src/tla_regions_sum.py ${STUDY}
