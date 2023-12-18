#!/bin/bash

#SBATCH --output slurm-%j.out
#SBATCH --job-name=subread
#SBATCH --nodes=1
#SBATCH --mincpus=40
#SBATCH --mem=400G
#SBATCH --time=7-23:59:59

module load R subread

cd /bcst/YAMADA/jklai

srun Rscript ./rcode/subread.R