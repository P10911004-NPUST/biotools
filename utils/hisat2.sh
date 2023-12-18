#!/bin/bash

#SBATCH --output slurm-%j.out
#SBATCH --job-name=hisat2
#SBATCH --nodes=1
#SBATCH --mincpus=40
#SBATCH --mem=400G
#SBATCH --time=7-23:59:59

module load R hisat2 samtools

cd /bcst/YAMADA/jklai

srun Rscript ./rcode/hisat2.R