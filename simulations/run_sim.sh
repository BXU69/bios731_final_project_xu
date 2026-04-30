#!/bin/bash
#SBATCH --array=1-18
#SBATCH --job-name=lcr_sim
#SBATCH --partition=wrobel
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=UNLIMITED
#SBATCH --output=logs/sim_%A_%a.out
#SBATCH --error=logs/sim_%A_%a.err

module purge
module load R

# Rscript to run an r script
# This stores which job is running (1, 2, 3, etc)
JOBID=$SLURM_ARRAY_TASK_ID
Rscript simulations/run_sim.R $JOBID

