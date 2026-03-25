#!/bin/bash
#SBATCH -A hpc2n2026-137
#SBATCH -t 24:00:00
#SBATCH -N 1
#SBATCH -c 1
#SBATCH -o output_%j.log
#SBATCH -e error_%j.log

ml purge >/dev/null 2>&1
ml GCC/14.2.0 R/4.5.1

cd $SLURM_SUBMIT_DIR

# Force renv to use project library
export RENV_PROJECT="$SLURM_SUBMIT_DIR"
export RENV_PATHS_LIBRARY="$SLURM_SUBMIT_DIR/renv/library"

Rscript scripts/DataGen.R