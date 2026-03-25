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

# Run with renv fully activated
Rscript -e 'setwd(Sys.getenv("SLURM_SUBMIT_DIR")); source("renv/activate.R"); renv::load(); source("scripts/DataGen.R")'