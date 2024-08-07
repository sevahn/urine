#!/bin/bash
#SBATCH --job-name=all_sed_supt.sh
#SBATCH --output=all_sed_supt.sh.out
#SBATCH --error=all_sed_supt.sh.err
#SBATCH --qos=normal
#SBATCH --ntasks=1
#SBATCH --partition=quake,owners,normal
#SBATCH --cpus-per-task=1
#SBATCH --mem=10000
#SBATCH --time=8:00:00


source activate snakemake

python3 bootstrap_ci.py all_sed_supt 20231228_all_sed_supt_BULKURINE_urineONLY_deg.csv False
