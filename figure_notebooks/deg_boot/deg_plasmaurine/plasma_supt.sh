#!/bin/bash
#SBATCH --job-name=plasma_supt.sh
#SBATCH --output=plasma_supt.sh.out
#SBATCH --error=plasma_supt.sh.err
#SBATCH --qos=normal
#SBATCH --ntasks=1
#SBATCH --partition=quake,owners,normal
#SBATCH --cpus-per-task=1
#SBATCH --mem=10000
#SBATCH --time=8:00:00



source activate snakemake

python3 bootstrap_ci.py plasma_supt 20231228_plasma_supt_all_bioivt_urine_deg.csv False
