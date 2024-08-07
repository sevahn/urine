#!/bin/bash
#SBATCH --job-name=plasma_sediment.sh
#SBATCH --output=plasma_sediment.sh.out
#SBATCH --error=plasma_sediment.sh.err
#SBATCH --qos=normal
#SBATCH --ntasks=1
#SBATCH --partition=quake,owners,normal
#SBATCH --cpus-per-task=1
#SBATCH --mem=10000
#SBATCH --time=8:00:00



source activate snakemake

python3 bootstrap_ci.py plasma_sediment 20231228_plasma_sediment_all_bioivt_urine_deg.csv False
