#!/bin/bash
#SBATCH --job-name=normal_plasma_sediment.sh
#SBATCH --output=normal_plasma_sediment.sh.out
#SBATCH --error=normal_plasma_sediment.sh.err
#SBATCH --qos=normal
#SBATCH --ntasks=1
#SBATCH --partition=quake,owners,normal
#SBATCH --cpus-per-task=1
#SBATCH --mem=10000
#SBATCH --time=8:00:00



source activate snakemake

python3 bootstrap_ci.py plasma_sediment 20231228_plasma_sediment_normal_bioivt_urine_deg.csv True
