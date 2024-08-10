#!/bin/bash
#SBATCH --job-name=deconsample1cfRNA_deconv
#SBATCH --output=deconsample1cfRNA_deconv.out
#SBATCH --error=deconsample1cfRNA_deconv.err
#SBATCH --qos=normal
#SBATCH --ntasks=1
#SBATCH --partition=quake,owners,normal
#SBATCH --cpus-per-task=1
#SBATCH --mem=20000
#SBATCH --time=24:00:00



source activate deconv

python3 /home/users/hagop/cellfracker_pub2/core/deconv_wrapper.py --basis-matrix-file /home/users/hagop/cellfracker_pub2/tutorial/data/tsp_v1_basisMatrix.txt --mixture-path /home/users/hagop/cellfracker_pub2/tutorial/data/samples_test.csv --biolog-rep-name sample1 --out-path /home/users/hagop/cellfracker_pub2/tutorial/outputs/ --cpm-threshold 0 --deconv-method nusvr --save-predictions --no-cpm-normalization