"""
Generate individual slurm jobs for each sample to be bootstrapped! 
"""


import pandas as pd
import numpy as np
import os

########## ENTER ARGUMENTS HERE ###############

# path to your CPM-normalized counts
# IMPORTANT: first row should be gene names, column names should be samples

env = "snakemake" # name of the conda environment 

partition = "quake,owners,normal" # partition on supercomputer for jobs to be run 

### other params to tweak (need be) 
memory = "10000" # job memory
time = "8:00:00" # job run time (for jobs with large nu values they take longer) 
scriptName = "bootstrap_ci.py" # script name is deconvolve.py, change if different for you

########## DO NOT EDIT PAST THIS LINE ###########
 
def write_sh_file(fname_sh, scriptName, group1_group2, deg_file, is_normal):
    shname = fname_sh + ".sh"
    file = open(shname, "w+")
    file.write("#!/bin/bash")
    file.write("\n")
    file.write("#SBATCH --job-name=" + shname + "\n")
    file.write("#SBATCH --output=" + shname + ".out" + "\n")
    file.write("#SBATCH --error=" + shname + ".err" + "\n")
    file.write("#SBATCH --qos=normal\n")
    file.write("#SBATCH --ntasks=1\n")
    file.write("#SBATCH --partition=" + partition + "\n")
    file.write("#SBATCH --cpus-per-task=1\n")
    file.write("#SBATCH --mem=" + str(memory) + "\n")
    file.write("#SBATCH --time=" + time + "\n")
    file.write("\n\n\n")
    file.write("source activate " + env + "\n")
    
    file.write("\n")
    print(type(is_normal), 'is_normal', is_normal)
    # line to call nu-SVR deconvolution 
    file.write('python3 ' + scriptName + ' ' + group1_group2 + ' ' + deg_file + ' ' + str(is_normal))
    file.close()

for f in os.listdir():
    if "deg.csv" in f:
        deg_comparison = f.split("_")
        comparison = deg_comparison[1] + "_" + deg_comparison[2]
        is_normal = False
        print(comparison)
        if 'normal' in deg_comparison:
           print(deg_comparison)
           print(f)
           job_name = 'normal_' + comparison
           is_normal = True
        else: job_name = comparison 
        write_sh_file(job_name, scriptName, comparison, f, is_normal) 
