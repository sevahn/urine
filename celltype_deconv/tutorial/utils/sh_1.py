"""
Generate individual slurm jobs for each sample to be deconvolved!
"""


import pandas as pd
import numpy as np
import os
########## ENTER ARGUMENTS BELOW ###############

# path to your CPM-normalized counts
# IMPORTANT: first row should be gene names, column names should be samples
sampPath = "/home/users/hagop/cellfracker_pub2/tutorial/data/samples_test.csv"
basis_mat_path = "/home/users/hagop/cellfracker_pub2/tutorial/data/tsp_v1_basisMatrix.txt"
savePredictions = True 
doCpmNorm = False
cpmThresh = 0 # minimum gene expression level in sample (can increase/decrease need be)
deconMethod = "nusvr" # specify to "nuSVR" or "nnls" or "QP" (whichever kind of deconvolution you'd like to perform)
out_path = "/home/users/hagop/cellfracker_pub2/tutorial/outputs/"
env = "deconv" # name of the conda environment 
sh_save_dir = "/home/users/hagop/cellfracker_pub2/tutorial/shfiles/" # folder to save generated .sh files to 
deconv_wrapper_path = "/home/users/hagop/cellfracker_pub2/core/deconv_wrapper.py"
partition = "quake,owners,normal" # partition on supercomputer for jobs to be run 

### other params to tweak (need be) 
memory = "20000" # job memory
time = "24:00:00" # job run time (for jobs with large nu values they take longer) 

######## ENTER ARGUMENTS ABOVE ##################



########## DO NOT EDIT PAST THIS LINE ###########
importName = "deconv"
pre = "decon"
hyper = "cfRNA_deconv"

# read in list of samples, where the columns are the sample names
samps = list(pd.read_csv(sampPath, index_col = 0).columns)

biologRepSRR = {}
for i in samps:
    biologRepSRR[i] = [i]

for bioRep in biologRepSRR: 
    srr = biologRepSRR[bioRep]
    print(srr)
    shname = bioRep + hyper + ".sh"
    shname = os.path.join(sh_save_dir, shname)
    file = open(shname, "w+")
    file.write("#!/bin/bash")
    file.write("\n")
    file.write("#SBATCH --job-name=" + pre + bioRep + hyper + "\n")
    file.write("#SBATCH --output=" + pre + bioRep + hyper + ".out" + "\n")
    file.write("#SBATCH --error=" + pre + bioRep + hyper + ".err" + "\n")
    file.write("#SBATCH --qos=normal\n")
    file.write("#SBATCH --ntasks=1\n")
    file.write("#SBATCH --partition=" + partition + "\n")
    file.write("#SBATCH --cpus-per-task=1\n")
    file.write("#SBATCH --mem=" + str(memory) + "\n")
    file.write("#SBATCH --time=" + time + "\n")
    file.write("\n\n\n")
    file.write("source activate " + env + "\n")
    
    file.write("\n")

    # line to call nu-SVR deconvolution 
    python_call = []
    python_call.append("python3")
    python_call.append(f"{deconv_wrapper_path}")
    python_call.append(f"--basis-matrix-file {basis_mat_path}")
    python_call.append(f"--mixture-path {sampPath}")
    python_call.append(f"--biolog-rep-name {bioRep}")
    

    if not(os.path.exists(out_path)):
        os.mkdir(out_path)
    python_call.append(f"--out-path {out_path}")
    python_call.append(f"--cpm-threshold {cpmThresh}")
    python_call.append(f"--deconv-method {deconMethod}")

    if savePredictions: 
        python_call.append("--save-predictions")
    
    if doCpmNorm:
        python_call.append("--do-cpm-normalization")
    else: 
        python_call.append("--no-cpm-normalization")
    
    file.write(" ".join(python_call))
    file.close()
