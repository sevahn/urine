import numpy as np
import pandas as pd
import scanpy as sc

import os
import random

import argparse
from datetime import date

parser = argparse.ArgumentParser(description = 'bootstrap CI for two groups')
parser.add_argument('group1_group2', type = str, help = 'format group1_group2 (ie "plasma_supt")')
parser.add_argument('deg_file', type = str, help = 'DEG file from voom')
parser.add_argument('is_normal', type = str, help = 'was the urine sample normal?')

args = parser.parse_args()
os.system('echo ' + args.group1_group2 + ' ' + args.deg_file)

countsf = "20231228_filtered_sed_supt_urineONLY.csv"
tmmf = "20231228_urineOnly_TMM.csv"

#comparison = args.group1_group2.split("_")
#group1, group2 = comparison[0], comparison[1] 

# for saving the bootstrapped DEGs
day_today = date.today()
fend = day_today.strftime("%Y%m%d") + ".csv"

# proceed with the bootstrapping...
random.seed(229)

cts_data = pd.read_csv(countsf,
                      sep = ",", index_col = (0, 1))
tmm_vals = pd.read_csv(tmmf,
                      sep = ",", index_col = 0)

tmm_vals.index = [i[1:] if i[0] == 'X' else i for i in tmm_vals.index]

if args.is_normal == 'True':
    os.system('echo NORMAL!')
    sed_id = [i for i in cts_data.columns if 'sediment' in i and 'N' in i]
    supt_id = [i for i in cts_data.columns if 'supt' in i and 'N' in i]
else:
    os.system('echo NOT NORMAL')
    sed_id = [i for i in cts_data.columns if 'sediment' in i]
    supt_id = [i for i in cts_data.columns if 'supt' in i]

plasma_id = [i for i in cts_data.columns if 'sed' not in i and 'supt' not in i]


samp_dict = {}
sed_id = [i for i in cts_data.columns if 'sediment' in i]
supt_id = [i for i in cts_data.columns if 'supt' in i]

samp_dict['sediment'] = sed_id
samp_dict['sedimentNormal'] = [i for i in sed_id if 'N' in i]
samp_dict['sedimentStone'] = [i for i in sed_id if 'N' not in i ]
samp_dict['supt'] = supt_id
samp_dict['suptStone'] = [i for i in supt_id if 'N' not in i]
samp_dict['suptNormal'] = [i for i in supt_id if 'N' in i]
samp_dict['plasma'] = plasma_id
samp_dict['urine'] = sed_id + supt_id

# get the two groups for the lfc comparison 
#group1_id = samp_dict[group1]
#group2_id = samp_dict[group2]

sed_leukneg_id = ["1747_sediment_repool", "1761_sediment", "N1_sediment", "N2_sediment", "1755_sediment_repool", "1756_sediment_repool", "1757_sediment_repool", "N4_sediment_repool", "N5_sediment_repool", "N6_sediment_repool"]
sed_leukpos_id = ["1741_sediment_repool", "1742_sediment", "1748_sediment_repool", "1754_sediment_repool", "221_sediment_repool", "1746_sediment_repool", "1749_sediment"]

supt_leukneg_id = ["1747_supt_purple_repool", "1757_supt_repool", "N1_supt_repool", "N2_supt_repool", "N3_supt_repool", "N5_supt_repool", "N6_supt_repool"]   
supt_leukpos_id = ["1742_supt_purple_repool", "1748_supt_purple_repool", "1749_supt_purple", "221_supt_purple_repool", "1741_supt_repool", "1746_supt_repool"]
print(len(sed_leukneg_id) + len(supt_leukneg_id) + len(sed_leukpos_id) + len(supt_leukpos_id))

if args.group1_group2 == "leukpos_sed_supt":
    group1_id = sed_leukpos_id 
    group2_id = supt_leukpos_id

elif args.group1_group2 == "leukneg_sed_supt":
    group1_id = sed_leukneg_id
    group2_id = supt_leukneg_id
 
elif args.group1_group2 == "all_sed_supt":
    group1_id = sed_leukpos_id + sed_leukneg_id
    group2_id = supt_leukpos_id + supt_leukneg_id
 
os.system('echo group1_id = ' + str(group1_id))
os.system('echo group2_id = ' + str(group2_id))
os.system('echo nous sommes ici!')

cpm = cts_data.div(cts_data.sum(axis = 0), axis = 1) * 10 ** 6

prior_count = 1 
# mimic cpm function from edgeR https://www.rdocumentation.org/packages/edgeR/versions/3.14.0/topics/cpm
cpm_tmm = cpm.div(tmm_vals['norm.factors'])
cpm_tmm = cpm_tmm.swaplevel()
cpm_tmm = np.log2(cpm_tmm + prior_count)


#os.system('echo cpm.head() ', cpm_tmm.head())
qThresh = 0.05
normal_degs = pd.read_csv(args.deg_file, index_col = 0)
normal_degs = normal_degs[normal_degs['adj.P.Val'] <= qThresh]


num_bootstrap = 1000
ci = 0.95
alpha = 1 - ci

failed_genes = []
passing_genes = []
pass_mask = [] # use this to mask the deg up/down
ci_vals = []
true_lfcs = []

for g in normal_degs.index:
    cts = cpm_tmm.loc[g].T
    
    group1_cts = cts.loc[group1_id] 
    group2_cts = cts.loc[group2_id] 
    
    # get the true lfc
    true_lfc = np.median(group1_cts) - np.median(group2_cts)
    
    all_lfc = []
    for boot in range(num_bootstrap):
        boot_sed = group1_cts.sample(replace = True, frac = 1)
        boot_supt = group2_cts.sample(replace = True, frac = 1)
        
        lfc = boot_sed.median().values[0] - boot_supt.median().values[0]
        all_lfc += [lfc]
    
    qHat_lower = np.quantile(np.sort(all_lfc), (1 - alpha / 2))
    qHat_upper = np.quantile(np.sort(all_lfc), (alpha / 2))
    gene_ci = (2 * true_lfc - qHat_lower, 2 * true_lfc - qHat_upper)
    
    if true_lfc > 0:
        if gene_ci[0] < 0 or gene_ci[1] < 0:
           pass_mask += [False]
        else: pass_mask += [True]  
    
    elif true_lfc < 0:
        if gene_ci[0] > 0 or gene_ci[1] > 0:
            pass_mask += [False]
        else: pass_mask += [True]          
    if true_lfc == 0:
        pass_mask += [np.nan]
        os.system('echo ' + g + ' true LFC is weird') 
    
    ci_vals += [gene_ci]
    true_lfcs += [true_lfc]

resDF = pd.DataFrame()
resDF['genes'] = normal_degs.index.tolist()
resDF['95ci_lower'] = [i[0] for i in ci_vals] 
resDF['95ci_upper'] = [i[1] for i in ci_vals]
resDF['passed'] = pass_mask
resDF['trueLFC'] = true_lfcs

os.system('echo done computing!')
os.system('echo ' + args.group1_group2 + '_all_DEG_95CI_' + fend)

fname = args.deg_file[:-4] + "DEG95CI.csv" 
resDF.to_csv(fname, header = True)