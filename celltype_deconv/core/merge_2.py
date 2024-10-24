# merge together all the csv files, one for support vectors and the other for coefficients
# this script goes into each of the deconvolved sample directories
# and merges the coefficient files and the support vector files
# the best coefficient will then be chosen locally 

import pandas as pd
import numpy as np
import scanpy as sc

import os

# dataset nam
datasetName = f"NCI_cfRNA" # name of dataset
fend = "_20231019" # name of file ending...
basePath = os.getcwd() + "/"

#### do not edit below this line

# read in samples
sampNames = [i for i in os.listdir()]

supVecBase = "_coarsegrain__supportVectors.csv"
coefBase = "_deconvolutionCoefs.csv"
predBase = '-NUSVRsvr_gene_preds.csv'

allSuppVec = pd.DataFrame()
allCoefs = pd.DataFrame()
allPreds = pd.DataFrame()

# concatenate all the samples over all hyperparam combinations
for i, bioRep in enumerate(sampNames):
    if i % 100 == 0:
        os.system('echo through ' + str(i) + ' samples out of ' + str(len(sampNames)) + '.')
    # check if the directory exists and its contents are nonempty (directory is created before files are written) 
    dirpath = f"{basePath}{bioRep}/"
    print(dirpath)
    if os.path.isdir(dirpath) and len(os.listdir(dirpath)) > 0:
        print(os.listdir(dirpath))
        files = os.listdir(f"{dirpath}")
        print(files)
        print([i for i in files if supVecBase in i])
        supp_vec = [i for i in files if supVecBase in i][0]
        coefs = [i for i in files if coefBase in i][0]
        preds = [i for i in files if predBase in i][0]
        sampSV = pd.read_csv(f"{basePath}{bioRep}/{supp_vec}", sep = ",", index_col = 0)
        sampCoef = pd.read_csv(f"{basePath}{bioRep}/{coefs}", sep = ",", index_col = 0)
        sampPreds = pd.read_csv(f"{basePath}{bioRep}/{preds}", sep = ",", index_col = 0)

        allCoefs = pd.concat([sampCoef, allCoefs], ignore_index = False)
        allSuppVec = pd.concat([sampSV, allSuppVec], axis = "columns", ignore_index = False)
        allPreds = pd.concat([sampPreds, allPreds], axis = 'columns', ignore_index = False)

# pick the best hyperparam combination per sample
bestCoef = pd.DataFrame()

for samp in sampNames:
    allHyperThisSamp = [i for i in allCoefs.index if samp in i]
    hyperCombosThisSamp = allCoefs.loc[allHyperThisSamp]
    best = hyperCombosThisSamp[hyperCombosThisSamp['rmse'] == np.min(hyperCombosThisSamp['rmse'])]
    bestCoef = pd.concat([bestCoef, best])

performance = bestCoef.iloc[:,-2:].T
bestCoef = bestCoef.iloc[:,:-2]

# send negative coefs to zero
bestCoef[bestCoef < 0] = 0 # this is a samp x cells matrix
bestCoef = bestCoef.T

# normalize by total to get the fractions of cell type specific RNA
fracs = bestCoef.div(bestCoef.sum(axis = 0), axis = 1) 

fracs = pd.concat([fracs, performance], axis = "index")

# get the support vectors corresponding to the best coefficient pair
suppvecs = allSuppVec[bestCoef.columns.tolist()]
preds = allPreds[bestCoef.columns.tolist() + [i for i in allPreds.columns if 'round' in i]]

# strip the hyperparameter information so it's just the sample names
bestCoef.columns = [i.split("-NUSVR")[0] for i in bestCoef.columns.tolist()]
suppvecs.columns = [i.split("-NUSVR")[0] for i in bestCoef.columns.tolist()]

# write out the support vectors and the fractions
suppvecs.to_csv(f'{datasetName}_support_vectors_{fend}.csv', sep = ",", header = True, index = False)
fracs.to_csv(f'{datasetName}_fractions_{fend}.csv', sep = ",", header = True, index = True)
preds.to_csv(f'{datasetName}_preds_{fend}.csv', sep = ",", header = True, index = True)
