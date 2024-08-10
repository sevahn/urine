# each time you run this script
# change the seed if you want different values
# check num boot iters
# check nuVals and cVals

# sevahn
# deconv code v3 
# introduce bootstrapping to get a CI on a given cell type fraction
# TO DO - introduce the grid search
# crash on line 412 if the row space of the basis matrix is zero
# this script is from /home/groups/quake/sevahn/pee_pilot/deconvolve_peev3_wpreds_20220627/deconvolve_v2_20220622.py
# renamed to deconv_wpreds.py

# objective: get gene predictions to see what's going on with the urine

# the essentials
import numpy as np
import pandas as pd

# helper functios
from preprocess import processMixture, scale 
from deconvolution import deconvolve

# for file writing and timing deconvolution
import os

# for parsing commandline args
import argparse

######### no need to edit below this line #########

def deconvolutionWrapper(countDF, cpmThresh, basisMatrix, biologRepName, deconMethod, doCpmNorm, savePredictions, outPath): 
    """
    countDF: pd DataFrame  
       df of the raw counts 

    cpmThresh: double
        CPM values less than this are set to zero
    
    basisMatrix: np.ndarray
       The TSP basis matrix. 
         
    biologRepName: str
        Name of biological replicate that will be deconvolved for labeling in the output support vector/coefficient csv files

    deconMethod: str (either, 'qp', 'nnls', or 'nusvr')
        Deconvolution method to perform on samples 

    doCpmNorm : boolean
        Set to true to perform CPM normalization, false otherwise. 
    
    savePredictions: boolean
        Set to true to save gene predictions to a CSV file. 
    
    outPath : str 
        Path to folder where output files are saved. 

    Wrapper to the deconvolve function that performs the appropriate preprocessing/scaling prior to the deconvolve call. 
    """

    thuyCPM = processMixture(countDF, cpmThresh, doCpmNorm)

    # deconvolution can only be performed on the gene intersection of what's in the mixture and the
    # row space of the basis matrix (e.g. genes) -- reindex both to the gene intersection only 
    intersection = np.intersect1d(thuyCPM.index, basisMatrix.index)
    relThuy = thuyCPM.loc[intersection, :]
    relBM = basisMatrix.loc[intersection, :]
    
    # drop gene duplicates if present 
    val, counts = np.unique(relThuy.index, return_counts = True)
    dups = val[counts > 1].tolist()
    relThuy = relThuy.drop(dups)
    relBM = relBM.drop(dups)

    # perform preprocessing about the index to improve runtime performance
    sigMat = relBM
    mixture = relThuy
    
    # scale along the columns (e.g. a given sample or cell type)
    scaledSigMat = scale(sigMat)
    # also scale the mixture separately to zero mean and unit variance
    scaledMixture = scale(mixture)

    # it's deconvolution time!
    deconvolve(scaledMixture, scaledSigMat, deconMethod, biologRepName, savePredictions, outPath)

def main(biologRepName, mixturePath, cpmThresh, deconMethod, basisMatrixPath, doCpmNorm, savePredictions, outPath):
    """
    Main deconvolution function that calls necessary helper functions to perform deconvolution
    Parameters
    ---------
    biologRepName : str 
        Name of sample that will be deconvolved

    mixturePath: str
        Filepath to counts matrix of samples for deconvolution 

    cpmThresh: int TODO: should this be float/double?
        CPM values less than this are set to zero

    deconMethod: str (either: 'nusvr', 'qp', 'nnls'; case insensitive)
        Deconvolution method to apply on samples
    
    doCpmNorm : boolean
        Set to true to perform CPM normalization. Otherwise, CPM normalization is not performed.

    Returns
    --------
    Nothing, the functions that it calls return/save files as needed
    """
    # read in the basis matrix 
    basisMatrix = pd.read_csv(basisMatrixPath, sep = "\t", index_col = 0) 

    # read in the bulk RNA counts table for deconvolution 
    # drop the zero genes
    countDF = pd.read_csv(mixturePath, sep = ",", index_col = 0)
    countDF = countDF[biologRepName].to_frame() 
    countDF = countDF.loc[(countDF != 0).any(axis = 1)]
    
    mixture = countDF 
    samp_name = biologRepName
    deconvolutionWrapper(mixture, cpmThresh, basisMatrix, samp_name, deconMethod, doCpmNorm, savePredictions, outPath)

if __name__ == "__main__": 
    """
    Parses various arguments.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("--basis-matrix-file", type=str, help='Path to basis matrix file.')
    parser.add_argument("--cpm-threshold", type=float, default=0)
    # Justification for CPM normalization flags: https://stackoverflow.com/questions/15008758/parsing-boolean-values-with-argparse - HJC
    parser.add_argument('--do-cpm-normalization', dest='do_cpm_norm', action='store_true')
    parser.add_argument('--no-cpm-normalization', dest='do_cpm_norm', action='store_false')
    parser.set_defaults(do_cpm_norm=None)
    # https://stackoverflow.com/questions/2020598/in-python-how-should-i-test-if-a-variable-is-none-true-or-false
    parser.add_argument('--save-predictions', action='store_true')
    parser.add_argument('--out-path', type=str)
    parser.add_argument('--mixture-path', type=str) 
    parser.add_argument('--deconv-method', type=str)
    parser.add_argument('--biolog-rep-name', type=str) 

    args = parser.parse_args()
    if (args.do_cpm_norm is None):
        exit("Must specify whether or not to apply CPM normalization using the --do-cpm-normalization or --no-cpm-normalization flags, respectively.")
        
    main(args.biolog_rep_name, args.mixture_path, args.cpm_threshold, args.deconv_method, args.basis_matrix_file, args.do_cpm_norm, args.save_predictions, args.out_path)
