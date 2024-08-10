from sklearn.svm import NuSVR # for nu-SVR
from scipy.optimize import nnls # for NNLS

# deconvolution functions
from cvxopt import matrix # for QP 
from cvxopt import solvers # for QP

import os 
import time 
import numpy as np


def nuSVR(sigMat, mixtures, nuVal, cVal):
    """
    Perform nu-SVR on mixtures given reference signature matrix and mixtures for deconvolution.
    Code inspired from AutoGeneS.

    Parameters
    ----------
    sigMat : pd.DataFrame (genes x cellTypes)
        Reference signature matrix for deconvolution (genes x cell types) in CPM space.
        Rows are genes and cols are the different cell/tissue types
    mixtures : pd.DataFrame (genes x mixtures)
        cfRNA mixtures for deconvolution in CPM space.
        Rows are genes and cols are the different mixtures.
    nuVal: float
        nu-SVR hyperparameter, is the upper bound on margin errors and lower bound on fraction of support vectors 
    cVal: float
        nu-SVR hyperparameter, is the regularization strength 


    Returns
    -------
    coefs : numpy array of (Mixtures x cellTypes)
        nu-SVR hyperplane coefficients >= 0 (else normalized to zero)
    svDict : dictionary
        keys = mixtures and values are lists of genes used to define the hyper place
    """
    clf = NuSVR(C = cVal, nu = nuVal, kernel = "linear")
    
    # coefs = (mixture x cell types matrix)
    # possessing the hyperplane coefficients of the sigMatrix columns 
    coefs = np.zeros((mixtures.shape[1], sigMat.shape[1]))
    svDict = {}
    fileprefix = ""
    for m in range(mixtures.shape[1]):
        fileprefix += mixtures.columns[m]
        thisMix = mixtures.iloc[:, m]
        t0 = time.time()  
        clf.fit(sigMat, thisMix)
        os.system("echo done fitting -- total fit time was " + str(time.time() - t0))
        coefficients = clf.coef_
        supportVecs = clf.support_ # mask of the indices of the support vectors
        
        # get the genes that were used in the deconvolution
        goodGenes = sigMat.index[supportVecs].to_list()
       
        sampName = mixtures.columns[m]
        os.system("echo sampName in nuSVR " + sampName)
        svDict[sampName] = goodGenes
        coefs[m] = coefficients[0]

        # make predictions on mixture
        preds = clf.predict(sigMat)
        intercept = clf.intercept_[0] * np.ones(sigMat.shape[0])
 
    fileprefix += "_coarsegrain_" 
    return(coefs, svDict, fileprefix, preds, intercept) 

def NNLS(sigMat, mixture, sumToOne = False):
    # call as is
    mixture = np.squeeze(mixture)
    coefs, resid = nnls(sigMat, mixture, maxiter = 10 ** 3)
    if sumToOne == True:
        # add sum to one constraint
        print("sum to one in NNLS not currently supported")
    
    # save the coefs/resid as a dataframe
    coefs = coefs.reshape((1, len(coefs)))
    return(coefs, resid) 

def QP(sigMat, mixture):
    """
    Perform quadratic programming. Referenced guide here for implementation:
    https://courses.csail.mit.edu/6.867/wiki/images/a/a7/Qp-cvxopt.pdf

    Parameters:
    -----------
    @param sigMat = pd.DataFrame of the scaled basis matrix
    @param mixture = pd.DataFrame of the scaled mixture

    Returns:
    -----------
    coefs, a numpy array of the learned coefficient weights
    """
    # set up matrix product
    P = matrix(sigMat.values.T.dot(sigMat.values))
    q = matrix(-1 * sigMat.values.T.dot(mixture.values))
    
    # constraint of coefficients summing to one
    dimMix = mixture.shape[0]
    numCells = sigMat.shape[1] 
    A = matrix(np.ones(numCells, dtype = "float").reshape((1, numCells)))
    b = matrix(np.ones(1, dtype = "float").reshape((1,1)))
    
    # constraint of coefficients each greater than or equal to zero (nonnegativity)
    # eg product with negative identity matrix is >= 0 
    G = matrix(-1.0 * np.identity(numCells))
    h = matrix(0.0, (numCells, 1))
    
    # solve the problem
    soln = solvers.qp(P, q, G, h, A, b) 
    coefs = soln["x"]
    coefs = np.asarray(coefs)
    coefs = coefs.reshape((1, len(coefs)))
    return(coefs)
