# essentials
import pandas as pd

# scale mixture/matrix to zero mean/unit var
from sklearn import preprocessing


def processMixture(rawCountsDF, cpmThresh, cpmNormalize):
    """
    Read in mixtures, eliminate genes that are all zero counts in sample, 
    normalize to CPM.

    Parameters
    ----------
    rawCountsDF : pd.DataFrame() 
       df of the raw counts 
        
    cpmThresh: float
        cpmVals less than this are sent to zero 
    
    cpmNormalize: bool
        perform cpmNormalization if True, else False
    
    Returns
    -------
    thuyCPM: pd.DataFrame (genes x mixtures)
        dataframe with CPM counts as values in matrix
    """
    
    if cpmNormalize:
        thuyCPM = rawCountsDF.div(rawCountsDF.sum(axis = 0), axis = 1) * (10 ** 6)
    else: thuyCPM = rawCountsDF
            
    # set everything less than the CPM thresh to zero
    thuyCPM[thuyCPM < cpmThresh] = 0
    thuyCPM = thuyCPM[(thuyCPM.T != 0).any()] # TODO: replace with drop for improved computation
        
    return(thuyCPM)


def scale(df):
    """
    Scale data in DF to zero mean and unit variance to improve runtime performance

    Parameters
    ---------
    sigMat: pd.DataFrame 
        Basis matrix that will be used to deconvolve samples

    mixture: pd.DataFrame 
        Mixture that will be deconvolved
 
    Returns
    --------
    scaledDF: pd.DataFrame
        Data scaled to zero mean and unit variance    
    """

    scaledDF = preprocessing.scale(df.values)
    scaledDF = pd.DataFrame(data = scaledDF, index = df.index, columns = df.columns) 
    return scaledDF 
