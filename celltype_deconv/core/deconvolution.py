# the essentials
import numpy as np
import pandas as pd

# for computing RMSE 
from sklearn.metrics import mean_squared_error

# for pearson correlation  
from scipy import stats

# for file writing 
import os

from decon_methods import NNLS, nuSVR, QP

def deconvolve(scaledMixtures, scaledSigMat, deconMethod, biologRepName, savePredictions, outPath):
    """
    Perform deconvolution and return performance metrics (pearson R and RMSE)

    Parameters
    ----------
    
    scaledMixtures: pd.DataFrame
        Mixture scaled to zero mean and unit variance

    scaledSigMat : pd.DataFrame
        Basis matrix scaled to zero mean and unit variance

    deconMethod: str (either, 'qp', 'nnls', or 'nusvr')
        Deconvolution method to perform on samples 
        
    biologRepName: str
        Name of biological replicate that will be deconvolved for labeling in the output support vector/coefficient csv files
    
    savePredictions: boolean
        Set to true to save gene predictions to a CSV file. 
    
    outPath : str 
        Path to folder where output files are saved. 

    Returns
    -------
    nothing: outputs relevant files along the way. 
             (e.g. support vectors per nu/C combo + corresponding learned coefs)
    """
    ######### DECONVOLUTION TIME!!! ##########
    srrName = biologRepName
    savePath = os.path.join(outPath, srrName)
    if not os.path.exists(savePath):
        os.mkdir(savePath)
    savePath += "/"

    coefFile = srrName + "_deconvolutionCoefs.csv" 
    srrName += "-" + deconMethod.upper()
 
    ### INITIALIZE FILE THAT WILL BE SAVED ##
    predDF = pd.DataFrame(index = scaledSigMat.index.tolist())
    predDF['ground_truth_' + srrName] = scaledMixtures # since mixtures were prefiltered before passing to deconvolution fns 
    if os.path.exists(savePath + coefFile):
        aggCoefDF = pd.read_csv(savePath + coefFile, sep = ",", index_col = 0)
    else:
        os.system("creating coefficient file")
        aggCoefDF = pd.DataFrame(columns = scaledSigMat.columns.tolist() + ["r", "rmse"]) 
   
    if deconMethod.upper() == "QP":
        os.system("echo Performing Quadratic Programming")
        coefs = QP(scaledSigMat, scaledMixtures)
        coefDF = pd.DataFrame(data = coefs, index = [srrName], columns = scaledSigMat.columns)        
        groundTruth = np.asarray(scaledMixtures.values, dtype = "float").squeeze()
       
        # non-negativity & sum-to-one is implicit in QP function, use coefs as-is 
        preds = np.asarray(scaledSigMat.values.dot(coefs.T), dtype = "float").squeeze()
        pearsonR_QP, pearsonP_QP = stats.pearsonr(groundTruth, preds)
        rmse_QP = np.sqrt(mean_squared_error(preds, groundTruth))

        # concat to coefficient dataframe
        coefDF["r"] = [pearsonR_QP]
        coefDF["rmse"] = [rmse_QP]
        aggCoefDF = pd.concat([aggCoefDF, coefDF])
        predDF[srrName + '-QP_predictions'] = preds

    if deconMethod.upper() == "NNLS":
        os.system("echo performing NNLS")
        coefs, resid = NNLS(scaledSigMat, scaledMixtures, sumToOne = False)
        dfSampName = srrName
        coefDF = pd.DataFrame(data = coefs, index = [dfSampName], columns = scaledSigMat.columns)
        
        groundTruth = np.asarray(scaledMixtures.values, dtype = "float").squeeze()
        
        normCoefs = coefDF # make a copy of coefs
        normCoefs[normCoefs < 0] = 0 # set < 0 -> 0, implicit in NNLS
        
        # enforce sum to one
        normCoefs = normCoefs / normCoefs.sum(axis = "columns").values[0]
        
        # get dot product on predictions
        preds_norm = np.asarray(scaledSigMat.values.dot(normCoefs.T), dtype = "float").squeeze()
        pearsonR_NNLS, pearsonP_NNLS = stats.pearsonr(groundTruth, preds_norm)
        coefDF["rmse"] = [np.sqrt(mean_squared_error(preds_norm, groundTruth))]
        coefDF['r'] = [pearsonR_NNLS] 
         
        aggCoefDF = pd.concat([aggCoefDF, coefDF])
        predDF[srrName + '-NNLS_predictions'] = preds_norm

    if deconMethod.upper() == "NUSVR":
        os.system("echo Performing nu-SVR")
        nuVals = [0.05, 0.1, 0.15, 0.25, 0.5, 0.6, 0.75, 0.8, 0.85]
        cVals = [0.1, 0.25, 0.5, 0.75, 1, 10]
        
        coefDF, svDictDF, allCorrDF = pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
 
        # do a grid search for C/nu combinations
        scaledMixtureList = np.ndarray.flatten(np.asarray(scaledMixtures.values.tolist())).tolist()
        for cVal in cVals:
            for nuVal in nuVals:
                os.system("echo nu = " +  str(nuVal))
                os.system("echo C = " +  str(cVal))
                coefs, svDict, fileprefix, clfPred, intercept = nuSVR(scaledSigMat, scaledMixtures, nuVal, cVal)
                nuSampName = [srrName + "-nu=" + str(nuVal) + "-C=" + str(cVal)] 
                
                # iteratively write to the coefficient file
                coefUpdateDF  = pd.DataFrame(data = coefs, columns = scaledSigMat.columns, 
                                index = nuSampName)
                coefDF = pd.concat([coefDF, coefUpdateDF])

                # account for non-uniform length of values in svDict
                # update the key in the dictionary to account for the nu value and save it
                
                # update the keys of the dictionary to the full sample name
                print('srrName = ', srrName)
                print(srrName)
                print(svDict.keys())
                svDict[nuSampName[0]] = svDict.pop(biologRepName) 
                svUpdate = pd.DataFrame(dict([ (k, pd.Series(v)) for k, v in svDict.items()]))
                svDictDF = pd.concat([svDictDF, svUpdate], axis = "columns")

                # get the correlations
                corrDF = pd.DataFrame(index = nuSampName, columns = ["r", "rmse"])
                
                # get the normalized coefs and compute RMSE on predictions
                normCoef = coefs # make copy of coefficients
                normCoef[normCoef < 0] = 0 # set the coefs < 0 -> 0
                normCoef = normCoef / normCoef.sum() # normalize to get the fractional contribs # i bet this is a row so you had to make it a vector
                predExpr_normCoefs  = np.asarray(np.asmatrix(scaledSigMat).dot(normCoef.T)).squeeze() # use the whole basis matrix
                
                rmse_nuSVR = np.sqrt(mean_squared_error(scaledMixtureList, predExpr_normCoefs))
                pearsonR_nuSVR, pearsonP_nuSVR = stats.pearsonr(scaledMixtureList, predExpr_normCoefs)
                os.system('echo rmse = ' +  str(rmse_nuSVR))
                corrDF['rmse'] = [rmse_nuSVR]
                corrDF['r'] = [pearsonR_nuSVR]

                clfPred = clfPred.tolist() # so this is the classifier prediction
   				
                predDF[nuSampName[0]] = predExpr_normCoefs
                allCorrDF = pd.concat([corrDF, allCorrDF])
        
        coefDF = pd.concat([coefDF, allCorrDF], axis = "columns") # append along columns

        aggCoefDF = pd.concat([aggCoefDF, coefDF])
        # put everything in its own folder
        svDictDF.to_csv(savePath + fileprefix + "_supportVectors.csv", sep = ",",
                 header = True)   
    if savePredictions: 
        predDF.to_csv(savePath +  srrName + 'svr_gene_preds.csv', sep = ",", header = True, index = True)
    aggCoefDF.to_csv(savePath + coefFile, sep = ",", header = True, index = True) 

