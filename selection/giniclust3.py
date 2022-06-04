# -*- coding: utf-8 -*-
# @Time : 2022/6/3 14:42
# @Author : Tory Deng
# @File : giniclust3.py
# @Software: PyCharm
import math

import numpy as np
from statsmodels.nonparametric.smoothers_lowess import lowess
from scipy import stats
from scipy.interpolate import interp1d
from sklearn import preprocessing


def giniIndexCalculation(array):
    plusArray = array + 0.0000001
    plusArray = np.sort(plusArray)
    index = np.arange(1, plusArray.shape[0] + 1)
    n = plusArray.shape[0]
    indexValue = (np.sum((2 * index - n - 1) * plusArray)) / (n * np.sum(plusArray))
    return indexValue


def giniIndex(filteredArray):
    funcGeneGini = []
    funcGeneMax = []
    for i in range(0, len(filteredArray)):
        exp = filteredArray[i]
        maxExp = np.max(exp)
        giniIndexValue = giniIndexCalculation(exp)
        funcGeneGini.append(giniIndexValue)
        funcGeneMax.append(float(maxExp))
    return funcGeneGini, funcGeneMax


def loessRegression(funcGeneGini, funcLogGeneMax, funcGene):
    """
    Calculate per gene p-vals. Copied and polished from the source code of GiniClust3.

    Parameters
    ----------
    funcGeneGini
      Per gene gini values
    funcLogGeneMax
      log-transformed values
    funcGene
      Gene names

    Returns
    -------
    funcSigGiniGeneGini, funcSigGiniGenePvalue
      per gene gini and p-vals in dicts
    """
    dictGini = {}
    for i in range(len(funcGene)):
        dictGini[funcGene[i]] = funcGeneGini[i]
    fit = lowess(funcGeneGini, funcLogGeneMax, frac=0.9)
    f = interp1d(list(zip(*fit))[0], list(zip(*fit))[1], bounds_error=False)
    giniFit = f(funcLogGeneMax)
    residue = funcGeneGini - giniFit
    posRes = [residue[i] for i in range(len(residue))]
    residueSort = sorted(posRes)
    quarter = int(len(residueSort) * 3 / 4)
    cutoff = residueSort[quarter]

    quantileGini, quantileLogMax, quantileGene, outlierGini, outlierLogMax, outlierGene = [], [], [], [], [], []
    for i in range(0, len(residue)):
        if residue[i] <= cutoff:
            quantileGene.append(funcGene[i])
            quantileGini.append(funcGeneGini[i])
            quantileLogMax.append(funcLogGeneMax[i])
        else:
            outlierGini.append(funcGeneGini[i])
            outlierGene.append(funcGene[i])
            outlierLogMax.append(funcLogGeneMax[i])

    fit2 = lowess(np.array(quantileGini), np.array(quantileLogMax), frac=0.9)
    reFit = interp1d(list(zip(*fit2))[0], list(zip(*fit2))[1], bounds_error=False)
    quantileGiniReFitPredict = reFit(np.array(quantileLogMax))
    dictFunction = {}
    for i in range(len(quantileLogMax)):
        dictFunction[quantileLogMax[i]] = quantileGiniReFitPredict[i]
    uniqQuantileLogMax = list(dictFunction.keys())
    sortQuantileLogMax = sorted(uniqQuantileLogMax, reverse=True)
    k = (dictFunction[sortQuantileLogMax[1]] - dictFunction[sortQuantileLogMax[2]]) / (
            sortQuantileLogMax[1] - sortQuantileLogMax[2])
    b = dictFunction[sortQuantileLogMax[1]] - k * sortQuantileLogMax[1]

    outlierGiniReFitPredict = reFit(outlierLogMax)
    # remove value NaN
    for i in range(0, len(outlierGiniReFitPredict)):
        if math.isnan(outlierGiniReFitPredict[i]):
            outlierGiniReFitPredict[i] = outlierLogMax[i] * k + b

    newGene = quantileGene + outlierGene
    newGini = np.array(list(quantileGini) + list(outlierGini))
    newFitGini = np.array(list(quantileGiniReFitPredict) + list(outlierGiniReFitPredict))
    newResidualGini = newGini - newFitGini

    pvalue = stats.norm.cdf(-abs(preprocessing.scale(newResidualGini)))
    funcSigGiniGeneGini, funcSigGiniGenePvalue = {}, {}
    for i in range(0, len(pvalue)):
        funcSigGiniGeneGini[newGene[i]] = dictGini[newGene[i]]
        funcSigGiniGenePvalue[newGene[i]] = pvalue[i]
    return funcSigGiniGeneGini, funcSigGiniGenePvalue
