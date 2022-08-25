import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from importlib import reload
import sys
import warnings
from scipy.stats import norm
from scipy.stats import chi2

def nll_data(betas,ses,alpha, uniqueness = True):
    #assumes independence, but we can generalize
    #print(alpha)
    #alpha[0] = 1
    alpha = np.abs(np.array(alpha))
    if uniqueness:
        alpha[0] = 1
    alpha[alpha < .01] = .01
    alpha[alpha > 10] = 10
    #compute meta_means
    #print(alpha)
    betas = betas @ np.diag(alpha)
    ses = ses @ np.diag(alpha)
    #print(betas.head())
    weights = np.power(ses, -2)
    weights = np.divide(weights, weights.sum(axis = 1).values.reshape(-1,1))
    meta_means = np.multiply(betas, weights).sum(axis = 1)
    #print(meta_means.head())
    quadratic = np.power(np.divide(np.array(betas) - np.array(meta_means).reshape(-1,1), np.array(ses)),2).sum().sum()
    return .5*quadratic

def nll_search2D(betas, ses, alphas = np.linspace(.1,5,500)):
    df_results = pd.DataFrame(columns = ['Alpha','NLL'])
    for i, alpha in enumerate(alphas):
        nll = nll_data(betas, ses, [1, alpha])
        df_results.loc[i] = [alpha, nll]
    return df_results

def inv_var_meta(betas, ses):
    inv_var = np.power(ses, -2)
    num = np.multiply(betas, inv_var).sum(axis = 1)
    denominator = np.sqrt(inv_var.sum(axis = 1))
    meta = np.divide(num, denominator)
    return meta


def compute_meta_scores(betas, ses, alpha, norm_alpha = True):
    if norm_alpha:
        alpha[0] = 1
    alpha = np.abs(alpha)
    betas = betas @ np.diag(alpha)
    ses = ses @ np.diag(alpha)
    meta = inv_var_meta(betas, ses)
    return meta

def compute_meta_scores_mean(betas, ses, alpha, norm_alpha = True):
    meta = compute_meta_scores(betas, ses, alpha, norm_alpha = norm_alpha)
    return np.mean(np.abs(meta))
