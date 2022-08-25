import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from importlib import reload
import sys
import warnings
from tqdm import tqdm_notebook as tqdm
#from scipy import minimize
from scipy.stats import norm
#import sympy as sym
from scipy.stats import chi2




def compute_meta_score(betas, ses, alpha, norm_alpha = True):
    if norm_alpha:
        alpha[0] = 1
    alpha = np.abs(alpha)
    betas = betas @ np.diag(alpha)
    ses = ses @ np.diag(alpha)
    meta = inv_var_meta(betas, ses)
    return np.abs(meta).mean()

def inv_var_meta(betas, ses):
    #alpha[0] = 1
    #betas = betas @ np.diagonal(alpha)
    #ses = ses @ np.diagonal(alpha)
    inv_var = np.power(ses, -2)
    #return a vector
    #df_meta['NUM'] = df_meta['INV_VAR_' + diluted_key]*(df_meta['BETA_' + diluted_key]) + (df_meta['INV_VAR_' + base_key])*(df_meta['BETA_' + base_key])
    num = np.multiply(betas, inv_var).sum(axis = 1)
    denominator = np.sqrt(inv_var.sum(axis = 1))
    meta = np.divide(num, denominator)
    return meta

#write functions for other meta analysis techniques
def bolormaa_meta(betas, ses):
    #requires inputs be dfs
    dof = betas.shape[1]
    z_scores = np.divide(betas, ses)
    corr_matrix = z_scores.corr()
    inv_corr = np.linalg.inv(corr_matrix)
    chi_sq_stats = np.multiply(z_scores.values.transpose(), inv_corr @ z_scores.values.transpose())
    chi_sq_stats = chi_sq_stats.sum(axis = 0)
    return 1 - chi2.cdf(chi_sq_stats, dof)

def compute_eff_sample_size(cases, controls):
    inv_cases = 1/cases
    inv_controls = 1/controls
    return 2/(inv_cases + inv_controls)

#need to test these functions line by line to make sure they work as intended
def cpassoc_meta(betas, ses, eff_sample_sizes):
    dim_beta = betas.shape[1]
    z_scores = np.divide(betas, ses)
    weights = np.sqrt(np.diag(eff_sample_sizes))
    corr_matrix = z_scores.corr()
    inv_corr = np.linalg.inv(corr_matrix @ weights)
    chi_sq_denominator = np.ones(dim_beta) @ np.linalg.inv(weights @ corr_matrix @ weights) @ np.ones(dim_beta)
    chi_sq_num = (inv_corr @ z_scores.values.transpose()).sum(axis = 0)
    chi_sq_num = np.power(chi_sq_num, 2)
    chi_sq_stats = chi_sq_num/chi_sq_denominator
    return 1 - chi2.cdf(chi_sq_stats, 1)

def cpassoc_meta_filter(betas, ses, eff_sample_sizes, threshold = 1.96):
    dim_beta = betas.shape[1]
    z_scores = np.divide(betas, ses)
    weights = np.sqrt(np.diag(eff_sample_sizes))
    #compute corr matrix
    z_scores_filtered = z_scores.loc[np.abs(z_scores.max(axis = 1)) < threshold]
    corr_matrix = z_scores_filtered.corr()
    inv_corr = np.linalg.inv(corr_matrix @ weights)
    chi_sq_denominator = np.ones(dim_beta) @ np.linalg.inv(weights @ corr_matrix @ weights) @ np.ones(dim_beta)
    chi_sq_num = (inv_corr @ z_scores.values.transpose()).sum(axis = 0)
    chi_sq_num = np.power(chi_sq_num, 2)
    chi_sq_stats = chi_sq_num/chi_sq_denominator
    return 1 - chi2.cdf(chi_sq_stats, 1)


def weighted_z_meta(p_vals, sign_beta, eff_sample_sizes):
    #sign beta can be z-scores or betas
    signs = np.sign(sign_beta)
    z_scores = norm.ppf(1 - (p_vals/2))
    z_scores = np.multiply(z_scores, signs)
    weights = np.sqrt(eff_sample_sizes)
    denominator = np.sqrt(np.power(weights, 2).sum())
    numerator = np.multiply(z_scores, weights).sum(axis = 1)
    meta_z_scores = numerator/denominator
    return 2*norm.cdf(-np.abs(meta_z_scores))
