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

def nll_data_sample_overlap(betas,ses,alpha, covar_matrix_params, uniquenenss = True):
    #assumes independence, but we can generalize
    #print(alpha)
    if uniqueness:
        alpha[0] = 1
    ses = np.array(ses)
    #keep definition of alpha consistent with other function
    alpha = np.power(alpha, -1)
    betas = np.array(betas)
    #first compute sigma inverse
    #reshape into n x n matrix
    covar_matrix = np.array(covar_matrix_params).reshape((betas.shape[1],betas.shape[1]))
    np.fill_diagonal(covar_matrix, 1)
    #symmetrize covar_matrix
    covar_matrix = .5*(covar_matrix + covar_matrix.transpose())
    covar_matrix[covar_matrix > 1] = .95
    covar_matrix[covar_matrix < -1] = -.95
    #print(covar_matrix)
    alpha = np.abs(np.array(alpha))
    alpha[alpha < .01] = .01
    alpha[alpha > 100] = 100

    #compute covar matrix for each SNP
    covar_matrix_list = np.zeros((betas.shape[0], betas.shape[1], betas.shape[1]))
    meta_means = np.zeros(betas.shape[0])
    nll = 0
    covar_matrix_inv = np.linalg.inv(covar_matrix)
    inv_ses = np.power(ses, -1)
    covar_matrix_inv_full = np.repeat(covar_matrix_inv.reshape(1,betas.shape[1],betas.shape[1]),
                                      betas.shape[0], axis = 0)


    covar_matrix_total = np.einsum('ij,ijk,ik -> ijk',inv_ses, covar_matrix_inv_full,inv_ses)
    meta_means = np.einsum('j,ijk,ik -> i', alpha, covar_matrix_total, betas)
    meta_means_den = np.einsum('j,ijk,k -> i', alpha, covar_matrix_total, alpha)
    meta_means = np.divide(meta_means, meta_means_den)
    delta_beta = betas - np.einsum('j,i -> ij', alpha, meta_means) #alpha @ meta_means
    nll = .5*np.einsum('ij,ijk,ik -> i', delta_beta, covar_matrix_total, delta_beta).sum()
    nll += .5*betas.shape[0]*np.log(np.linalg.det(covar_matrix))

    return nll
