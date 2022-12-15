import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from importlib import reload
import sys
import warnings
#from tqdm import tqdm_notebook as tqdm
#from scipy import minimize
#from snp_sim_fix import percent_nulls
from scipy.stats import norm
#import sympy as sym
from scipy.stats import chi2
import meta_simple as meta
from scipy.optimize import minimize
from sklearn.metrics import average_precision_score as compute_aps
#decide range of values for p_effect
#decide range of values for sample_sizes
#decide range of values for dilution

def generate_betas(rng, effect_size = .02, n_trials = 1e5, p_effect = .1):
    betas = np.zeros(int(n_trials))
    number_of_effects = int(np.ceil(n_trials*p_effect))
    #scale is SD
    betas[0:number_of_effects] = rng.normal(loc = 0, scale = effect_size, size = number_of_effects)
    return betas

def generate_gwas_stats(rng,effect_size, n_trials, p_effect, dilution, sample_sizes):
    #dilution contains list of dilution values
    dilution = np.power(dilution, -1)
    n_trials = int(n_trials)
    n_studies = len(sample_sizes)
    beta_vars = ['BETA_' + str(i+1) for i in range(n_studies)]
    se_vars = ['SE' + str(i+1) for i in range(n_studies)]
    z_vars = ['Z' + str(i+1) for i in range(n_studies)]
    p_vars = ['P' + str(i+1) for i in range(n_studies)]
    gwas_cols = beta_vars + se_vars + ['TrueEffectSize']
    df_stats = pd.DataFrame(columns = gwas_cols)

    #print(df_stats.head())
    effects = generate_betas(rng, effect_size = effect_size, n_trials = n_trials, p_effect = p_effect)
    df_stats['TrueEffectSize'] = effects
    se_sample_sizes = np.power(sample_sizes, -1/2)
    #may need broadcast se_sample_sizes
    noise = rng.normal(loc = 0, scale = 1, size = (n_trials, n_studies))
    betas = effects.reshape(-1,1)@dilution.reshape(1,-1)
    noise =  noise @ np.diag(se_sample_sizes)

    for i, var in enumerate(beta_vars):
        df_stats[var] = betas[:,i] + noise[:,i]
        df_stats[z_vars[i]] = np.divide(df_stats[var], se_sample_sizes[i])
        df_stats[p_vars[i]] = 2*norm.cdf(-1*np.abs(df_stats[z_vars[i]]))
    df_stats[se_vars] = se_sample_sizes
    df_stats['TrueEffect'] = 1*(effects != 0)
    return df_stats

def nll_data(betas,ses,alpha):
    #assumes independence, but we can generalize
    #print(alpha)
    #alpha[0] = 1
    alpha = np.abs(np.array(alpha))
    alpha[alpha < .1] = .1
    alpha[alpha > 10] = 10
    #normalize
    alpha[0] = 1
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
