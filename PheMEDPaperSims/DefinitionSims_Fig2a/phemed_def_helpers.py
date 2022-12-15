import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from importlib import reload
import sys
from scipy.optimize import minimize
from scipy.stats import norm
from numpy.random import Philox, Generator

def sim_case_status_table(n_cases, n_controls, ppv, npv):
    #returns a dataframe breaking down the counts of labeled status, true stats
    df_stats = pd.DataFrame(columns = ['Labeled Status','True Status','Count'])
    df_stats.loc[0] = ['Labeled Case', 'True Case', np.round(n_cases*ppv)]
    df_stats.loc[1] = ['Labeled Case', 'True Control', np.round(n_cases*(1 - ppv))]
    df_stats.loc[2] = ['Labeled Control', 'True Control', np.round(n_controls*(npv))]
    df_stats.loc[3] = ['Labeled Control', 'True Case', np.round(n_controls*(1 - npv))]

    return df_stats

def generate_probs(rng, n_tests = int(1e5), true_effect_size = .03, p_effects = .1):
    n_tests = int(np.round(n_tests))
    maf_control = rng.uniform(.05,.5,int(n_tests))
    effect_size = np.zeros(int(n_tests))
    effect_size[0:int(p_effects*n_tests)] = rng.normal(0,true_effect_size, int(p_effects*n_tests))
    #effect = delta_p*(1/p + 1/(1 - p))
    hmean = np.divide(1, maf_control) + np.divide(1, 1 - maf_control)
    delta_p = np.divide(effect_size, hmean)
    maf_case = maf_control + delta_p
    return maf_case, maf_control

def generate_snp_counts(rng, df_counts, maf_case, maf_control):
    df_snp_list = []
    df_counts['index'] = 0
    for row in range(df_counts.shape[0]):
        df_counts.loc[row, 'index'] = row
        count = df_counts['Count'].values[row]
        status = df_counts['True Status'].values[row]
        if status == 'True Case':
            df_temp = pd.DataFrame(rng.binomial(count, maf_case)).transpose()
        else:
            df_temp = pd.DataFrame(rng.binomial(count, maf_control)).transpose()
        #print(df_temp.shape)
        df_temp.columns = ['SNP_' + str(i) for i in range(len(df_temp.columns))]
        #print(df_temp)
        df_temp['index'] = row
        df_snp_list.append(df_temp)
    df_snp_counts = pd.concat(df_snp_list)
    #print(df_snp_counts.shape)
    df_counts = pd.merge(df_counts, df_snp_counts)

    return df_counts

def compute_lor_se(df_counts):
    df_counts['Count'] = pd.to_numeric(df_counts['Count']) #weird bug that its an object type
    df_sum = df_counts.groupby(by = 'Labeled Status').sum().reset_index()
    df_sum = df_sum.sort_values(by = 'Labeled Status') #ensures cases start first
    snp_cols = [var for var in df_sum.columns if var.startswith('SNP')]
    df_stats = df_sum[snp_cols].transpose()
    df_stats.columns = ['LCaseSNPCount', 'LControlSNPCount']
    case_count = df_sum['Count'].values[0]
    control_count = df_sum['Count'].values[1]
    df_stats['LCaseNoSNPCount'] = case_count - df_stats['LCaseSNPCount']
    df_stats['LControlNoSNPCount'] = control_count - df_stats['LControlSNPCount']
    ses = np.sqrt(np.sum(np.power(1.0*df_stats.values,-1), axis = 1))
    betas = np.log(np.divide(df_stats['LCaseSNPCount'],
                         df_stats['LCaseNoSNPCount'])) - np.log(np.divide(df_stats['LControlSNPCount'],
                                                        df_stats['LControlNoSNPCount']))
    df_stats['Beta'] = betas
    df_stats['SE'] = ses
    return df_stats

def extract_betas_ses(diluted_stats, gold_stats):
    ses_diluted = diluted_stats['SE'].values
    ses_gold = gold_stats['SE'].values
    betas_diluted = diluted_stats['Beta'].values
    betas_gold = gold_stats['Beta'].values
    df_betas = pd.DataFrame(columns = ['Beta1','Beta2'])
    df_betas['Beta1'] = betas_gold
    df_betas['Beta2'] = betas_diluted
    df_ses = pd.DataFrame(columns = ['SE1','SE2'])
    df_ses['SE1'] = ses_gold
    df_ses['SE2'] = ses_diluted

    return df_betas, df_ses

def extract_betas_ses3(diluted_stats, gold_stats, best_stats):
    ses_diluted = diluted_stats['SE'].values
    ses_gold = gold_stats['SE'].values
    ses_best = best_stats['SE'].values
    betas_diluted = diluted_stats['Beta'].values
    betas_gold = gold_stats['Beta'].values
    betas_best = best_stats['Beta'].values
    df_betas = pd.DataFrame(columns = ['Beta1','Beta2','Beta3'])
    df_betas['Beta1'] = betas_gold
    df_betas['Beta2'] = betas_diluted
    df_betas['Beta3'] = betas_best
    df_ses = pd.DataFrame(columns = ['SE1','SE2','SE3'])
    df_ses['SE1'] = ses_gold
    df_ses['SE2'] = ses_diluted
    df_ses['SE3'] = ses_best

    return df_betas, df_ses

def nll_data(betas,ses,alpha, uniqueness = True):
    #assumes independence, but we can generalize
    #print(alpha)
    #alpha[0] = 1
    alpha = np.abs(np.array(alpha))
    if uniqueness:
        alpha[-1] = 1
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
