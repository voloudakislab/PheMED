import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from importlib import reload
import sys
import warnings
#from tqdm import tqdm_notebook
from scipy.optimize import minimize
from scipy.stats import norm
from scipy.stats import normaltest
from scipy.stats import shapiro
from scipy.stats import anderson
#p.seterr(all='warn')

def df_row_shift(df,start_index):
    if start_index > 0:
        new_start_index = start_index
    else:
        new_start_index = -1*start_index #as per definition in paper
        #new_start_index = df.shape[0] - np.abs(start_index)
    df_start = df.loc[new_start_index:df.shape[0]]
    #df_end = df.loc[0:new_start_index - 1]
    #all wrap around values should be replaced with mean so as not to contribute to R(k)
    #for var in df_end.columns:
    #    if df_end[var].dtypes == 'float64':
    #        df_end[var] = df[var].mean()
    #df_final = pd.concat([df_start,df_end])
    return df_start

def compute_autocorr(df, var, shift_range = 1000, start = -1):
    #will be buggy if you chooose too few rows for autocorr
    if start != 0:
        start = -1*shift_range
    df_corr = pd.DataFrame(columns = ['Shift','CovRatio','BaseCov'])
    mean_df = np.abs(df[var]).mean()
    cov_denominator = (1/df.shape[0])*np.dot(np.abs(df[var]) - mean_df,
                                          np.abs(df[var]) - mean_df)
    for i in np.arange(start,shift_range):
        shift = i + 1

        #print(df.shape)
        #print(shift)
        df_shifted = df_row_shift(df, shift)
        #print(df_shifted.shape)
        cov = (1/df.shape[0])*np.dot(np.abs(df_shifted[var]) - mean_df,
                                          np.abs(df.loc[0:df.shape[0] - np.abs(shift) - 1,var]) - mean_df) #use normalization constant 1/N as per paper
        rsq = cov/cov_denominator
        base_cov = cov
        df_corr.loc[i] = [shift, rsq,base_cov]
    return df_corr

def compute_block_size(df_corr, df_meta):
    solution = 'Not Found'
    start_index = 0
    df_corr_mini = df_corr.loc[df_corr['Shift'] >= 0].reset_index(drop = True)
    #print(df_corr_mini.head())
    while solution == 'Not Found':
        window = int(np.ceil(np.max([5, np.sqrt(np.log10(df_meta.shape[0]))])))
        #print(window)
        indep_threshold = 2*np.sqrt(np.log10(df_meta.shape[0])/df_meta.shape[0])
        score = (df_corr_mini.loc[start_index + 1: start_index + window + 1,'CovRatio']  < indep_threshold).mean()
        if score == 1:
            solution = 'Found'
        else:
            start_index += 1
    return start_index

def trap_function(x, total_distance = 10):
    if np.abs(x) > total_distance:
        return 0
    elif np.abs(x) <= total_distance/2:
        return 1
    else:
        return 2*(1 - np.abs(x/total_distance))

def compute_g_and_d(df_corr, block_size):
    df_corr['AbsShift'] = np.abs(df_corr['Shift'])
    df_corr_mini = df_corr.loc[df_corr['AbsShift'] <= block_size].copy() #fixes loc issue
    df_corr_mini['Lambda'] = df_corr_mini['Shift'].apply(lambda x: trap_function(x, total_distance = block_size))
    df_corr_mini['LambdaG'] = np.multiply(df_corr_mini['Lambda'], df_corr_mini['AbsShift'])
    #print(df_corr_mini.head())
    return [np.dot(df_corr_mini['LambdaG'], df_corr_mini['BaseCov']),
            (4/3)*np.dot(df_corr_mini['Lambda'], df_corr_mini['BaseCov'])**2]

def compute_final_block_size(df_meta, g, d):
    return int(np.ceil((df_meta.shape[0]*2*(g**2)/d)**(1/3)))

def circular_block_bootstrap(rng,df, block_size):
    #find how many blocks we are grabbing
    n_samples = int(np.floor(df.shape[0]/block_size))
    #compute starting point of blocks
    rand_rows = rng.integers(df.shape[0], size = (n_samples,1))
    #compute entire block explicitly
    selected_rows = rand_rows + np.multiply(np.ones((n_samples,block_size)), np.arange(0,block_size))
    #keep things in mod arithmetic
    selected_rows = np.remainder(selected_rows, df.shape[0])
    #choose selected rows
    df_sample = df.loc[selected_rows.flatten()].reset_index(drop = True)
    return df_sample

#functions for computing p-values
def simple_p_test(df):
    #check n_successes, n_failures
    high_dilution_shape = df.loc[df['Alpha'] > 1].shape[0]
    low_dilution_shape = df.loc[df['Alpha'] <= 1].shape[0]
    valid_test_flag = False
    if high_dilution_shape >= 10 and low_dilution_shape >= 10:
        valid_test_flag = True

    min_shape = np.min([high_dilution_shape, low_dilution_shape])
    return 2*min_shape/df.shape[0], valid_test_flag

def norm_alpha(df):
    if df['Alpha'].mean() > 1:
        df['Alpha'] = np.power(df['Alpha'], -1)
    return df

def compute_normal_p(data):
    mean = data['Alpha'].mean()
    sd = data['Alpha'].std()
    #print(mean, sd)
    z = (1 - mean)/sd
    return 2*(1 - norm.cdf(z))

def normality_check(data):
    passed_checks = True
    stat, p = shapiro(data['Alpha'])
    #print('Shapiro', stat, p)
    if p < .05:
        passed_checks = False
    stat, p = normaltest(data['Alpha'])
    #print('DAgostinos Normal Test', stat, p)
    if p < .05:
        passed_checks = False
    result = anderson(data['Alpha'])
    #print('Anderson Statistic: %.3f' % result.statistic)
    sl, cv = result.significance_level[2], result.critical_values[2]
    if result.statistic >= result.critical_values[2]:
        passed_checks = False
        #print('%.3f: %.3f, data looks normal (fail to reject H0)' % (sl, cv))
        #print('%.3f: %.3f, data does not look normal (reject H0)' % (sl, cv))
    return passed_checks

#EVT
def grab_top_values(df_results, top_n = 250):
    #make the extreme values increasing towards the null if consistently above
    if df_results['Alpha'].mean() > 1:
        df_temp = df_results.sort_values(by = 'Alpha', ascending = True).reset_index(drop = True)
        df_temp['Alpha'] = np.power(df_temp['Alpha'], -1)
    else:
        df_temp = df_results.sort_values(by = 'Alpha', ascending = False).reset_index(drop = True)

    base_alpha = df_temp.loc[top_n, 'Alpha']
    return df_temp.loc[0:top_n - 1,'Alpha'] - base_alpha, base_alpha

def compute_covar(a, k, n):
    return ((1-k)/n)*np.array([[2*(a**2), a],[a,(1 - k)]])

def gpd_tail_p(x,k,a):
    term = k*x/a
    if term > 1:
        return 0
    else:
        return 1 - (1 - (1 - term)**(1/k))

def compute_full_mle(alpha_vals, vals):
    k = vals[0]
    a = vals[1]
    a = np.max([1e-8,a, k*np.max(alpha_vals) + 1e-8])
    n = alpha_vals.shape[0]
    ratio = k/a
    mle = -n*np.log(a) - (1 - (1/k))*np.log(1 - ratio*alpha_vals).sum()
    return mle

def compute_anderson_darling(alpha_vals, k, a):
    #alpha_vals = alpha_vals[:-1]
    p_vals = alpha_vals.apply(lambda x: 1 - gpd_tail_p(x, k, a)).values
    #print(p_vals[0:5])
    reversed_p_vals = p_vals[::-1]
    #print(reversed_p_vals)
    n = p_vals.shape[0]
    #print(n)
    log_p_sum = np.log(reversed_p_vals) + np.log(1 - p_vals)
    factors = 2*(np.arange(n) + 1) - 1
    log_p_sum = np.multiply(factors, log_p_sum)
    #print(log_p_sum)
    stat = -n -(1/n)*np.sum(log_p_sum)
    return stat, p_vals


def grab_threshold(k):
    #grabs cutoff value at p = .05
    rounded_k = k #np.round(k, 1)
    if rounded_k < -.9:
        return .771
    elif rounded_k < -.5:
        return .83
    elif rounded_k < -.2:
        return .903
    elif rounded_k < -.1:
        return .935
    elif rounded_k < 0:
        return .974
    elif rounded_k < .1:
        return 1.02
    elif rounded_k < .2:
        return 1.074
    elif rounded_k < .3:
        return 1.14
    elif rounded_k < .4:
        return 1.221
    else:
        return 1.321


def fit_gpd(df, initial_vals = np.array([.2,.01])):
    fit_found = False
    top_n = 250

    #df_output = pd.DataFrame(columns = [''])
    while fit_found == False and top_n > 20:
        #print(top_n)
        alpha_vals, base_alpha = grab_top_values(df, top_n = np.min([top_n,df.shape[0] - 1]))
        opt = minimize(lambda x: -1*compute_full_mle(alpha_vals, x), initial_vals,
                            method = 'Nelder-Mead', options = {'maxiter':10000})
        k = opt.x[0]
        a = opt.x[1]
        #print('Optimization', k, a)
        anderson_threshold = grab_threshold(k)
        a = np.max([1e-8,a, k*np.max(alpha_vals) + 1e-6])
        anderson_stat, _ = compute_anderson_darling(alpha_vals, k, a)
        #print(anderson_stat)
        if anderson_stat < anderson_threshold:
            fit_found = True
            #print('Fit Found', top_n, anderson_stat, anderson_threshold)
        else:
            top_n = top_n - 10
    return top_n, k, a, alpha_vals, base_alpha, fit_found

def compute_evt_p(rng, df, n_tail, a, k, n_sim = int(1e4)):
    covar_matrix = compute_covar(a,k,n_tail)
    #print(a, k, covar_matrix, n_sim)
    sim_vars = rng.multivariate_normal([a,k], compute_covar(a,k,n_tail), size = n_sim)
    p_sims = np.zeros(n_sim)
    for i in range(n_sim):
        p_sims[i] = 2*(n_tail/df.shape[0])*gpd_tail_p(1, sim_vars[i,1], sim_vars[i,0])
    return np.mean(p_sims)
