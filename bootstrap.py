import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from importlib import reload
import sys
import warnings
from tqdm import tqdm_notebook
from scipy.optimize import minimize
from scipy.stats import norm



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
    for i in tqdm_notebook(np.arange(start,shift_range)):
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
    df_corr_mini = df_corr.loc[df_corr['AbsShift'] <= block_size]
    df_corr_mini['Lambda'] = df_corr_mini['Shift'].apply(lambda x: trap_function(x, total_distance = block_size))
    df_corr_mini['LambdaG'] = np.multiply(df_corr_mini['Lambda'], df_corr_mini['AbsShift'])
    print(df_corr_mini.head())
    return [np.dot(df_corr_mini['LambdaG'], df_corr_mini['BaseCov']),
            (4/3)*np.dot(df_corr_mini['Lambda'], df_corr_mini['BaseCov'])**2]

def compute_final_block_size(df_meta, g, d):
    return int(np.ceil((df_meta.shape[0]*2*(g**2)/d)**(1/3)))

def circular_block_bootstrap(df, block_size):
    #find how many blocks we are grabbing
    n_samples = int(np.floor(df.shape[0]/block_size))
    #compute starting point of blocks
    rand_rows = np.random.randint(df.shape[0], size = (n_samples,1))
    #compute entire block explicitly
    selected_rows = rand_rows + np.multiply(np.ones((n_samples,block_size)), np.arange(0,block_size))
    #keep things in mod arithmetic
    selected_rows = np.remainder(selected_rows, df.shape[0])
    #choose selected rows
    df_sample = df.loc[selected_rows.flatten()].reset_index(drop = True)
    return df_sample
