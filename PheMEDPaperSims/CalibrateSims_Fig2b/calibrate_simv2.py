#depends on vmupa, boots
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from importlib import reload
import bootstrap as boots
import validating_mupav4 as vmupa
import warnings
from scipy.optimize import minimize
import argparse
from numpy.random import Philox, Generator


parser = argparse.ArgumentParser()
parser.add_argument('--ntests',  type=float, default = 1e5)
parser.add_argument('--true_p_effects', type=float, default = .1)
#parser.add_argument('--dilution', type=float, nargs = '+')
#parser.add_argument('--sample_sizes', type=float, nargs = '+')
parser.add_argument('--job_id', type=int)
parser.add_argument('--effect_size', default = .02, type = float)
args = parser.parse_args()

dilution = np.array([1.0,1.0]) #np.array(args.dilution)
true_p_effects = args.true_p_effects
sample_sizes = np.array([2e4,2e4]) #np.array(args.sample_sizes)
n_tests = args.ntests
job_id = args.job_id
random_seed = job_id
effect_size = args.effect_size


key = 2**96 + 2**33 + 2**17 + 2**9
rng = Generator(Philox(key = key + random_seed))
df_final = pd.DataFrame(columns = ['Trial', 'SimP','SimValid','NormalP','NormalValid','EstimatedPheMED',
                                   'Converged','BootsMessage','BootstrapConverged','MedianBoots'])
df_results_list = []
for trial in range(10): #100
    jobid = 100*job_id + trial
    gwas_stats = vmupa.generate_gwas_stats(rng, effect_size, n_tests, true_p_effects, dilution, sample_sizes)
    beta_vars = [var for var in gwas_stats.columns if 'BETA' in var]
    sd_vars = [var for var in gwas_stats.columns if 'SE' in var]
    betas = gwas_stats[beta_vars]
    ses = gwas_stats[sd_vars]
    optimizer = minimize(lambda x: vmupa.nll_data(betas,ses, x), np.array((1,1)),
                        options={'maxiter': 300}, method = 'Nelder-Mead')
    estimated_phemed = optimizer.x[1]
    converged = optimizer.message
    df_results = pd.DataFrame(columns = ['JobId','BootstrapTrial','BootsPheMED', 'Converged','EstimatedPheMED'])
    for i in range(2000):
        indices = rng.integers(0, gwas_stats.shape[0], gwas_stats.shape[0])
        new_gwas_stats = pd.DataFrame(gwas_stats.loc[indices])
        betas = new_gwas_stats[beta_vars]
        ses = new_gwas_stats[sd_vars]
        optimizer = minimize(lambda x: vmupa.nll_data(betas,ses, x), np.array((1,1)),
                        options={'maxiter': 300}, method = 'Nelder-Mead')
        df_results.loc[i] = [jobid, i, optimizer.x[1], optimizer.message, estimated_phemed]

    df_results = df_results.rename(columns = {'BootsPheMED': 'Alpha'})
    simple_p, valid_flag = boots.simple_p_test(df_results)
    df_results = boots.norm_alpha(df_results)
    normal_p = boots.compute_normal_p(df_results)
    normal_valid = boots.normality_check(df_results)
    boots_message = df_results['Converged'].value_counts().index[0]
    converged_count = df_results['Converged'].value_counts().max()
    median_boots = df_results['Alpha'].median()
    df_results_list.append(df_results)
    df_final.loc[trial] = [jobid, simple_p, valid_flag, normal_p, normal_valid,
                           estimated_phemed, converged, boots_message, converged_count,median_boots]
    if trial % 10 == 0:
        if trial != 0:
            df_final.to_csv('data/pvalue_calibrate' + str(job_id) + '.csv', index = False)
df_final.to_csv('data/pvalue_calibrate' + str(job_id) + '.csv', index = False)
