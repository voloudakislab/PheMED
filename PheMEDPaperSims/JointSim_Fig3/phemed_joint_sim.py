import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from importlib import reload
import sys
import warnings
from scipy.stats import norm
from scipy.stats import chi2
from scipy.optimize import minimize
import argparse
from numpy.random import Philox, Generator
import phemed_def_helpers as phe_help

parser = argparse.ArgumentParser()
parser.add_argument('--ntests',  type=float, default = 1e5)
parser.add_argument('--true_p_effects', type=float, default = .1)
parser.add_argument('--ppv_gold', type=float)
parser.add_argument('--npv_gold', type=float)
parser.add_argument('--sample_size_gold', type=float, default = 2e4)
parser.add_argument('--sample_size_best', type=float, default = 2e4)
parser.add_argument('--ppv_diluted', type=float)
parser.add_argument('--npv_diluted', type=float)
parser.add_argument('--ppv_best', type=float, default = 1)
parser.add_argument('--npv_best', type=float, default = 1)
parser.add_argument('--sample_size_diluted', type=float, default = 2e4)
parser.add_argument('--job_id', type=int)
parser.add_argument('--effect_size', default = .025, type = float)
parser.add_argument('--n', default = 1e4, type = float)
args = parser.parse_args()

ppv_gold = args.ppv_gold
npv_gold = args.npv_gold
ppv_diluted = args.ppv_diluted
npv_diluted = args.npv_diluted
ppv_best = args.ppv_best
npv_best = args.npv_best
n_gold = args.n #args.sample_size_gold
n_diluted = args.n #args.sample_size_diluted
n_best = args.n #args.sample_size_best
p_effects = args.true_p_effects
n_tests = int(args.ntests)
job_id = args.job_id
random_seed = job_id
effect_size = args.effect_size
#now use Philox
key = 2**96 + 2**33 + 2**17 + 2**9
rng = Generator(Philox(key = key + random_seed))

df_alpha = pd.DataFrame(columns = ['Trial','Alpha','Method','ExpectedAlpha','Diluted_PPV','Diluted_NPV','Gold_PPV','Gold_NPV',
                                  'N_Gold','N_Diluted', 'N_Tests','P_Effects','Message','UsedBrute','JointAnalysis'])
alpha_index = 0
expected_alpha = (ppv_gold + npv_gold - 1)/(ppv_diluted + npv_diluted - 1)
df_diluted = phe_help.sim_case_status_table(n_diluted, n_diluted, ppv_diluted, npv_diluted)
df_gold = phe_help.sim_case_status_table(n_gold, n_gold, ppv_gold, npv_gold)
df_best = phe_help.sim_case_status_table(n_best, n_best, ppv_best, npv_best)
index = 0
for trial in range(2000): #halve samples
    maf_case, maf_control = phe_help.generate_probs(rng, p_effects = p_effects, n_tests = n_tests, true_effect_size = effect_size)
    df_diluted_counts = phe_help.generate_snp_counts(rng, df_diluted, maf_case, maf_control)
    df_gold_counts = phe_help.generate_snp_counts(rng, df_gold, maf_case, maf_control)
    df_best_counts = phe_help.generate_snp_counts(rng, df_best, maf_case, maf_control)
    diluted_stats = phe_help.compute_lor_se(df_diluted_counts)
    gold_stats = phe_help.compute_lor_se(df_gold_counts)
    best_stats = phe_help.compute_lor_se(df_best_counts)
    #Compute Dilution
    df_betas, df_ses = phe_help.extract_betas_ses3(diluted_stats, gold_stats, best_stats)
    optimizer = minimize(lambda x: phe_help.nll_data(df_betas,df_ses, x, uniqueness = True), np.array((1,1,1)),
                        options={'maxiter': 1000}, method = 'Nelder-Mead')
    dilution = 1/optimizer.x[0]
    used_brute = False
    if dilution < .1 or dilution > 10:
        rranges = (slice(.1,10,.1),slice(.1,10,.1), slice(1,1.1,0.1))
        optimizer_brute = brute(lambda x: phe_help.nll_data(df_betas,df_ses, x, uniqueness = True), rranges)
        dilution = 1/optimizer_brute[0]
        used_brute = True
    df_alpha.loc[index] = [trial, dilution, 'LL', expected_alpha, ppv_diluted, npv_diluted, ppv_gold, npv_gold, n_gold,
                                 n_diluted, n_tests, p_effects, optimizer.message, used_brute, True]
    index += 1
    df_betas, df_ses = phe_help.extract_betas_ses(diluted_stats, best_stats)
    optimizer = minimize(lambda x: phe_help.nll_data(df_betas,df_ses, x, uniqueness = True), np.array((1,1)),
                        options={'maxiter': 1000}, method = 'Nelder-Mead')
    dilution = 1/optimizer.x[0]
    used_brute = False
    if dilution < .1 or dilution > 10:
        rranges = (slice(.1,10,.1), slice(1,1.1,0.1))
        optimizer_brute = brute(lambda x: phe_help.nll_data(df_betas,df_ses, x, uniqueness = True), rranges)
        dilution = 1/optimizer_brute[0]
        used_brute = True
    df_alpha.loc[index] = [trial, dilution, 'LL', expected_alpha, ppv_diluted, npv_diluted, ppv_gold, npv_gold, n_gold,
                                 n_diluted, n_tests, p_effects, optimizer.message, used_brute, False]
    index += 1
    if trial % 100 == 0:
        df_alpha.to_csv('data/jointphemed_validate_results' + str(job_id) + '.csv', index = False)
df_alpha.to_csv('data/jointphemed_validate_results' + str(job_id) + '.csv', index = False)
