import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from importlib import reload
import sys
import warnings
from scipy.stats import norm
from scipy.stats import chi2
import phemed_helpers as phe
from scipy.optimize import minimize
import argparse
import logging
import bootstrap as boots

#assumes betas and ses are last columns of df
#ID cols ... Beta Cols ... SE Cols

#add max iter parameter for optimize

#To Do:
#Add CIs, p-values
def boolean_string(s):
    if s not in {'False', 'True'}:
        raise ValueError('Not a valid boolean string for argument merge_snps')
    return s == 'True'


parser = argparse.ArgumentParser()
parser.add_argument("--sum_stats", type = str,
            help = "Path to csv containing all summary stats")
parser.add_argument("--n_studies", type = int,
            help = "Number of studies being analyzed")
parser.add_argument("--out", type = str, help = "name of output of PheMED")
parser.add_argument("--snp_list", type = str, default = "data/random_EUR_AFR.clumped",
help = "List of approximately independent SNPs to compute effective dilution")
parser.add_argument("--merge_snps", type = boolean_string, default = True,
help = "Only retain SNPs that are approximately independent as defined by snp_list.")
parser.add_argument("--compute_cis", type = boolean_string, default = True,
help = "Compute confidence intervals for estimates.")
parser.add_argument("--seed", type = int, default = 18,
            help = "Seed for computing CIs")
parser.add_argument("--n_CIs", type = int, default = 2000,
            help = "Number of bootstrap samples to compute CIs")

args = parser.parse_args()
stats_path = args.sum_stats
snp_path = args.snp_list
n_studies = args.n_studies
merge_snps = args.merge_snps
out_file = args.out
compute_cis = args.compute_cis
seed = args.seed
n_trials = args.n_CIs

logging.basicConfig(filename= out_file + ".log",
                    format='%(asctime)s~%(levelname)s~%(message)s',
                    filemode='w')

logger = logging.getLogger()
logger.setLevel(logging.DEBUG)
logger.info("Running PheMED")
command_string = """ phemed.py --out_file {out_file}
                             --sum_stats {sum_stats}
                             --n_studies {n_studies}
                             --snp_list {snp_list}
                             --merge_snps {merge_snps}
                             --seed {seed}
                             --compute_cis {compute_cis}
                             --n_CIs {n_trials}""".format(out_file = out_file,
                             sum_stats = stats_path, n_studies = n_studies,
                            snp_list = snp_path, merge_snps = merge_snps, seed = seed,
                            compute_cis = compute_cis, n_trials = n_trials)
logger.info(command_string)
try:
    df_stats = pd.read_csv(stats_path)




    #read in SNPs
    with open(snp_path) as f:
        lines = f.readlines()
    valid_snps = ['rs' + line.split('rs')[-1].split(' ')[0] for line in lines]

    initial_shape = df_stats.shape[0]
    logger.info("Number of rows in summary stats: {n_rows}".format(n_rows = initial_shape))
    if merge_snps:
        df_stats = df_stats.loc[df_stats['SNP'].isin(valid_snps)].reset_index(drop = True)
    post_merge_shape = df_stats.shape[0]
    logger.info("Number of valid SNPs in summary stats: {n_rows}".format(n_rows = post_merge_shape))
    if post_merge_shape < 1e5:
        logger.warning("Number of valid SNPs in summary stats appears low")

    beta_vars = [var for i, var in enumerate(df_stats.columns) if i < df_stats.shape[1] - n_studies
    and i >= df_stats.shape[1] - 2*n_studies]
    se_vars = [var for i, var in enumerate(df_stats.columns) if i >= df_stats.shape[1] - n_studies]

    #addressing nulls
    df_stats[beta_vars] = df_stats[beta_vars].fillna(0)
    df_stats[se_vars] = df_stats[se_vars].fillna(1000)

    #later add qc for betas, ses
    betas = df_stats[beta_vars]
    ses = df_stats[se_vars]

    optimizer = minimize(lambda x: phe.nll_data(betas,ses, x), np.ones(n_studies),
                        options={'maxiter': 300}, method = 'Nelder-Mead')
    #code normalizes first entry to be 1
    message = optimizer.message
    logger.info("PheMED " + message)
    #dilution values
    optimizer.x[0] = 1
    logger.info("Effective dilution values are : " + str(list(np.round(optimizer.x, 4))))

    if compute_cis:
        df_stats = df_stats.sort_values(by = ['CHR','POS']).reset_index(drop = True)
        #compute block size
        block_size_list = []
        for beta_var in beta_vars:
            corr_stats = boots.compute_autocorr(df_stats, beta_var, shift_range = 500)
            block = boots.compute_block_size(corr_stats, df_stats)
            g, d = boots.compute_g_and_d(corr_stats, 2*block)
            beta_block_size = boots.compute_final_block_size(df_stats, g, d)
            block_size_list.append(beta_block_size)
        block_size_list.append(1)
        final_block_size = np.max(block_size_list)

        np.random.seed(seed)
        columns = ['Sample'] + ['PheMed_' + str(i + 1) for i in range(n_studies) if i != 0] + ['Convergence']
        df_results_sim = pd.DataFrame(columns = columns)
        logger.info("Computing confidence intervals")
        for i in range(n_trials):
            df_results_sample = boots.circular_block_bootstrap(df_stats,final_block_size)
    #df_meta_sample = df_meta.sample(n = df_meta.shape[0], replace = True)
            betas = df_results_sample[beta_vars]
            ses = df_results_sample[se_vars]
            optimizer = minimize(lambda x: phe.nll_data(betas,ses, x), np.ones(n_studies),
                        options={'maxiter': 300}, method = 'Nelder-Mead')
            df_results_sim.loc[i] = [i] + list(optimizer.x[1:n_studies]) +  [optimizer.message]

        df_ci = pd.DataFrame(df_results_sim.quantile([.025,.50,.975]))
        logger.info("Saving confidence intervals")
        df_ci.to_csv(out_file + '_CIs.csv', index = False)

    #compute P-values


except Exception as e:
        #ex_type, ex, tb = sys.exc_info()
    logger.error(str(e))
    raise
finally:
    logger.info('Analysis finished')


#print results
