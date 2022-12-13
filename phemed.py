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
from numpy.random import Philox, Generator

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
parser.add_argument("--max_iters", type = int, default = 300,
            help = "Maximum iterations for optimization")
parser.add_argument("--optimizer_method", type = str, default = "Nelder-Mead",
            help = "Algorithm for optimization, see scipy.minimize for valid choices")
parser.add_argument("--block_window_initialize", type = int, default = 2000,
                        help = "Default value for finding block size for blocked bootstrap")
parser.add_argument("--bruteforce_start", type = boolean_string, default = True,
            help = "Use bruteforce to choose initial guess for MLE for estimating p-values using extreme value theory")
if __name__ == '__main__':

    args = parser.parse_args()
    stats_path = args.sum_stats
    snp_path = args.snp_list
    n_studies = args.n_studies
    merge_snps = args.merge_snps
    out_file = args.out
    compute_cis = args.compute_cis
    seed = args.seed
    n_trials = args.n_CIs
    max_iters = args.max_iters
    bruteforce = args.bruteforce_start
    optimizer_method = args.optimizer_method
    shift_range = args.block_window_initialize
    key = 2**96 + 2**33 + 2**17 + 2**9
    rng = Generator(Philox(key = key + seed))
    logging.basicConfig(filename= out_file + ".log",
                    format='%(asctime)s~%(levelname)s~%(message)s',
                    filemode='w')
    logging.captureWarnings(True)
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG) #lowest level, DEBUG, NOTSET
    logger.info("Running PheMED")
    #np.seterrcall(logger)
    #np.seterr(all='warn')
    command_string = """ phemed.py --out {out_file}
                             --sum_stats {sum_stats}
                             --n_studies {n_studies}
                             --snp_list {snp_list}
                             --merge_snps {merge_snps}
                             --seed {seed}
                             --compute_cis {compute_cis}
                             --n_CIs {n_trials}
                             --block_window_initialize {block_window}""".format(out_file = out_file,
                             sum_stats = stats_path, n_studies = n_studies,
                            snp_list = snp_path, merge_snps = merge_snps, seed = seed,
                            compute_cis = compute_cis, n_trials = n_trials, block_window = shift_range)
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
                        options={'maxiter': max_iters}, method = optimizer_method)
    #code normalizes first entry to be 1
        message = optimizer.message
        logger.info("PheMED " + message)
    #dilution values
        optimizer.x[0] = 1
        logger.info("Effective dilution values are : " + str(list(np.round(optimizer.x, 4))))
        if min(optimizer.x) < .2 or max(optimizer.x) > 5:
            logger.warning("Extreme dilution values detected.  Confirm that effect alleles are aligned across studies.")
        if compute_cis:
            if n_trials < 1000:
                logger.warning("Number of Bootstrap Samples appears small.  This can result in noisy CI and p-value outputs.")
            df_stats = df_stats.sort_values(by = ['CHR','POS']).reset_index(drop = True)
        #compute block size
            block_size_list = []
            for beta_var in beta_vars:
                corr_stats = boots.compute_autocorr(df_stats, beta_var, shift_range = shift_range)
                block = boots.compute_block_size(corr_stats, df_stats)
                if block >= shift_range:
                    logger.warning("""Initial block size provided is too small for convergence.
                    Consider rerunning with a larger block size with the block_window_initialize option.""")
                g, d = boots.compute_g_and_d(corr_stats, 2*block)
                beta_block_size = boots.compute_final_block_size(df_stats, g, d)
                block_size_list.append(beta_block_size)
            block_size_list.append(1)
            final_block_size = np.max(block_size_list)

            #np.random.seed(seed)
            columns = ['Sample'] + ['PheMed_' + str(i + 1) for i in range(n_studies) if i != 0] + ['Convergence']
            df_results_sim = pd.DataFrame(columns = columns)
            logger.info("Computing confidence intervals")
            for i in range(n_trials):
                df_results_sample = boots.circular_block_bootstrap(rng, df_stats,final_block_size)
    #df_meta_sample = df_meta.sample(n = df_meta.shape[0], replace = True)
                betas = df_results_sample[beta_vars]
                ses = df_results_sample[se_vars]
                optimizer = minimize(lambda x: phe.nll_data(betas,ses, x), np.ones(n_studies),
                        options={'maxiter': max_iters}, method = optimizer_method)
                df_results_sim.loc[i] = [i] + list(optimizer.x[1:n_studies]) +  [optimizer.message]

            df_ci = pd.DataFrame(df_results_sim.quantile([.025,.50,.975]))
            df_ci = df_ci.reset_index().rename(columns = {'index':'quantile'})
            logger.info("Saving confidence intervals")
            df_ci.to_csv(out_file + '_CIs.csv', index = False)


    #compute P-values

        index = 0
        if compute_cis:
            logger.info("Computing p-values")
            phemed_sim_vars = [var for var in df_results_sim.columns if 'PheMed_' in var]
            df_pvalue = pd.DataFrame(columns = ['PheMedStudy', 'PValue','Method','PassedQC'])
            for var in phemed_sim_vars:
                df_sim_mini = df_results_sim[['Sample',var]].rename(columns = {var:'Alpha'})
                #naive P-value computation
                p_naive, valid_test = boots.simple_p_test(df_sim_mini)
                df_pvalue.loc[index] = [var, p_naive, 'NaiveCount', valid_test]
                index += 1

                #normal P-value computation
                df_sim_mini = boots.norm_alpha(df_sim_mini)
                p_normal = boots.compute_normal_p(df_sim_mini)
                valid_test = boots.normality_check(df_sim_mini)
                df_pvalue.loc[index] = [var, p_normal, 'Normal', valid_test]
                index += 1

                #Extreme Value Theory
                top_n, k, a, alpha_vals, base_alpha, fit_found = boots.fit_gpd(df_sim_mini, bruteforce_initializer = bruteforce)
                valid_test = True
                if p_naive > top_n/df_sim_mini.shape[0]:
                    valid_test = 'P-Value is not an Extreme Value'
                    p_evt = np.nan
                    min_eval = 0 #default value to avoid bug
                elif not fit_found:
                    valid_test = 'Failed to optimize fit'
                    p_evt, min_eval = boots.compute_evt_p(rng, logger, df_sim_mini, top_n, a, k, n_sim = int(1e4))
                else:
                    p_evt, min_eval = boots.compute_evt_p(rng, logger, df_sim_mini, top_n, a, k, n_sim = int(1e4))

                if min_eval < 0:
                    valid_test += " Covariance matrix not positive definite, check logs"
                df_pvalue.loc[index] = [var, p_evt, 'Extreme Value Theory', valid_test]
                index += 1
                #export pvalues
            logger.info("Saving p-values")
            df_pvalue.to_csv(out_file + '_PVals.csv', index = False)

    except Exception as e:
        #ex_type, ex, tb = sys.exc_info()
        logger.error(str(e))
        raise
    finally:
        logger.info('Analysis finished')


#print results
