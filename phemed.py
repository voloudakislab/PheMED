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
parser.add_argument("--sample_sizes", type = str, default = "None Provided",
            help = "csv containing cases and controls for each study with columns Cases and Controls for each study.")
parser.add_argument("--eff_sample_sizes", type = str, default = "None Provided",
            help = "csv containing effective sample sizes for each study, with column N to denote the sample size for each study.")
parser.add_argument("--out", type = str, help = "name of output of PheMED")
parser.add_argument("--snp_list", type = str, default = "data/random_EUR_AFR.clumped",
help = "List of approximately independent SNPs to compute effective dilution.  The default value is to used the clumped SNPs in the data subdirectory.")
parser.add_argument("--merge_snps", type = boolean_string, default = True,
help = "Only retain SNPs that are approximately independent as defined by snp_list.  The default value is True.")
parser.add_argument("--compute_cis", type = boolean_string, default = True,
help = "Compute confidence intervals for estimates.  The default value is True.")
parser.add_argument("--seed", type = int, default = 18,
            help = "Seed for computing CIs.  The default value is 18.")
parser.add_argument("--n_CIs", type = int, default = 2000,
            help = "Number of bootstrap samples to compute CIs.  The default value is 2,000.")
parser.add_argument("--max_iters", type = int, default = 300,
            help = "Maximum iterations for optimization.  The default value is 300.")
parser.add_argument("--optimizer_method", type = str, default = "Nelder-Mead",
            help = "Algorithm for optimization, see scipy.minimize for valid choices.  The default value is Nelder-Mead.")
parser.add_argument("--block_window_initialize", type = int, default = 2000,
                        help = "Window for finding block size for blocked bootstrap.  The default value is 2,000.")
parser.add_argument("--dilution_limit", type = float, default = 10,
                        help = "Maximum allowed dilution value, must be larger than 1.")
parser.add_argument("--bruteforce_start", type = boolean_string, default = True,
            help = "Use bruteforce to choose initial guess for MLE for estimating p-values using extreme value theory.  The default value is True.")
if __name__ == '__main__':

    args = parser.parse_args()
    stats_path = args.sum_stats
    sample_sizes_path = args.sample_sizes
    eff_sample_sizes_path = args.eff_sample_sizes
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
    dilution_limit = args.dilution_limit
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
    #initialize outputs df
    df_outputs_summary = pd.DataFrame()
    #np.seterrcall(logger)
    #np.seterr(all='warn')
    sample_size_argument = "--sample_sizes"
    if sample_sizes_path == "None Provided":
        sample_size_argument = "--eff_sample_sizes"
        sample_size_path_to_use =  eff_sample_sizes_path
    else:
        sample_size_path_to_use =  sample_sizes_path
    command_string = """ phemed.py --out {out_file}
                             --sum_stats {sum_stats}
                             --n_studies {n_studies}
                             {sample_size_argument} {sample_size_path_to_use}
                             --snp_list {snp_list}
                             --merge_snps {merge_snps}
                             --seed {seed}
                             --compute_cis {compute_cis}
                             --n_CIs {n_trials}
                             --dilution_limit {dilution_limit}
                             --block_window_initialize {block_window}""".format(out_file = out_file,
                             sum_stats = stats_path, n_studies = n_studies, sample_size_argument = sample_size_argument,
                             sample_size_path_to_use = sample_size_path_to_use,
                            snp_list = snp_path, merge_snps = merge_snps, seed = seed,
                            compute_cis = compute_cis, n_trials = n_trials, dilution_limit = dilution_limit, block_window = shift_range)
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

        optimizer = minimize(lambda x: phe.nll_data(betas,ses, x, dilution_limit = dilution_limit), np.ones(n_studies),
                        options={'maxiter': max_iters}, method = optimizer_method)
    #code normalizes first entry to be 1
        message = optimizer.message
        logger.info("PheMED " + message)
    #dilution values
        optimizer.x[0] = 1
        #test if optimizer values exceed dilution limit
        if (np.min(optimizer.x) <= 1/dilution_limit) or (np.max(optimizer.x) >= dilution_limit):
             logger.warning("Very high (or low) values for the effective dilution were found.  Possibly consider different values for the dilution_limit parameter.")

        optimizer.x = np.clip(optimizer.x, 1/dilution_limit, dilution_limit)

        #Saving Dilution Vals
        df_dilution = pd.DataFrame(columns = ['StudyId','PheMED'])
        df_dilution['StudyId'] = np.arange(len(optimizer.x)) + 1
        df_dilution['PheMED'] = np.array(optimizer.x)

        df_outputs_summary['Study'] = df_dilution['StudyId'].apply(lambda x: 'Study ' + str(x))
        df_outputs_summary['Reference'] = False
        df_outputs_summary.loc[0, 'Reference'] = True
        df_outputs_summary['PheMED'] = df_dilution['PheMED']

        #df_dilution.to_csv(out_file + '_DilutionVals.csv', index = False)

        dilution_vals = np.round(optimizer.x, 4)
        logger.info("Effective dilution values are : " + str(list(dilution_vals)))
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
            columns = ['Sample'] + ['Study ' + str(i + 1) for i in range(n_studies) if i != 0] + ['Convergence']
            df_results_sim = pd.DataFrame(columns = columns)
            logger.info("Computing confidence intervals")
            for i in range(n_trials):
                df_results_sample = boots.circular_block_bootstrap(rng, df_stats,final_block_size)
    #df_meta_sample = df_meta.sample(n = df_meta.shape[0], replace = True)
                betas = df_results_sample[beta_vars]
                ses = df_results_sample[se_vars]
                optimizer = minimize(lambda x: phe.nll_data(betas,ses, x, dilution_limit = dilution_limit), np.ones(n_studies),
                        options={'maxiter': max_iters}, method = optimizer_method)
                df_results_sim.loc[i] = [i] + list(optimizer.x[1:n_studies]) +  [optimizer.message]

            df_ci = pd.DataFrame(df_results_sim.quantile([.025,.975]))
            df_ci = df_ci.reset_index().rename(columns = {'index':'quantile'})
            #logger.info("Saving confidence intervals")
            #df_ci.to_csv(out_file + '_CIs.csv', index = False)


            #save CIs to df_outputs_summary
            df_outputs_summary['CI_.025'] = np.nan
            df_outputs_summary['CI_.975'] = np.nan
            for var in df_ci.columns:
                if var != 'quantile':
                    low_ci = df_ci[var].values[0]
                    high_ci = df_ci[var].values[-1]
                    df_outputs_summary.loc[df_outputs_summary['Study'] == var, 'CI_.025'] = low_ci
                    df_outputs_summary.loc[df_outputs_summary['Study'] == var, 'CI_.975'] = high_ci
    #compute P-values

        index = 0
        if compute_cis:
            logger.info("Computing p-values")
            phemed_sim_vars = [var for var in df_results_sim.columns if 'Study ' in var]
            df_pvalue = pd.DataFrame(columns = ['Study', 'PValue','Method','PassedQC'])
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
            #logger.info("Saving p-values")
            #df_pvalue.to_csv(out_file + '_PVals.csv', index = False)
            #save df_pvalue to df_outputs_summary

            #get largest P-value that passes QC.
            df_top_p = df_pvalue.loc[df_pvalue['PassedQC'] == True].groupby('Study')[['PValue']].max().reset_index().rename(columns = {'index': 'Study'})
            df_outputs_summary['P'] = np.nan
            for study in df_top_p['Study'].unique():
                p_value_temp = df_top_p.loc[(df_top_p['Study'] == study), 'PValue'].values[0]
                df_outputs_summary.loc[df_outputs_summary['Study'] == study, 'P'] = p_value_temp

            for method in df_pvalue['Method'].unique():
                df_outputs_summary[method] = np.nan
                df_outputs_summary[method + '_qc'] = True
                for study in df_pvalue['Study'].unique():
                    p_value_temp = df_pvalue.loc[(df_pvalue['Study'] == study) & (df_pvalue['Method'] == method), 'PValue'].values[0]
                    qc_temp = df_pvalue.loc[(df_pvalue['Study'] == study) & (df_pvalue['Method'] == method), 'PassedQC'].values[0]
                    #print(study, p_value_temp, qc_temp)
                    df_outputs_summary.loc[df_outputs_summary['Study'] == study, method] = p_value_temp
                    df_outputs_summary.loc[df_outputs_summary['Study'] == study, method + '_qc'] = qc_temp
            #get best P-Value



    #use dilution values to compute effective sample size
        logger.info("Computing Dilution Adjusted Effective Sample Sizes")
        if sample_size_path_to_use == "None Provided":
            logger.info("No sample sizes provided, dilution adjusted sample size could not be computed.")

        else:
            logger.info("Normalizing Dilution Values so Smallest Dilution Value is 1")
            dilution_vals_norm = dilution_vals/min(dilution_vals)
            df_samples = pd.read_csv(sample_size_path_to_use)
            if sample_size_argument != "--eff_sample_sizes":
                #compute effective sample sizes
                df_samples['Den_Eff_N'] = np.divide(1, df_samples['Cases']) + np.divide(1, df_samples['Controls'])
                df_samples['N'] = np.divide(4, df_samples['Den_Eff_N'])

            df_sample_sizes = df_samples['N'].values
            logger.info("Unadjusted Effective Sample Sizes are: " + str(list(np.round(df_sample_sizes, 2))))
            dilution_adj_sample_sizes = np.divide(df_sample_sizes, np.power(dilution_vals_norm,2))
            dilution_adj_sample_sizes = np.round(dilution_adj_sample_sizes, 2)
            logger.info("Dilution Adjusted Effective Sample Sizes are: " + str(list(dilution_adj_sample_sizes)))
            df_outputs_summary['DilutionAdjNEff'] = dilution_adj_sample_sizes
        df_outputs_summary.to_csv(out_file + '_Summary.csv', index = False)

    except Exception as e:
        #ex_type, ex, tb = sys.exc_info()
        logger.error(str(e))
        raise
    finally:
        logger.info('Analysis finished')


#print results
