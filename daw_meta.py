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

parser = argparse.ArgumentParser()
parser.add_argument("--sum_stats", type = str,
            help = "Path to csv containing all summary stats")
parser.add_argument("--n_studies", type = int,
            help = "Number of studies being analyzed")
parser.add_argument("--dilution_weights", type = str,
            help = "Dilution weights file")
parser.add_argument("--out", type = str, help = "MetaAnalysis Output")

if __name__ == '__main__':

    args = parser.parse_args()
    stats_path = args.sum_stats
    n_studies = args.n_studies
    out_file = args.out
    weights_file = args.dilution_weights
    logging.basicConfig(filename= out_file + "_meta.log",
                    format='%(asctime)s~%(levelname)s~%(message)s',
                    filemode='w')

    logging.captureWarnings(True)
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG) #lowest level, DEBUG, NOTSET
    logger.info("Running Dilution Adjusted Weights Meta-Analysis")

    command_string = """ daw_meta.py --out {out_file}
                             --sum_stats {sum_stats}
                             --n_studies {n_studies}
                             --dilution_weights {weights}
                             """.format(out_file = out_file,
                             sum_stats = stats_path, n_studies = n_studies, weights = weights_file)
    logger.info(command_string)
    try:
        df_stats = pd.read_csv(stats_path)

        beta_vars = [var for i, var in enumerate(df_stats.columns) if i < df_stats.shape[1] - n_studies and i >= df_stats.shape[1] - 2*n_studies]
        se_vars = [var for i, var in enumerate(df_stats.columns) if i >= df_stats.shape[1] - n_studies]
        #addressing nulls
        df_stats[beta_vars] = df_stats[beta_vars].fillna(0)
        df_stats[se_vars] = df_stats[se_vars].fillna(1000)

        #later add qc for betas, ses
        betas = df_stats[beta_vars].values
        ses = df_stats[se_vars].values

        df_weights = pd.read_csv(weights_file)
        #multiply betas, SEs by dilution values
        weights = df_weights['PheMED'].values
        logger.info("Performing Re-Weighting")
        betas = betas @ np.diag(weights)
        ses = ses @ np.diag(weights)

        def inv_var_meta(betas, ses, get_betas_only = False):
            #alpha[0] = 1
            #betas = betas @ np.diagonal(alpha)
            #ses = ses @ np.diagonal(alpha)
            inv_var = np.power(ses, -2)
            #return a vector
            #df_meta['NUM'] = df_meta['INV_VAR_' + diluted_key]*(df_meta['BETA_' + diluted_key]) + (df_meta['INV_VAR_' + base_key])*(df_meta['BETA_' + base_key])
            num = np.multiply(betas, inv_var).sum(axis = 1)
            if get_betas_only:
                denominator = inv_var.sum(axis = 1)
            else:
                denominator = np.sqrt(inv_var.sum(axis = 1))
            meta = np.divide(num, denominator)
            return meta
        logger.info("Performing Meta-Analysis")
        df_stats['Z_META'] = inv_var_meta(betas, ses)
        df_stats['BETA_META'] = inv_var_meta(betas, ses, get_betas_only = True)
        df_stats['P_META'] = 2*norm.cdf(-np.abs(df_stats['Z_META']))
        output_vars = [var for var in df_stats.columns if var not in beta_vars + se_vars]
        df_output = df_stats[output_vars]
        logger.info("Saving Results")
        df_output.to_csv(out_file + "_meta_results.csv", index = False)
    except Exception as e:
        #ex_type, ex, tb = sys.exc_info()
        logger.error(str(e))
        raise
    finally:
        logger.info('Analysis finished')
