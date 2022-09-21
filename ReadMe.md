# PheMED
Phenotypic Measurement of Effective Dilution

### Dependencies
After installing [Anaconda](https://store.continuum.io/cshop/anaconda/), running the following commands will create an environment suitable to run PheMED.
```
conda env create --file environment.yml
source activate phemed
```
### Running PheMED
To run PheMED on the sample data, run the following command:
```
python phemed.py --sum_stats data/sim_data.csv --n_studies 2 --out output/local_test
```
where --sum_stats identifies the csv with the merged log odds ratios and standard errors for each study.  --n_studies denotes the number of studies being analyzed and --out denotes the name of the output file.  For the sum_stats csv, phemed assumes that when analyzing N studies, the last N columns refer to the standard errors of the respective ordered studies and the N columns preceding those refer to the log odds ratios. Furthermore, the first study listed is used as the reference study when measuring effective dilution. (See the sim_data.csv in the data directory for an example.)

PheMED expects the input csv file to have columns: SNP, CHR, and POS corresponding to the rsid, chromosome and base pair position of the SNP.  

__Understanding Outputs__: By default, PheMED outputs a csv ```<out>_CI.csv``` estimating percentiles 2.5%, 50% and 97.5% for the effective dilution for each study except for the first study in the list.  (e.g. In the CI file, PheMed_2 refers to the effective dilution for the second GWAS study listed in the input file.)  PheMED  also produces a log file ```<out>.log``` and a p-value file ```<out>_PVals.csv```.  For examples of output, see the output directory.   

__P Values__: For the p-value file, results are printed in tidy format, where we leverage three different p-value methodologies.  For each of the methodologies, there is a PassedQC column that indicates if it is appropriate to use that methodology to compute p-value.  As such, we do not require that all three methodologies pass the QC check; instead, only one such methodology needs to pass the QC check to estimate the p-value.

Nevertheless, for naive count, if the p-value does not pass the QC check (e.g. the p-value is very small and becomes hard to estimate from the bootstrap simulation), if the number of bootstrap samples is sufficiently large (e.g. = 2000, the default), we can still use the confidence intervals to infer if p < .05.  See our paper for details.  

__Brief Overview of Other Arguments__
Computing the CIs when measuring the effective dilution between two GWA studies can take around an hour.  If you would like to run the pipeline without computing the CIs, use ```--compute_cis False ``` or the ```--n_CIs <smaller integer>``` to reduce the number of bootstrap samples; however decreasing the number of samples below the default value can result in noisy p-values.

For additional arguments run
```
python phemed.py -h
```

__FAQ__: For FAQ and troubleshooting, please see the [FAQ here](https://github.com/voloudakislab/phemed/tree/main/faq)

__Citing our Paper__: TBA.  
