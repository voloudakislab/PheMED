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
python phemed.py --sum_stats data/sim_data.csv --n_studies 2 --out_file output/local_test
```
where --sum_stats identifies the csv with the merged log odds ratios and standard errors for each study.  --n_studies denotes the number of studies being analyzed and --out_file denotes the name of the output file.  For the sum_stats csv, phemed assumes that when analyzing N studies.  The last N columns refer to the standard errors of the respective ordered studies and the N columns preceding those refer to the log odds ratios.  (See the sim_data.csv in the data directory for an example.)

PheMED expects the input csv file to have columns: SNP, CHR, and POS corresponding to the rsid, chromosome and base pair position of the SNP.  

__Outputs__: By default, PheMED outputs a csv ```<out>_CI.csv``` estimating percentiles for the effective dilution for each study except for the first study in the list.  PheMED also produces a log file ```<out>.log```

For additional arguments run
```
python phemed.py -h
```
