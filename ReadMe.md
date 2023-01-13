# PheMED
Phenotypic Measurement of Effective Dilution

### Dependencies
After installing [Anaconda](https://store.continuum.io/cshop/anaconda/) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html), running the following commands will create and activate an environment suitable to run PheMED.
```
conda env create --file environment.yml
source activate phemed
```
### Running PheMED
To run PheMED on the sample data, run the following command:
```
python phemed.py --sum_stats data/sim_data.csv --n_studies 2 --out output/local_test
```
- `--sum_stats` identifies the csv with the merged log odds ratios and standard errors for each study
- `--n_studies` denotes the number of studies being analyzed
- `--out` denotes the name of the output file  

For the sum_stats csv, PheMED expects the input csv file to have columns: `SNP`, `CHR`, and `POS` corresponding to the rsid, chromosome and base pair position of the SNP. The columns that follow are assumed to contain effect sizes {log(OR)<sub>1</sub>,log(OR)<sub>2</sub>,...,log(OR)<sub>n</sub>} followed by the SEs {SE<sub>1</sub>,SE<sub>2</sub>,...,SE<sub>3</sub>}, when analyzing n studies. Furthermore, the first study listed is used as the reference study when measuring effective dilution. (See the sim_data.csv in the data directory for 2 studies or an example with three studies below.)

| SNP | CHR | POS | STUDY1 | STUDY2  | STUDY3 | SE1   | SE2  | SE3   |
|-----|-----|-----|--------|---------|--------|-------|------|-------|
| rs1 | 1   | 1   | 0.069  | 0.010   | 0.040  | 0.044 | 0.01 | 0.027 |
| rs2 | 1   | 2   | 0.038  | 0.015   | 0.026  | 0.044 | 0.01 | 0.027 |
| rs3 | 1   | 3   | 0.045  | -0.0004 | 0.022  | 0.044 | 0.01 | 0.027 |
| rs4 | 1   | 4   | -0.074 | -0.0107 | -0.042 | 0.044 | 0.01 | 0.027 |

__Understanding Outputs__: By default, PheMED outputs a csv ```<out>_CI.csv``` estimating percentiles 2.5%, 50% and 97.5% for the effective dilution for each study except for the first study in the list.  (e.g. In the CI file, PheMed_2 refers to the effective dilution for the second GWAS study listed in the input file.)  PheMED  also produces a log file ```<out>.log```, a p-value file ```<out>_PVals.csv``` and a file containing effective dilution estimates ```<out>_DilutionVals.csv```.  For examples of output, see the output directory.   

__P Values__: For the p-value file, results are printed in tidy format, where we leverage three different p-value methodologies.  For each of the methodologies, there is a PassedQC column that indicates if it is appropriate to use that methodology to compute p-value.  As such, we do not require that all three methodologies pass the QC check; instead, only one such methodology needs to pass the QC check to estimate the p-value.

Nevertheless, for naive count, if the p-value does not pass the QC check (e.g. the p-value is very small and becomes hard to estimate from the bootstrap simulation), if the number of bootstrap samples is sufficiently large (e.g. = 2000, the default), we can still use the confidence intervals to infer if p < .05.  See our paper for details.  

__Brief Overview of Other Arguments__
Computing the CIs when measuring the effective dilution between two GWA studies can take around an hour.  If you would like to run the pipeline without computing the CIs, use ```--compute_cis False ``` or the ```--n_CIs <smaller integer>``` to reduce the number of bootstrap samples; however decreasing the number of samples below the default value can result in noisy p-values.  While sample sizes are not required to run PheMED, if you would like to compute a dilution adjusted effective sample size, you will need a csv indicating the sample sizes for each study maintained in the same order as the input file and add either the argument ```--eff_sample_sizes <path_to_file>``` or ```--sample_sizes <path_to_file>```, where the former argument expects a CSV with a column labeled N detailing the effective sample sizes and the later argument expects a CSV with columns labeled Cases and Controls, that detail the number of cases and controls for each study.  For examples of these files, see ```data/eff_sample_size_sim.csv``` and  ```data/sample_size_sim.csv``` respectively.  

For additional arguments run
```
python phemed.py -h
```
### Performing Dilution Adjusted Weights (DAW) Meta-Analysis with PheMED
After running the phemed script above, users can run a dilution adjusted weights meta-analysis with the following command
```
python daw_meta.py --sum_stats <path to summary stats file>
                  --n_studies <number of studies>
                  --out <prefix for output files>
                  --dilution_weights <path to DilutionVals csv outputted from the phemed script>
```
The script will then output a csv containing meta-analyzed Z-scores, p-values and effect sizes as measured according to the reference study.
Example output below: 
| SNP | CHR | POS | Z_META             | BETA_META            | P_META               |
|-----|-----|-----|--------------------|----------------------|----------------------|
| rs1 | 1   | 1   | 2.216573171014446  | 0.04237889772244553  | 0.026652272616020767 |
| rs2 | 1   | 2   | 2.001733259609271  | 0.03827135065330354  | 0.04531342729038265  |
| rs3 | 1   | 3   | 0.8086750926075479 | 0.015461144927879799 | 0.41870205861603993  |

__FAQ__: For FAQ and troubleshooting, please see the [FAQ here](https://github.com/voloudakislab/phemed/tree/main/faq)

__Support__: If you have any outstanding questions not addressed by the docs or FAQ, feel free to post your question as a [GitHub Issue here](https://github.com/voloudakislab/PheMED/issues) or send an email to david.burstein2 {at} mssm.edu

__Citing our Paper__: TBA.  
