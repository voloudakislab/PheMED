# __FAQ for PheMED__
Phenotypic Measurement of Effective Dilution

### Q1: Why did I get a effective dilution value smaller than 1?
This means that your reference GWAS, corresponding to the log odds ratios that appears first in your input does not have the highest deltaP value = NPV + PPV - 1.  The analysis is correct, but if you want all of your phenotypic dilution values to be greater (or equal to) 1, this issue can be remedied by changing your reference GWAS to the GWAS that achieves the smallest effective dilution value.
