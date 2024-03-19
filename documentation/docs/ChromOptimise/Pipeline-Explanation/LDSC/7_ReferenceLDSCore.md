---
title: 7_RunLDSC
description: "Calculates the ldscores for annotation files."
sidebar_position: 1
---

# 7_RunLDSC

## Explanation

This is a script that will run 2 R scripts and also [LDSC](https://github.com/bulik/ldsc), which outputs information useful for determining if a ChromHMM model has biologically relevant states (not just statistically relevant).

Note that you still need to input the bin size, sample size and the number of models learned here. The reason for this is because the file structure is designed such that multiple runs of the same dataset can be analysed concurrently.

The script achieves this using [partitioned heritability](https://www.nature.com/articles/ng.3404) and requires reference files in order to work properly. These files can be downloaded from [this online repository](https://zenodo.org/records/10515792). The reference files do not have to be from the 1000 genomes project (just any large collection of SNPs will do). The required files are:

- 1000 genomes (or similar) PLINK files (for your genome build)
- 1000 genomes (or similar) weight files (built from ldsc, best to get these from the online repository)
- Collection of GWAS traits in sumstats format (again, best to get these from the online repository)

## Example usage

```bash
# Generates annotation and ldscore for each chromosome using the model with
# 5 states in the model directory
sbatch 7_RunLDSC.sh \
--config="path/to/configuration/directory" \
--state=5 \
--binsize=200 \
--samplesize=75 \
--nummodels=6 \
```

```bash
# Generates annotation and ldscore for each chromosome using the model with
# the optimum number of states (as decided by 6_OptimalNumberOfStates) 
sbatch 7_RunLDSC.sh \
--config="path/to/configuration/directory" \
--gwas="height" \
--binsize=200 \
--samplesize=75 \
--nummodels=6 \
```
