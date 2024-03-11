---
title: 7_RunLDSC
description: "The script used to integrate LDSC."
sidebar_position: 1
---

# 7_RunLDSC

## Explanation

This is a script that will run 2 R scripts and also [LDSC](https://github.com/bulik/ldsc), which outputs information useful for determining if a ChromHMM model has biologically relevant states (not just statistically relevant).

The script achieves this using [partitioned heritability](https://www.nature.com/articles/ng.3404) and requires reference files in order to work properly. These files can be downloaded from [this online repository](https://zenodo.org/records/10515792). The reference files do not have to be from the 1000 genomes project (just any large collection of SNPs will do). The required files are:

- 1000 genomes (or similar) PLINK files (for your genome build)
- 1000 genomes (or similar) weight files (built from ldsc, best to get these from the online repository)
- Collection of GWAS traits in sumstats format (again, best to get these from the online repository)

## Example usage

```bash
# Generates heatmap of enrichments and bar plots for all GWAS-traits for the 
# model with 5 states in the model directory
sbatch 7_RunLDSC.sh \
--config="path/to/configuration/directory" \
--state=5
```

```bash
# Generates heatmap of enrichments and bar plots for all height based traits for  
# the model with 6 states in the model directory
sbatch 7_RunLDSC.sh \
--config="path/to/configuration/directory" \
--state=6 \
--gwas="height"
```