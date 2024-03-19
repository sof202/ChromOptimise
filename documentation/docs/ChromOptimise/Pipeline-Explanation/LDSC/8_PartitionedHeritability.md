---
title: 8_PartitionedHeritability
description: "Visualises the partitioned heritability (enrichment)"
sidebar_position: 3
---

# 8_PartitionedHeritability

## Explanation

This script will run LDSC again, but now with the `-h2` flag. This signifies
that partitioned heritability is wanted. Specifically, what we want to extract
is the enrichment for each of our states (and the other baseline annotations).
Enrichment is described as the proportion of heritability divided by the
proportion of SNPs in the annotated regions in 
[this paper](https://www.biorxiv.org/content/10.1101/014241v1.full.pdf).

After extracting this information, the script will call the Rscript
[HeritabilityPlots.R](./HeritabilityPlots.md) to help with visualisation.


## Usage

Typically, this script will be ran after being called by 
[7_ReferenceLDSCore.sh](./7_ReferenceLDSCore.md). However, one may want to use
this script in isolation (perhaps for testing a different set of gwas traits).

The reason why this script is split away from the previous script in the
pipeline is because the previous script is ran as an array. The memory assigned
to each task in the array is split evenly and is not dynamic. As a result,
unless your machine/system has RAM in the 100s of GBs, the partitioned
heritability step will return out of memory errors. As we only need to compute
the plots once, we just use a separate script to run these analyses.

## Example usage

```bash
# Generates heatmap of enrichments and bar plots for all GWAS-traits for the 
# model in the BinSize_200_SampleSize_75_6 directory
sbatch 7_RunLDSC.sh \
--config="path/to/configuration/directory" \
--binsize=200 \
--samplesize=75 \
--nummodels=6 \
```

