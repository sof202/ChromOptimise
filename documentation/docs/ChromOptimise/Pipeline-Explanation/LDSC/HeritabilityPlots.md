---
title: HeritabilityPlots
description: "The script used to create plots of partitioned heritability."
sidebar_position: 4
---

# HeritabilityPlots

## Explanation

This script will take in the `.results` files that are outputted from ldsc's 
partitioned heritability and visualise them for you. The two output files are:

1) A heatmap of the enrichment (Proportion of heritability for state / 
proportion of SNPS for state) for each state against each GWAS trait selected
2) A barplot for the p-values of the enrichment scores for each state for each 
GWAS trait selected 

This is the main output of runLDSC.sh that enables the user to interpret the 
biological relevance of their states. If lots of states are showing very 
little enrichment for any of the traits that you are considering, this 
suggests that the model is too specific still (or there is a problem with your 
input data).

## Warning message

Sometimes, this script will output the file WARNING.txt into the root of the
plots folder. This file will tell the user that their enrichment heatmap has
a high proportion of negative enrichments. Mathematically a negative enrichment
is nonsensicle (it is defined as the quotient of two scrictly positive values).
However, ldsc uses an unbiased estimator of $r^2$ due to the biases when
running over a large number of SNPS (which can be negative). The creators of
ldsc state that:

> If more than a few percent of all SNPs have negative LD Scores, this 
probably means that something is wrong, either that the sample size used for 
estimating LD Scores is too small or the window size used is too large. 

If this does happen to you, please look at the enrichment heatmap and check the
following:

- Is there a specific category that appears to be more often negatively 
enriched than not?
  - If so, then add the category to the list of blacklisted categories in
  `CreateAnnotationFile.R` (the vector can be found at the end of the file).
- Are the negative enrichment values somewhat evenly distributed throughout the
heatmap?
  - If so, then make the window size smaller in `7_ReferenceLDSCore.sh`. Change
  the flag `--ld-wind-cm 1` to have a smaller value than 1.
  - This line can be found in the reference LD scores section of the script

If problems continue to persist, there may be a problem with the dataset that
you are calculating LD scores for. Create an issue for this repository or
[ldsc](https://github.com/bulik/ldsc/issues) for further assistance.
