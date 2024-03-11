---
title: HeritabilityPlots
description: "The script used to create plots of partitioned heritability."
sidebar_position: 3
---

# HeritabilityPlots

## Explanation

This script will take in the `.results` files that are outputted from ldsc's partitioned heritability and visualise them for you. The two output files are:

1) A heatmap of the enrichment (Proportion of heritability for state / proportion of SNPS for state) for each state against each GWAS trait selected
2) A barplot for the p-values of the enrichment scores for each state for each GWAS trait selected 

This is the main output of runLDSC.sh that enables the user to interpret the biological relevance of their states. If lots of states are showing very little enrichment for any of the traits that you are considering, this suggests that the model is too specific still (or there is a problem with your input data).