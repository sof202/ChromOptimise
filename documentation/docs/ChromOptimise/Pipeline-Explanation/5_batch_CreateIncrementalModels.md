---
title: 5_batch_CreateIncrementalModels
description: "The script used to create/learn models using ChromHMM."
sidebar_position: 7
---

# 5_batch_CreateIncrementalModels

## Explanation

This uses ChromHMM's `LearnModel` command to learn multiple hidden Markov models. It also obtains the [estimated log likelihood](/ChromOptimise/ChromHMM-overview.md#estimated-log-likelihood) value of each model.

This script will always generate a model with 2 states. Inspecting the emission parameters for this simple model is a good way of validating your data.
\
Generally, most of your genomic data will be non-coding regions. The 2 state model should reflect this property in the overlap files produced.

### Maximum number of states

This script uses the 'information' initialisation method for the starting parameter set for the model. Due to this, the number of states in the model cannot exceed the total number of combinations of marks in your dataset.
\
One may not know this number before running the script. However, the maximum number of mark combinations one can have for a binary file will always be less than or equal to $2^k$, where $k$ is the total number of marks. If one exceeds the maximum number of states permitted for this initialisation technique, ChromHMM will tell the user the maximum number of states allowed in the associated error message.



# Example usage

```bash
# Generates 6 models using the hg19 assembly (sizes 2-7) using 
# the binarized files with bin size 200 and sample size 75%
# --array=1-3 is used so that each array element learns 2 models (optimal).
sbatch --array=1-3 5_batch_CreateIncementalModels.sh \
--config="path/to/configuration/directory" \
--nummodels=6 \
--increment=1 \
--binsize=200 \
--samplesize=75 \
--assembly=hg19
```

:::warning[consistency]
Ensure that the bin size, sample size and assembly are consistent across the pipeline.
:::