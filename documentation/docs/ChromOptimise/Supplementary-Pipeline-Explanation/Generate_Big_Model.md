---
title: Generate_Big_Model
description: "The script that produces a single large model file."
sidebar_position: 2
---

# Generate_Big_Model

## Explanation

This uses the binarized data produced by 4_BinarizeBamFiles.sh to produce a
single very complex model with a high number of user specified states. The
model produced will use a random initialisation of emission/transition
parameters so that an arbitrary number of states can be used in the model (the
default 'information' initialisation method does not allow for this). 

Ideally the model size used will be large enough that the model produced is
'guaranteed' to have lots of redundant states. We appreciate that the model
learning time for this script can be incredibly long, but this stage is
important for determining good values for thresholds.

:::note[Baum-Welch] 
ChromHMM uses the Baum-Welch algorithm to train the model.
This algorithm finds a **local** maximum of the 
[estimated log likelihood](/ChromOptimise/ChromHMM-overview.md#estimated-log-likelihood),
**not** the **global** maximum.

Therefore, because the initialisation is random, the models produced by this
script can take variable lengths of time to train (for different seeds) and may
be highly innacurate (despite their complexity).
:::

## Example usage

```bash
# Generates a model with 50 states.
sbatch Generate_Big_Model.sh \
path/to/configuration/directory \
50
```
