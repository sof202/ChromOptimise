---
sidebar_position: 4
---

# Supplementary pipeline explanation

There are [additional scripts](https://github.com/sof202/ChromOptimise/tree/main/supplementary) that are given with ChromOptimise that are not directly a part of the main pipeline. However, these scripts proved useful in constructing the pipeline and so were elected to be included in this repository.

## CompareModels.sh

This uses ChromHMM's `CompareModels` command to generate comparisons between the models produced in the model learning step of the main pipeline.

This script will compare each model with all of the models that are less complex than it (so a model with 8 states will be compared against only those models that have fewer than 8 states).

This script is useful in assessing which states are described well in less complex models and which are lost when reducing model complexity.

**Note**: `CompareModels` only looks at the emission parameters for each state in the models inputted. Therefore, this doesn't capture redundant states.

Example:

```shell
# Compares all models in ${MODEL_DIR} to one another
sbatch CompareModels.sh
```

## Redundandancy_Threshold_Optimisation

In [config.R](./Configuration-Files-Setup.md#configr), the user needs to decide upon the threshold parameters used in identifying redundant states in models. It is tempting to pick an arbitrary small number here. However, this small pipeline aids the user in choosing more informed values for these parameters.

For a visual representation of the pipeline, please consult this [schematic representation](/pipelines/Supplementary_Pipeline.pdf).

### Generate_Big_Model.sh

This uses the binarized data produced by [4_BinarizeBamFiles.sh](./Pipeline-Explanation.md#4_binarizefilessh) to produce a single very complex model with a high number of user specified states.
The model produced will use a random initialisation of emission/transition parameters so that an arbitrary number of states can be used in the model (the default 'information' initialisation method does not allow for this).

**Note**: ChromHMM uses the Baum-Welch algorithm to train the model. This algorithm finds a **local** maximum of the [estimated log likelihood](./ChromHMM-overview.md#estimated-log-likelihood), **not** the **global** maximum.
Therefore, because the initialisation is random, the models produced by this script can take very different lengths of time to train and may be highly innacurate (despite their complexity).

Example:

```shell
# Generates a model with 50 states using a random seed of 1 using assembly 
# hg19 for the binary data with bin size 200 and sample size 75%
sbatch Generate_Big_Model.sh 50 1 200 75 hg19
```

### Generate_Redundancy_Metrics_Plots.sh

This runs two R scripts on a specified model (recommended to be a very complex model) so that the user can make informed decisions on the thresholds used in determining redundant states in 6_OptimalNumberOfStates.sh.

HistogramPlotForEuclideanDistances.R -> Creates a histogram for the Euclidean distances between pairs of states in the model provided.

ScatterPlotForTransitionMaxima.R -> Creates a scatter plot for the maximum transition probability towards each state in the model provided.

The user should input a very complex model as this will result in the two plots possessing two obvious groups, redundant states and useful states. The gaps between the obvious groups will provide a the user with a good idea of where to draw the line between 'useful state' and 'redundant state'.

For Euclidean distances -> Redundant state pairs have very low Euclidean distance, whilst non-redundant state paris will have a noticeably higher Euclidean distance
\
For Transition probabilities -> Redundant state candidates score very low whilst useful state candidates score closer to 1 (with a noticeable gap between these groups).

This script is only optional as the user is free to set their own redundancy determining thresholds if they have a good idea for the values from other sources.

Example:

```shell
# Generates redundancy metrics plots for the model with 50 states that used
# a random seed of 1 in its initialisation
sbatch Generate_Redundancy_Metrics_Plots.sh 50 1
```
