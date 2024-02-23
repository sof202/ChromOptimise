---
title: "Overall explanation of pipeline"
description: "Read this first to gain more context"
sidebar_position: 1
---

# Overall explanation of pipeline

The main pipeline of ChromOptimise is split into 3 main steps:

1) Manipulate the dataset using user defined parameters
2) Use ChromHMM's pipeline to obtain multiple hidden Markov models on the processed data
3) Use statistical analysis to evaluate the model with the optimum number of states

The pipeline expects aligned `.bam` files as its input. One may want to download these files from EGA, this can be done using [0_EGADownloading.sh](./0_EGADownloading.md).

For a visual representation of the pipeline, please consult this [schematic representation](/pipelines/Optimal_States_Analysis_Pipeline.pdf).

After completion, alongside the optimal number of states to use with the dataset, ChromOptimise will output the following:

- A file listing the determined redundant states for each model that was rejected
- A plot (and csv file) of the bayesian information criterion for each model learned
  - For help with interpretting this plot, head to [this page]
- A plot of the estimated log likelihood for each model learned

Alongside the above, the following information is outputted for each rejected model:

- The [isolation scores](./OptimalNumberOfStates/IsolationScores) for each state in the model
- The [most likely flanking states](./OptimalNumberOfStates/FlankingStates) for each state in the model
- The [Euclidean distances](./OptimalNumberOfStates/SimilarEmissions) between each pair of states in the model
- The collated information from the above points in a more digestible format

For an example output of the pipeline, see [these pages].