---
title: Overall Explanation
description: "The purpose of the supplementary pipeline."
sidebar_position: 1
---

# Overall Explanation

There are 
[additional scripts](https://github.com/sof202/ChromOptimise/tree/main/supplementary) that
are given with ChromOptimise that are not directly a part of the main pipeline.
However, these scripts proved useful in constructing the pipeline and so were
elected to be included in this repository.

## Explanation

It became abundantly clear during the production of this pipeline that no
matter the method for determining the optimal number of states to use with a
model, a threshold had to be created. The concept of a threshold is very common
in statistics, one needs to draw the line somewhere for significance. One sees
these thresholds all the time in statistical analysis, perhaps in a p-value or
a confidence interval. However, it is less common that a compelling explanation
behind *why* a threshold is actually given. Thresholds for p-values are usually
chosen just because 'they are really small'. This isn't really good enough as
it provides an avenue for biases to enter the analysis (whether intentionally
or not). This obviously defeats the point in making a threshold in the first
place.

This supplementary pipeline was generated to provide a starting point for
defining such thresholds. The main idea behind the pipeline is: In a model with
an abundance of states (way too many for the given dataset), most states will
be redundant. Only a few states in an exceedingly large model will actually be
relevant. Hopefully, there will be a noticable difference between these types
of states ('obviously redundant' and 'obviously relevant') under the metrics 
that are used by the main pipeline. 

This pipeline will create these large models and showcases the metrics used in
the main pipeline for said models. 
