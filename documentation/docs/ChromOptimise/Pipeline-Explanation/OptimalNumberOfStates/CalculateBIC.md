---
title: CalculateBIC
description: "The script used to plot the relative BICs of the models."
sidebar_position: 7
---

# CalculateBIC

## Explanaiton

The Bayesian information critereon (BIC) is a heuristic for measuring the
"information" captured by a model. In contrast with the Akaike information
criterion (AIC), BIC is much more punishing on the number of parameters in the
model when the dataset is very large (which it will be in the case of genomic
datasets). In general AIC is easier to interpret. However, with large datasets,
increasing the number of states is almost always improves the value of this
heuristic (as the estimated log likelihood is so low).

## Definition
The Bayesian information critereon is given by: 

$$
BIC = ln(n)k - 2ln(L)
$$

Where $n$ is the total number of observations, $k$ is the number of parameters
in the model and $L$ is the estimated likelihood.

## Relative BIC

These BIC values are genrally very large and hard to compare. As a result the
script will find the BIC of each model relative to the minimum BIC of all
models (this time using a simple ratio as there is no merit to using $exp$ with
this heuristic).


A smaller value for BIC generally indicates a better model. It is important to
note that BIC is a heuristics, not a metric. It merely approximates how good a
model is relative to its complexity. As such, the plot and .csv files that come
from this script is only given to provide further information to the user. It
isn't directly implemented into the pipeline as it isn't concrete.

You might use BIC in the following scenario: Perhaps a model with 6 states is
better under this metric than a model with 8 states, despite the more complex
model having no 'redundant states'. Considering the model with 6 states has a
lower BIC than the model with 8, the user may wish to save on computing power
and use the 6 state model in further downstream analysis.
