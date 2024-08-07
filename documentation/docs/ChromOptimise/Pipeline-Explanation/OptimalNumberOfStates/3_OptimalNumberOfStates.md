---
title: 3_OptimalNumberOfStates
description: "The script used to determine the optimum number of states."
sidebar_position: 1
---

# 3_OptimalNumberOfStates

## Explanation

This is the wrapper script that will parse your arguments into mutliple R 
scripts in order to determine redundant states. It starts with the most 
complex model and works its way down until a model with no redundant states is 
found. This model is classed as being optimal.

:::note[Max model number]
There is a possible source of error in the above logic where the output will 
not actually be the optimal number of states to use. Suppose the maximum 
number of states out of any model is only 3 (but the number of marks is much 
larger than this), then it is likely that this model has no redundant states. 
This doesn't actually suggests that the optimal number of states is 3, just 
that it is 3 or greater. The additional information of a model with 4 states 
having redundant states is what confirms the model with 3 states is optimal.
:::

## Thresholds

The major factor that affects the output of this script lies in the thresholds
defined in the configuration file for R Config.R. It can be difficult to
determine suitable values for these thresholds without prior knowledge. To help
with this, the supplementary pipeline was created. Running the supplementary
pipeline can help the user in identifying 'good' values for these thresholds
based off of the information gained from much larger models (that undoubtedly
have redundant states). To learn more about this pipeline, consult 
[these pages](/category/supplementary-pipeline---usage-and-explanation).

### The curse of dimensionality

The threshold for isolation scores is likely to be similar across datasets, but
the value for `emissions_threshold` will not. This threshold determines when
two states are 'too similar' in terms of their emission parameters (according
to the Euclidean distance metric).

Different datasets with differing numbers of marks will certainly require
different values for this threshold. The reason for this is as follows: The
number of marks is the dimension of the emission parameter vectors (the number
of columns in the emission files). As you add more dimensions to a vector
space, the distance between two random vectors is also expected to increase. It
becomes increasingly difficult for two random vectors to be 'very close' to one
another.

Consider a line of length 1 meter. Pick two random points on the line and at a
maximum they can be 1 meter from one another. The expected distance between the
points will be 0.5 meters (hopefully this is obvious). Now suppose you have a
square with side lengths 1. When you place two points inside this square at
random, then the expected distance between them is now $\frac{2}{3}$ meters. In
any one direction, the max distance is still 1 meter, but now there is a new
dimension that the points can differ in, meaning the distance between them can
be even larger than in the 1 dimensional case. This continues as you think
about cubes and hyper cubes (*etc.*).

As such, no value is given for the `emission_threshold`, this value should be
changed depending on the number of marks you have included in your dataset.

