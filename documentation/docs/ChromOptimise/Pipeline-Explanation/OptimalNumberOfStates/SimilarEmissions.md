---
title: SimilarEmissions
description: "The script used to calculate Euclidean distances."
sidebar_position: 2
---

# SimilarEmissions

## Explanation

This script is rather simple, it takes each pair of states from the emission file parsed and calculates the Euclidean distance between them. The script has an optional flag that causes the script to create a heatmap and a histogram of all of the Euclidean distances calculated. These plots can be useful in identifying a suitable threshold parameter to use in the configuration file [config.R](/ChromOptimise/Configuration-Files-Setup.md#configr) when used on exceptionally large models in the supplementary pipeline.
