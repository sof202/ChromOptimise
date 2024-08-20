---
title: Recreate_Transition_Heatmap
description: "A script to recreate a transition parameter heatmap, but better."
sidebar_position: 4
---

# Recreate_Transition_Heatmap

This script will create a transition heatmap from a ChromHMM transitions data
file but removes the main diagonal. The main diagonal's values are usually all
close to 1. This means you can't make out the smaller differences between the
rest of the transition parameters (due to the colouring method for each box).
Removing the main diagonal allows the user to gain more information at a
glance. This is because the more interesting state transition probabilities
(not between non-identical states) are now more vibrant and contrasted with one
another.


## Example usage

```bash
sbatch Recreate_Transition_Heatmap.sh \
    path/to/configuration/directory \
    path/to/transition_x.txt \ # Exists alongside other ChromHMM files 
    path/to/output/directory
```

Alternatively, you can just run the Rscript yourself, the accompanying bash
script is there mainly for consistency:

```bash
Rscript \
    RecreateTransitionHeatmap.R \
    path/to/transition_x.txt \
    path/to/output/directory 
```
