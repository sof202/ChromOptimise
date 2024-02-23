---
title: Execute_Redundancy_Metrics
description: "The script that produces the redundancy metrics for a single model file."
sidebar_position: 3
---

# Execute_Redundancy_Metrics

## Explanation

This script runs the same R scripts that 6_OptimalNumberOfStates.sh does. Head over to the [relevant pages] for explanations on the individual scripts. The only difference is that a histogram and heatmap is generated for the Euclidean distances calculated (as it can be rather difficult to interpret these files when there are lots of states).


## Obtaining 'good' threshold parameters

Regardless of how these thresholds are chosen, they will always have an element of subjectivity to them. A threshold is flawed in this way. But hopefully, the outputs of this script (that can be found in subdirectories where the model file inputted resides) can assist in choosing the values for these parameters more intelligently.

### `emissions_threshold`

The `emissions_threshold` is the more important of the two thresholds. The histogram especially should assist in deciding a suitable value for this. Ideally your histogram should look something like this:

![Euclidean distances histogram](/euclidean-distances/Euclidean_distances_histogram_model-80.png)

From this plot, there is a very obvious gap at around x=0.5. Going below this value and you see a much larger number of state pairs that have a 'low' Euclidean distance between them. The existence of the gap is important as it allows one to separate the two groups with a hard line. In the worst case scenario the histogram will be completely uniform and it will be impossible to tell where the cut off point should be placed. The value of 0.5 is of course, still just a suggestion. One could be more or less stringent with the thresholds, just make sure that the value chosen has some meaning behind it.

### `isolation_threshold`

The `isolation_threshold` is less important but still has a large impact on your dataset. This value is unlikely to change with the number of marks in the dataset and so it will only need to be changed when a new dataset or chromosome is being inspected. The isolation scores of the model analysed will be in the same directory that the model file is in. Opening the text file will reveal something like:

```text title="Isolation_Scores_model-[num].txt"
"state" "isolation_score"
1 0.0491809430882105
2 4.51568957599878
3 5.73507906818568
4 0.801777139523403
...
```

This file is not converted into a plot as the isolation scores can be very large for some states (if they were only assigned a few times across the chromosome inspected). The general strategy is to see any big jumps in the isolation score. Common states will have very low isolation scores, whilst rarer states will have much higher isolation scores (or none at all). In this light, the process is very similar to deciding on `emissions_threshold`.


:::tip[consistency]
Make sure to be consistent with your threshold parameters. The threshold parameters will likely be ill advised if in your actual analysis you: include more marks, use different bin/sample sizes, inspect different chromosomes, use different Phred scores (*etc.*). 
:::

## Example Usage

```bash
# Generates redundancy metrics  for the model with 50 states that used
# a random seed of 1 in its initialisation
# The output paramter can generally be left blank as the script defaults to the big models directory.
sbatch Execute_Redundancy_Metrics.sh
--config="path/to/configuration/directory" \
--size=50 \
--seed=1 \
--chromosome=1 \
--output=${BIG_MODELS_DIR}
```