---
title: IsolationScores
description: "The script used to calculate the isolation scores for each state."
sidebar_position: 4
---

# IsolationScores

## Motivation

In the pursuit of attempting to determine the states that are spatially
insignificant, it can be tempting to look at which states are exceedingly rare.
However, this doesn't capture the whole picture. Transcriptional start sites
are rare in the context of the human genome, but that doesn't make them
insignificant.

Because of this problem, a different metric was created in the isolation score.

## Explanation

The isolation score for a given state is defined as:
> The number of bins that separate two bins with the same state assignment on
> average

The isolation score is calculated for a model on a single user specified
chromosome. The reason why we don't consider all chromosomes in this
calculation is that it is difficult to come up with a concrete way of tackling
the following scenario:

Suppose a state is highly clustered on one chromosome (very low isolation
score) but is not assigned at all on another chromosome. 

We have determined that it is best in this scenario to consider the two
outcomes separately instead of trying to merge them into one single
observation. 

Another problem that comes with inspecting all chromosomes comes with the
smaller chromosomes. If a state is only assigned to two bins (that are
adjacent) on a smaller chromosome, then the isolation score will be 0. This
doesn't really capture the essence behind the isolation score as the states are
still relatively isolated (but there happens to be two next to each other).
Instead of increasing the complexity of the definition for isolation score, we
just calculate the scores for the largest chromosome (chromosome 1). 

The user has the power to change the chromosome used in this analysis via the
`--chromosome` flag used when calling
[3_OptimalNumberOfStates.sh](./3_OptimalNumberOfStates.md). This could prove
useful if a particular chromosome is of interest.
