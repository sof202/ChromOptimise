---
title: FlankingStates
description: "The script used to calculate the most likely flanks."
sidebar_position: 3
---

# FlankingStates

## Explanation

The definition of 'flanking state' is as follows:
> For a bin that is assigned state $x$, the downstream flanking state is the state assignment for the next bin that is not also assigned state $x$. The upstream flanking state is the same but for the previous bin that is not also assigned $x$.

Consider the below example:

| Bin number | 1 | 2 | 3 | 4 | 5 | 6 | 7 |
| -- | -- | -- | -- | -- | -- | -- | -- |
| state assignment | 1 | 1 | 1 | 2 | 2 | 3 | 2
| | | | | | | | |

In this example, lets find the flanking states for bin number 4. The next bin (bin 5) is also assigned state 2, so we go to the next bin (and so on) until we reach a bin that is not assigned 2. In this case, we go to bin 6 and find the downstream flank is 3. In the same way we can go backwards to find the upstream flank is 1 (assignment for bin 3).

The most likely flanking states are the expected up/downstream flanking states for any given instance of a state in the state assignment file. These most likely flanking states are what this script aims to find.

The most likely flanking states are found using the transition matrix for the model being inspected. ChromHMM doesn't just factor in the transition probabilities when finding the most likely state assignment. However, using the probabilities from these transition matrices is more generalisable to any set of observations (and also saves massively on processing time). 
