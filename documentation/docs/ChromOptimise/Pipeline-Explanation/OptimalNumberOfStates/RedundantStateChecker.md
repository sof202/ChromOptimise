---
title: RedundantStateChecker
description: "The script used to determine the redundant states in a model."
sidebar_position: 5
---

# RedundantStateChecker

## Explanation

This script uses the information gathered from SimilarEmissions.R, FlankingStates.R and IsolationScores.R to determine the existence of any redundant states in the model being inspected. 

The criterion used to determine if a state is redundant is given as:

- There exists another state with similar emission parameters (under Euclidean distance)
- At least one of the following:
  - The state is highly isolated (under isolation score) 
  - The similar states are likely to be flanked by the same states
  - At least one of the similar states is flanked on both sides by the other state (in the similar state pair)

This script will summarise the states and state pairs that exceed any thresholds given in the configuration file in its output. If any states are determined redundant, the model is rejected in favour of in the next less complex model.

## Reasoning behind criterion

The basic non-mathematical definition of a redundant state trying to be captured by the criterion is:
> A redundant state is a state that is already captured by a more spatially significant state

### Emission paramters

The "already captured" part is being checked with the first point in the criterion. If two states have similar emission parameters, then they capture the same information (when the spatial nature of the model is disregarded). Lots of other previous approaches with finding the optimum number of state focus on this quality of models. However, simply having two states that have similar emission parameters isn't enough, the spatial context that comes with them also matters. This context is captured by the second part of the criteron

### High isolation

If a state is highly isolated (and is already similar to another state), it implies that the model is picking up some of the random noise left in the data. In the extreme case, suppose a state is only assigned once on a chromosome. This not only applies that it is exceedingly rare, but also that when it does appear, it's just a 'blip' in the data. Under this light, the state is not spatially significant.

Being highly isolated is not good enough on its own however. Transcriptional start sites are highly isolated functional elements of the genome, but that doesn't make them 'random noise', quite the opposite. 

### Flanking states

Instead of thinking about the genome, think about a hidden markov model modelling the weather (for ease of explanation). You *could* have two states that are identical in terms of emission parameters (two states that appear cloudy). However, these states could still have temporal (we are now talking about time, but spatial and temporal are interchangable in this explanation) significance if they give you information about what happened in the past and what is about to happen. Suppose one type of cloudy is likely to preceed rain and the other clear sunny weather. Now the two states suddenly bring a great deal of information with them.

On the flip side, if the two types of cloudy both preceeded and exceeded the same types of weather, then there is very little setting them apart from one another (one of them is redundant). Also if one of the types of cloudy is likely to exceed and preceed the other type is cloudy, then once again, one of the states is redundant (the model is just capturing random noise).