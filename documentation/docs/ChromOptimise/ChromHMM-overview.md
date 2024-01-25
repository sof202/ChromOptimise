---
sidebar_position: 1
---

# ChromHMM overview

[ChromHMM](https://compbio.mit.edu/ChromHMM/#:~:text=ChromHMM%20is%20software%20for%20learning,and%20spatial%20patterns%20of%20marks.) is a tool created by Jason Ernst and Manolis Kellis. It uses the Baum-Welch algorithm to train a hidden Markov model to help characterise chromatin states for a genomic dataset (specifically one that can be peak called). This page covers what a hidden Markov model actually is and how the model is trained.

## Hidden Markov models
Hidden Markov models are a subset of Markov models. A common example used to describe how such models work involve weather prediction. Suppose there are three types of weather: Sunny, cloudy and rainy. Markov models could be applied to predict what the weather will be like on each day of a given week.

In a Markov model, one knows what the weather is like today and uses this information to make a prediction for the weather tommorow. It is assumed that there is a fixed rate of transitioning between each state (weather type) for each day, and the weather for tommorow only depends on today's weather ('memoryless property'). The model consists of transition probabilities between each weather state only. For example: if it is rainy today, it is unlikely to be sunny tommorow; if it is sunny today, it is likely to be sunny again tommorow; if it is cloudy today, it is equally likely be rainy, sunny or cloudy again tommorow (*etc.*). A useful property of Markov models is that they only rely on the previous state, therefore (in this example) the weather for next week (or next month *etc.*) can be predicted using today's weather only (albeit poorly).

In a *hidden* Markov model, one no longer knows what the weather is like today. Instead they recieve indirect observations that are assumed to be dependent on the state of the weather. These observations could be anything from the temperature to the clothes their friends are wearing that day. These models are called hidden markov models as the states (weather) is no longer directly observable. These models still use the transition probabilities between  states seen in normal markov models, but now have the added feature of emission probabilities. These emission probabilites are the probability of some obseravation *given* a certain state. In the weather example: Sunny days are likely to be hot, people are likely to wear coats when it rains (*etc.*). A combination of the observations and the previously assumed state are used to make predictions on the following state.

Hidden markov models therefore are very useful in the context of chromatin analysis. Instead of weather states, we now have specific chromatin states (such as genes, intronic regions, transcription end sites *etc.*). Instead of observations of temperature, clothing *etc.* we instead have combinations of epigenetic marks. Hidden Markov models can therefore make predictions on the genome that would be difficult to make otherwise.

## Estimated log likelihood 
In order to train the parameters for a hidden Markov model, ChromHMM uses the Baum-Welch algorithm. The aim of the Baum-Welch algorithm is to maximise the 'likelihood function'. The 'likelihood function' is a metric used to identify how accurate a model is when comparing it to the dataset. 

Informally, the likelihood function is 'the probability that this sequence of observations would be observed under the hidden Markov model'.

Formally, the likelihood funtion: $P(O | λ)$ is $P(O_1,O_2,\dots,O_T | λ)$ where:
- $O_i$ is the observation at time/position $i$
- $T$ is the total number of observations
- $λ$ is the parameters of the hidden Markov model:
    - The transition probabilities between states
    - The emission probabilities of each state
    - The initial state probabilities

The likelihood function is calculated using the [forwards-backwards algorithm](https://en.wikipedia.org/wiki/Forward–backward_algorithm) which is outlined in full in [this document](https://www.cs.ubc.ca/~murphyk/Bayes/rabiner.pdf).