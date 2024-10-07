---
title: ContiguousStateLengths
description: "The script used to determine the lengths of contiguous state assignments."
sidebar_position: 5
---

# Contiguous State Lengths

This script attempts to answer a similar question to [./IsolationScores.md]. 
At a base level, this question is trying to answer:

> Is this state useful?

It does this in a similar way to [./IsolationScores.md] by assessing whether or
not a state is 'isolated'. However it covers a certain edge case that cannot be
handled by [./IsolationScores.md]. One big problem with computing Isolation
scores is that it doesn't handle the cases like:

-----x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x---------------

(where "x" is the state in question and "-" is any other state).

The above shows that state "x" is not very stable. There is never a contiguous
assignment of the state of more than one bin. However, The isolation score
of this state (if this was the only manner in which it was assigned to the
genome by ChromHMM) would be 1 (relatively low).

To cover for this edge case (and a few others that look very similar). This
script was implemented. All the script does is look for a state assignment and
then report how many states downstream of it have the same assignment (until
this is no longer the case). If a state is likely to have a very low value 
(below a user defined threshold) for this metric (over the whole genome), then 
the state is flagged as potentially redundant.

The metric used for 'most likely' is a difficult one here however. After 
looking at the distributions of contiguous state lengths it is obvious that
there is a significant positive skew in the distribution. Due to this skew,
the mean is out of the window. It is also quite clear that ChromHMM is
generally quite unstable as the mode of the distribution (across all states 
over lots of different models) is 1 (or very close to 1). As such, we decided
that using the median was the best choice (to get around these shortcomings).

## Extra thoughts

Now you could argue that such a state isn't necessarily redundant. Such a state
pattern as shown above can certainly be relevant. The reason we remove this
case is because it isn't very interpretable in a biological context. How can
you actually show that such a state exists across datasets and is functionally
relevant when its structure is so 'fragile'? Maybe in the future, when 
sequencing has advanced even further, it could be interpretable. But not right
now (2024).
