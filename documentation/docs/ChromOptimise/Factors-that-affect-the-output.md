---
sidebar_position: 5
---

# Factors that affect the output

The optimum number of states to use with your dataset will likely change
depending on the quality of your data and the user inputs during the pipeline.
Below are some examples of this.
## Sample size
- A larger sample size will result in a higher level of completeness for your
  data, possibly increasing model accuracy
## Bin size
- This is used in the binarization process
- Smaller bin sizes will be more susceptible to random background noise in the 
data
- Larger bin sizes lose data precision
## Number of epigenetic marks
 - A large number of epigenetic marks is likely to require more states than a 
 small number of marks to achieve optimality
 ## Quality of data
- The random noise apparent in the raw sequencing data
    - This is directly linked to bin size
- The presence/absence of a control file
    - Scripts can be altered to allow for a control file to be used during the
      binarization process
    - Control files can vastly improve the algorithm's success rate in
      differentiating between background noise and actual peaks
