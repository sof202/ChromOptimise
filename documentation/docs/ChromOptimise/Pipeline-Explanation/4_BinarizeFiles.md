---
title: 4_BinarizeFiles
description: "The script used to binarize the dataset for ChromHMM."
sidebar_position: 6
---


# 4_BinarizeFiles

## Explanation

This uses ChromHMM's `BinarizeBam` command to binarize the subsampled `.bam` file produced in the previous step.

It is assumed at this point that the at least step 3 has been completed for all of the epigenetic marks that you aim to use in the hidden Markov models, and that the sample size used for each mark is the same.

:::note[Unwanted chromosomes]
In step 7, the tool ldsc is used to inspect biological relevance in the chromatin
states obtained from ChromHMM. This step specifically will only consider autosomal
chromosomes. As a result, sex (and mitochondrial) binary files are deleted at the
end of this script. Comment these lines out if you want to keep these chromosomes
for the model generation step (step 5).
:::

### Bin size

This step takes a user input for the bin size, this will have a noticable effect on the quality of your models. ChromHMM recommends a bin size of 200 and this is reflected in this pipeline. It is important to consider the implications of changing the bin size. A lower value will increase the granularity of the data, but increases the susceptibility to random noise in the dataset. A larger value will achieve the opposite. 

To learn more about how the bin size (and other factors) affect the final result of the pipeline, head to [this page](/ChromOptimise/Factors-that-affect-the-output.md).

# Example usage

```bash
# Takes all files in ${SUBSAMPLE_DIR} that used a sample size of 75% and
# binarizes them using a bin size of 200 and the hg19 assembly
sbatch 4_BinarizeFiles.sh \
--config="path/to/configuration/directory" \
--binsize=200 \
--samplesize=75 \
--assembly=hg19
```
