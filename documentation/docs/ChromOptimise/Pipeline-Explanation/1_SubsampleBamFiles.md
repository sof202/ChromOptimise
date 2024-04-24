---
title: 1_SubsampleBamFiles
description: "The script used to merge and subsample processed files."
sidebar_position: 5
---

# 1_SubsampleBamFiles

## Explanation

This merges all of the processed `.bam` files from your set of processed files
using `samtools merge` and subsequently subsamples this file using 
`samtools view -s`. Please ensure that your files can be found in 
PROCESSED_DIR (whatever you set this as in the config file) in mark denoted
folders.

The reason for merging files and subsampling instead of simply subsampling the 
original files is that it is unlikely that every processed `.bam` file is of 
the same size. If the `.bam` files are of different sizes, then different 
subsamples with the same number of files can have different sizes.
\
This results in the different samples no longer being directly comparable as 
they will likely have different sizes.

:::note[Large number of files]
If you are working with a lot of files, this step can use a large amount of 
memory due to the way samtools merges files, see [memory 
profiling](/ChromOptimise/Memory-Profiling.md) for more information.
:::

:::warning[outdated samtools]
Samtools had a known bug in version 0.1.18 where subsampling would fail if the 
sample size was above 50%, ensure that your version of samtools is up to date.
:::

## Example usage

```bash
# Merges all files in ${PROCESSED_DIR}/H3k27me3 and subsamples this file
# with sample size 75%
sbatch 1_SubsampleBamFiles.sh \
--config="path/to/configuration/directory" \
--mark="H3K27me3" \
--samplesize=75
```
