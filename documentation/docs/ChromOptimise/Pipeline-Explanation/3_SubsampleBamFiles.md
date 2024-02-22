---
title: 3_SubsampleBamFiles
description: "The script used to merge and subsample processed files."
sidebar_position: 5
---

# 3_SubsampleBamFiles

:::warning[Skipping script 1]
If you skipped script 2 (2_batch_ProcessBamFiles.sh), please ensure that your files are separated into different folders (by epigenetic mark). The names of these folders doesn't matter, but we recommend naming them after the epignetic mark name (so help messages make sense).
:::

## Explanation

This merges all of the processed `.bam` files from the previous step using `samtools merge` and subsequently subsamples this file using `samtools view -s`.

The reason for merging files and subsampling instead of simply subsampling the original files is that it is unlikely that every processed `.bam` file is of the same size. If the `.bam` files are of different sizes, then different subsamples with the same number of files can have different sizes.
\
This results in the different samples no longer being directly comparable as they will likely have different sizes.

:::note[Large number of files]
If you are working with a lot of files, this step can use a large amount of memory due to the way samtools merges files, see [memory profiling] for more information.
:::

:::warning[outdated samtools]
Samtools had a known bug in version 0.1.18 where subsampling would fail if the sample size was above 50%, ensure that your version of samtools is up to date.
:::

# Example usage

```bash
# Merges all files in ${PROCESSED_DIR}/H3k27me3 and subsamples this file
# with sample size 75%
sbatch 3_SubsampleBamFiles.sh \
--config="path/to/configuration/directory" \
--mark="H3K27me3" \
--samplesize=75
```