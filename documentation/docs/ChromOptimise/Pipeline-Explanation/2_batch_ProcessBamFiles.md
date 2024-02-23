---
title: 2_batch_ProcessBamFiles
description: "The script used to process collated bam files."
sidebar_position: 4
---

# 2_batch_ProcessBamFiles

:::warning[Skipping script 1]
If you skipped script 1 (1_MoveFilesToSingleDirectory.sh), please ensure that your files are separated into different folders (by epigenetic mark). The names of these folders doesn't matter, but we recommend naming them after the epignetic mark name (so help messages make sense).
:::

## Explanation

:::info[blueprint specifics]
The processing steps taken in this script are highly specific to the ChIP-Seq data obtained from the blueprint consortium. In the event that your data is not similar to this dataset, we recommend that the user completes this step on their own accord and skips this script.
:::

The purpose of this step is to sort, index, filter and remove the duplicates of the aligned `.bam` files. This step also outputs some statistics on the aligned read files using `samtools stats` and `samtools idxstats` before and after processing (for additional user information). 

Due to the current specificity of this script, it will likely only work with data that satisfies the following:

- Aligned read files are single end (if paired end, you could discard one of the strands)
- Picard (or similar) has been used to mark duplicate reads in the data

## Example usage

```bash
# Processes all files in ${RAW_DIR}/H3K27me3 using a Phred score threshold
# of 30 in an array of size 4 (jobs split into 4 subprocesses)
sbatch --array=1-4 2_batch_ProcessBamFiles.sh \
--config="path/to/configuration/directory" \
--mark="H3K27me3" \
--phred=30
```