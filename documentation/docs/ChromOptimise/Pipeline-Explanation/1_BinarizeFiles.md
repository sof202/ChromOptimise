---
title: 1_BinarizeFiles
description: "The script used to binarize the dataset for ChromHMM."
sidebar_position: 2
---


# 1_BinarizeFiles

## Explanation

This uses ChromHMM's `BinarizeBam` and `BinarizeBed` commands to binarize the 
`.bam` and `.bed` files that are your inputs.

It is assumed that your bam/bed files are named as:

mark-name-1.bed
mark-name-2.bam
mark-name-3.bed

It is also assumed that you have one file per mark. If this is not the case,
consider using `bedtools` or `samtools` to merge your `bed`/`bam` files into one.

:::note[Unwanted chromosomes]
In step 4, the tool ldsc is used to inspect biological relevance in the chromatin
states obtained from ChromHMM. This step specifically will only consider autosomal
chromosomes. As a result, sex (and mitochondrial) binary files are deleted at the
end of this script. Comment these lines out if you want to keep these chromosomes
for the model generation step (step 2).
:::

## Bin size

This step takes a user input for the bin size, this will have a noticable 
effect on the quality of your models. ChromHMM recommends a bin size of 200 
and this is reflected in this pipeline. It is important to consider the 
implications of changing the bin size. A lower value will increase the 
granularity of the data, but increases the susceptibility to random noise in 
the dataset. A larger value will achieve the opposite. 

To learn more about how the bin size (and other factors) affect the final 
result of the pipeline, head to
[this page](/ChromOptimise/Factors-that-affect-the-output.md).


## ChromBinarize

ChromHMM's inbuilt binarization tools are fairly primitive and limited. You may
instead want to use a tool specifically for peak calling (like MACS) or have
data that cannot be peak called easily. For this reason, we have also developed
[ChromBinarize](https://github.com/sof202/ChromBinarize). This is a tool that
can convert ONT, BS-Seq and MACS peak called datasets into ChromHMM binary 
format.
