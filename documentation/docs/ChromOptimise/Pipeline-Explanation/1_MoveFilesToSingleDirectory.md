---
title: 1_MoveFilesToSingleDirectory
description: "The script used to collate and organise files."
sidebar_position: 3
---

# 1_MoveFilesToSingleDirectory

:::warning[skipping this script]
If this step is skipped, please ensure that the raw `.bam` files are organised into folders with memorable names (epigenetic mark name is reccomended).

Otherwise subsequent scripts will fail.
:::

## Explanation

If you have very few files then you can likely achieve the output of this script manually. This script is in place as the number of files you are working with could be in the hundreds if the dataset is large enough.

If you downloaded files from EGA, each file will be in a separate directory (as each directory has its own checksum file to validate the downloads). It can be tedious to move each of these files yourself.

:::warning[Important]
This script relies on the `.bam` files containing the epigenetic mark in the file names. Please ensure this is the case for your raw files. 

If this is not the case with your dataset: Instead of inputting the mark name, input a different string into the commandline that you are using to differentiate the files.
:::

## Example usage

```bash
# Moves all files containing 'H3K27me3' to: ${RAW_DIR}/H3k27me3
sbatch 1_MoveFilesToSingleDirectory.sh \
--config="path/to/configuration/directory" \
--mark="H3K27me3"
```
