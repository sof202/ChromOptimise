---
title: 0_EGADownloading
description: "The script used to download bam files from EGA"
sidebar_position: 2
---

# 0_EGADownloading

## Explanation

This is an artefact of the pipeline originally being built for the blueprint data from EGA.

The script will download files from a file of file names using the pyega3 python package. It is likely that this script is not relevant to you, and you can safely ignore it.


## Prerequisites

### Step one: Login credentials need to be saved in `.json` format:

These should be provided to you via EGA.

```json
{
    "username":"username@domain.com",
    "password":"password"
}
```

### Step Two: Create a conda environment with pyega3 installed:

Note down the location of this conda environment. It is to be used in the configuration file [FilePaths.txt].

```bash
# Change 'myenv' to a memorable name
conda create -n myenv pyega3
```

## Example usage

```bash
# Downloads all files listed in FileOfFileNames.txt from EGA
sbatch 0_EGADownloading.sh \
--config="path/to/configuration/directory" \
--file="FileOfFileNames.txt"
```

