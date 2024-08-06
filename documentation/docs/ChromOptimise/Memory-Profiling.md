---
sidebar_position: 7
---

# Memory profiling

To gain insights into the peak heap consumption of each script in the pipeline,
we utilized the tool Heaptrack. This allowed for the estimation of memory
consumption based on unit tests (albeit only with relatively small amounts of
data). Utilizing this information will empower the user to make well-informed
decisions regarding the allocation of memory for the scripts involved in the
pipeline and ChromHMM.

Heaptrack, produced by KDE, can be obtained through an appimage download or
local compilation. The GitHub repository for Heaptrack can be accessed
[here](https://github.com/KDE/heaptrack).

## 1_BinarizeBamFiles.sh
ChromHMM's `BinarizeBam` command is the main contributor to memory consumption
in this script.

In testing, the following was observed:
- For a directory with 16.8 GB of .bam files: peak heap memory consumption was
  40 MB

## 2_batch_CreateIncrementalModels.sh
ChromHMM's `LearnModel` command is the main contributor to memory consumption 
in this script.

In testing, the following was observed:
- Using binary files that were 1 MB in total:
    - 2 state model training had 181 MB peak memory consumption
    - 3 state model training had 145 MB peak memory consumption
    - 4 state model training had 165 MB peak memory consumption
    - 5 state model training had 128 MB peak memory consumption
    - 6 state model training had 165 MB peak memory consumption
    - 7 state model training had 168 MB peak memory consumption

Based on the provided information, it is evident that the memory consumption is
primarily dependent on the total size of the binary files rather than the
number of states used in the models.

Additionally, it's important to note that this script is designed to be
executed as an array through the SLURM workload manager. Therefore, the peak
memory consumption is anticipated to be proportional to the size of the largest
.bam files multiplied by the number of array elements.

## 3_OptimalNumberofStates.sh
This script utilizes Rscripts which entail minor computations, resulting in
very low memory consumption. Through testing, it has been observed that the
memory usage never surpassed a few kilobytes.

## 4_ReferenceLDSCore.sh
This script's biggest requirement of memory comes from the SNP assigment R
script that it calls. In testing, chromosomes 1 and 2 use about 800MB of
memory. Because the script is ran as an array that processes each chromosome at
the same time, 22 GB of memory is allocated.

## 5_PartitionedHeritability.sh
According to [ldsc's creators](https://github.com/bulik/ldsc/wiki/FAQ),
partitioned heritability should only take up ~8GB for 50 annotations. The number
of annotations in this script will be greater than 50. In testing, over 20GB is
required to stop out of memory errors for just 58 annotations. To account for
this, the SLURM directive is set to 50GB of memory.

## Generate_Big_Model.sh 
This script is similar to
[3_batch_CreateIncrementalModels.sh](#3_batch_createincrementalmodelssh),
suggesting that the memory consumption pattern is likely to align with the
conclusions drawn there. However, it's important to note that due to the
computational time required for this script, utilizing heaptrack to monitor
memory allocations was impractical. Consequently, the memory consumption for
notably large models (20 states or more) remains unknown and is not extensively
documented.

