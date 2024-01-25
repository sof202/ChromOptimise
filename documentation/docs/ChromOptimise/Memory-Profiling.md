---
sidebar_position: 7
---

# Memory profiling

To gain insights into the peak heap consumption of each script in the pipeline, we utilized the tool Heaptrack. This allowed for the estimation of memory consumption based on unit tests (albeit only with relatively small amounts of data). Utilizing this information will empower the user to make well-informed decisions regarding the allocation of memory for the scripts involved in the pipeline and ChromHMM.
\
\
Heaptrack, produced by KDE, can be obtained through an appimage download or local compilation. The GitHub repository for Heaptrack can be accessed [here](https://github.com/KDE/heaptrack).


## Contents
- [0_DownloadBluePrint.sh](#0_downloadblueprintsh)
- [1_MovingFilesToSingleDirectory.sh](#1_movingfilestosingledirectorysh)
- [2_batch_ProcessBamFiles.sh](#2_batch_processbamfilessh)
- [3_SubsampleBamFiles.sh](#3_subsamplebamfilessh)
- [4_BinarizeBamFiles.sh](#4_binarizebamfilessh)
- [5_batch_CreateIncrementalModels.sh](#5_batch_createincrementalmodelssh)
- [6_OptimalNumberOfStates.sh](#6_optimalnumberofstatessh)
- [Generate_Big_Model.sh](#generate_big_modelsh)
- [Generate_Redundancy_Metrics_Plots.sh](#generate_redundancy_metrics_plotssh)
- [CompareModels.sh](#comparemodelssh)

## 0_DownloadBluePrint.sh
The memory consumption for this script depends on the abundance and size of the files that are being downloaded from EGA repositories. 

## 1_MoveFilesToSingleDirectory.sh
The memory usage for this script is typically minimal, thanks to the way Linux manages file movements without transferring data in memory, but rather updating the pointer in the file management system. Consequently, even when dealing with large multi-gigabyte files, the peak memory consumption remains in the kilobyte range. Hence, this script requires very little memory.
\
\
However, if there is an attempt to move files between different drives, the above characteristics no longer hold true, and a much higher peak memory consumption is expected due to the need for redundancy plans to address potential corruption.

## 2_batch_ProcessBamFiles.sh
The primary source of memory consumption in this script is the execution of the `samtools sort` command. This operation requires significantly more memory than any other part of the script due to the large number of temporary allocations it performs.
\
In testing, the following was observed:
- For a 1.5 GB compressed .bam file:
    - `samtools index` had 22.2 MB peak memory consumption
    - `samtools idxstats` had 2.6 MB peak memory consumption
    - `samtools stats` had 5 MB peak memory consumption
    - `samtools view` had 4.5 MB peak memory consumption
    - `samtools sort` had 927 MB peak memory consumption

Based on this information, it can be inferred that the peak memory consumption will correspond to the size of the largest .bam file processed by the script, with the peak consumption expected to be slightly less than the size of the largest file.
\
Additionally, it's important to note that this script is designed to be executed as an array through the SLURM workload manager. Therefore, the peak memory consumption is anticipated to be proportional to the size of the largest .bam files multiplied by the number of array elements.

## 3_SubsampleBamFiles.sh
This script is designed to utilize minimal memory when merging and subsampling a small number of files. During the merging process, `samtools merge` only opens a small component of each file being merged at a time, instead of opening each file in its entirety. The subsampling is carried out using `samtools view`, which uses very little peak memory consumption relative to file size (see [2_batch_ProcessBamFiles.sh](#2_batch_processbamfilessh)).
\
However, challenges may arise when merging a significantly larger number of files. According to the [documentation](http://www.htslib.org/doc/samtools-merge.html) for `samtools merge`, the command utilizes ~1 MB of data from each file at a time. It's important to be mindful of this information when merging large numbers of files within the pipeline. In the event of memory-related errors with this script, it is recommended to execute `ulimit -n` in the terminal to check the limit set on how many files can be open at one time. If the number of files to merge exceeds the output of `ulimit -n`, it is advisable to consider modifying the script to merge the files in batches instead of all at once.

## 4_BinarizeBamFiles.sh
ChromHMM's `BinarizeBam` command is the main contributor to memory consumption in this script.
\
In testing, the following was observed:
- For a directory with 16.8 GB of .bam files: peak heap memory consumption was 40 MB

## 5_batch_CreateIncrementalModels.sh
ChromHMM's `LearnModel` command is the main contributor to memory consumption in this script.
\
In testing, the following was observed:
- Using binary files that were 1 MB in total:
    - 2 state model training had 181 MB peak memory consumption
    - 3 state model training had 145 MB peak memory consumption
    - 4 state model training had 165 MB peak memory consumption
    - 5 state model training had 128 MB peak memory consumption
    - 6 state model training had 165 MB peak memory consumption
    - 7 state model training had 168 MB peak memory consumption

Based on the provided information, it is evident that the memory consumption is primarily dependent on the total size of the binary files rather than the number of states used in the models.
/
Additionally, it's important to note that this script is designed to be executed as an array through the SLURM workload manager. Therefore, the peak memory consumption is anticipated to be proportional to the size of the largest .bam files multiplied by the number of array elements.

## 6_OptimalNumberofStates.sh
This script utilizes Rscripts which entail minor computations, resulting in very low memory consumption. Through testing, it has been observed that the memory usage never surpassed a few kilobytes.

## Generate_Big_Model.sh 
This script is similar to [5_batch_CreateIncrementalModels.sh](#5_batch_createincrementalmodelssh), suggesting that the memory consumption pattern is likely to align with the conclusions drawn there. However, it's important to note that due to the computational time required for this script, utilizing heaptrack to monitor memory allocations was impractical. Consequently, the memory consumption for notably large models (20 states or more) remains unknown and is not extensively documented.

## CompareModels.sh
Based on testing, this script incurs minimal memory consumption. The majority of the script is dedicated to fundamental file management tasks and the utilization of ChromHMM's `CompareModels` command. Notably, the latter operation only necessitates kilobytes of peak heap consumption when comparing 7 different emission files.