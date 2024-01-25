---
sidebar_position: 6
---

# Processing times

SLURM workload manager has the feature of introducing a maximum wall time into each script. The default maximum wall time is set in the `#SBATCH` section of each script and can be overridden when setting the `--time` option when using `sbatch`. If the time of execution for a script exceeds this maximum wall time, it will terminate early. This is useful as it helps with workload management on HPCs that a large number of users have access to, and also reduces the effects of non-halting programs. 
\
Using bash builtins and time utilities, we were able to create estimates on how long each script is expected to take depending on the inputs and complexity of information. 
## Contents
- [General CPU information which testing was conducted on](#general-cpu-information-which-testing-was-conducted-on)
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

## General CPU information which testing was conducted on
Note that processing times are likely to vary across systems, the scripts were tested on a high performance cluster (HPC) that is split into several compute nodes. Each script was ran on a compute node which have the following properties:
- CPUs -> 16
- Sockets -> 2
- Cores per socket -> 8
- Threads per core -> 1
- CPU clock speed -> 2.601 GHz (max)
- CPU model -> Intel(R) Xeon(R) CPU E5-2640 v3 
- Cache
    - L1i (information) -> 32K
    - L1d (data) -> 32K
    - L2 -> 256K
    - L3 -> 20 MB

## 0_DownloadBluePrint.sh
The time taken to download files using pyega3 will of course depend on internet speed, size of files and number of files being downloaded. The former of these three makes it difficult to create well informed estimates for the run time of this script. In testing, the internet speed did not appear throttled and was consistent across file downloads.

## 1_MovingFilesToSingleDirectory.sh
Moving files on Linux systems is generally very fast even when files are large or abundant. The main souce of processing time in this script is the `find` command. Depending on how many files exist in the main directory you are working with, the processing time for `find` can vary as more files need to be checked. 
\
Regardless, in previous tests moving 130 files (each being 1-2 GB) amongst ~1000 total files took roughly one minute. It is unlikely that the processing time will ever exceed the default maximum wall time set in the `#SBATCH` parameters at the top of this script.

## 2_batch_ProcessBamFiles.sh
The processing time for this script mainly depends on the sizes of the files that are being processed. In tests thus far, the main contributor to computational time comes from `samtools sort` and appears to have a linear relationship with file size.
\
In testing, the following was observed:
-  Processing a ~1.1 GB file took ~10 minutes
-  Processing a ~1.4 GB file took ~12 minutes
-  Processing a ~1.8 GB file took ~16 minutes
-  Processing a ~2.1 GB file took ~19 minutes

The conclusion was made that one should expect that the processing time would be roughly:
\
[9 * file size (in GB)] minutes.
\
Additionally, it's important to note that this script is designed to be executed as an array through the SLURM workload manager. Therefore, the processing time is likely to vary depending on the number of cores assigned to each task in the array. Further note that if the size of the array is larger than the number of files being processed, all files will be processed by the highest indexed array element (causing slow down).
\
When applying an array of size 4, the following was observed:
-  Processing ~272 GB of files took ~ 11 hours and 45 minutes
-  Processing ~228 GB of files took ~ 9 hours and 53 minutes
-  Processing ~263 GB of files took ~ 10 hours and 47 minutes

This to some extent supports the above claim for the processing time being linear. It is unlikely that using an array size of 4 will result in exactly 4 times faster processing time. Therefore, it is difficult to support the processing time equation previously given using these tests.  

## 3_SubsampleBamFiles.sh
This script's largest contributor to computational time is `samtools merge`. Depending on the number of files and size of said files, the time taken can vary dramatically.
\
In testing, the following was observed:

Using a 50% sampling rate:
| Number of files | Average file size (GB) | time (minutes) |
| --              | --                     | --             |
| 2               | 3.5                    | 4
| 2               | 6.5                    | 7
| 133             | 8.3                    | 906
| 171             | 6.2                    | 1010
| 176             | 7.6                    | 1271

Using a 100% sampling rate:
| Number of files | Average file size (GB) | time (minutes) |
| --              | --                     | --             |
| 133             | 8.3                    | 909
| 171             | 6.2                    | 1186
| 176             | 7.6                    | 1366

Currently this shows a linear relationship with total file size. There is not enough information to determine the effect that the sampling rate has on the processing time (though it is likely to be minimal).

## 4_BinarizeBamFiles.sh
The processing time for this script is of course mainly taken up by ChromHMM's `BinarizeBam` command. This command will generally take up more time if more data is inputted and if a smaller bin size is used.
\
In testing, the following was observed:
- Using 17 GB of .bam files:
    - Using a bin size of 200bp took ~4 minutes
- Using 3730 GB of .bam files:
    - Using a bin size of 200bp took 11 hours and 52 minutes
    - Using a bin size of 400bp took 11 hours and 47 minutes

## 5_batch_CreateIncrementalModels.sh
This script's largest contributor to computational time is ChromHMM's `LearnModel` command. This command will take up more time if the bin size was chosen to be smaller in 4_BinarizeBamFiles.sh (leading to a larger total size of binary signal files) and if the number of states to be used in the model increases.
\
From the [user manual](https://compbio.mit.edu/ChromHMM/ChromHMM_manual.pdf) and the [source code](https://github.com/jernst98/ChromHMM/blob/master/edu/mit/compbio/ChromHMM/ChromHMM.java), it is clear that the model training is completed via the forward-backwards algorithm and the Baum-Welch algorithm. These algorithms are standard in hidden Markov model training and likely take up the majority of processing time.
\
From a [resource](https://www.cs.ubc.ca/~murphyk/Bayes/rabiner.pdf) sent by one of the creators, Jason Ernst, these algorithms have a computational complexity of O(N<sup>2</sup>T). Where N is the number of states in the model and T is the total number of observations (proportional to sizes of binary files).

**Note**: Although the algorithms mentioned have computational complexity of O(N<sup>2</sup>T), this doesn't serve as the complete picture. By default, ChromHMM will terminate once the change in the likelihood function is below 0.001. The number of iterations it takes to reach this point can vary massively between models and datasets. The maximum number of iterations (by default) is 200. Keep this maximum iteration number into account when considering processing times.
\
In testing, the following was observed:

For 1.6 MB of binary data (400bp bins, 3 marks)
| Number of states | Average time per iteration (seconds) |
| - | - |
| 2 | 1.27 |
| 3 | 1.31 |
| 4 | 1.71 |
| 5 | 2.17 | 
| 6 | 2.66 |
| 7 | 3.49 |
| 8 | 4.76 |


It is important to note that this script is designed to be executed as an array through the SLURM workload manager. Therefore, the processing time is likely to vary depending on the number of cores assigned to each task in the array. Further note that if the size of the array is larger than the number of models to learn, all models will be learnt by the highest indexed array element (causing significant slow down).

## 6_OptimalNumberOfStates.sh
This script makes use of an R script that is inside of a loop. However, the Rscript has very little computational complexity and so the script is usually very fast. In previous tests with 8 state models, the script took ~14 seconds. It is unlikely that this script will ever exceed the default wall time set in the `#SBATCH` parameters at the top of the script.

## Generate_Big_Model.sh 
Computational time for this script follows the same logic as [5_batch_CreateIncrementalModels.sh](#5_batch_createincrementalmodelssh). The main contributor to computational time is ChromHMM's `LearnModel` command, only this time the models are likely to be much larger. As a result expect the processing time to be very long if a high number of states is to be used.
\
In testing, the following was observed:
- For ~1 MB of binary data (200bp bins)
    - 20 states took 1 hour 48 minutes
    - 30 states took 3 hours 48 minutes
    - 40 states took 7 hours 9 minutes
    - 60 states took 14 hours 36 minutes
    - 70 states took 19 hours 18 minutes
    - 80 states took 25 hours 34 minutes

## Generate_Redundancy_Metrics_Plots.sh
This script calls two R scripts that make one plot each. The amount of data that is parsed into them is very small (unless the models have hundreds of states). Therefore, despite R being relatively slow, this script should never take any longer than one minute. In testing, processing a model with 80 states took 12 seconds.

## CompareModels.sh
Most computational time is accredited to ChromHMM's `CompareModels` command. However, this command does not require a large amount of time to run. In tests, comparing models to a model of 8 states (then to 7 states, 6, 5 etc.) took only 18 seconds. Therefore it is unlikely that this script will ever exceed the default wall time set in the `#SBATCH` parameters at the top of the script. 
