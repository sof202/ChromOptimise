---
sidebar_position: 2
---

# Configuration files setup

You will need to create three configuration files for this pipeline to work:

- [FilePaths.txt](#filepathstxt)
- [config.R](#configr)
- [LogFileManagement.sh](#logfilemanagementsh)

These files are used by each of the R and bash scripts to aid in organisation of the scripts and avoid repetition.

If you want to just run one script that does everything for you (no need to understand how each script works), you can create one more configuration file in [ChromOptimiseConfig.txt](#chromoptimiseconfigtxt). This file will be sourced in ChromOptimise.sh and all arguments will be parsed for you.


:::info
The scripts in this pipeline do not create the directory structure themselves. This is to avoid large files being dumped in unwanted locations. Please check your file paths are correct and then run Create_File_Structure.sh.
:::

:::warning[EOL errors]
You must ensure that these files are written with EOL: \n (LF) and not EOL: \r\n (CRLF).
:::

After one creates each of these configuration files, place them in the 'configuration' directory. Then run the `setup` executable from the main directory.

:::note[blueprint]
The pipeline was completed with blueprint data in mind, if your data is already downloaded, processed, binarized etc. then the associated lines in the config files might not be required.
:::
## Data directory structure

A guide for the structure of the data directory is given below (You only need to create the directories starting with an integer):

```text
Main_Data_Directory
├── 0_Downloaded_Files
├── 1_Organised_Raw_Bam_Files
│   ├── Epigenetic_Mark_1
│   ├── ...
│   └── Epigenetic_Mark_n
├── 2_Processed_Bam_Files
│   ├── Epigenetic_Mark_1
│   ├── ...
│   └── Epigenetic_Mark_n
├── 3_Subsampled_Bam_Files
├── 4_Binary_Files
│   ├── BinSize_xxx_SampleSize_yyy
│   ├── ...
│   └── BinSize_zzz_SampleSize_www
├── 5_Model_Files
├── 6_Optimum_Number_Of_States
│   ├── Results_From_Run_1
│   ├── Results_From_Run_2
│   ├── ...
│   ├── Results_From_Run_n
│   └── Likelihood_Values_Of_Models
├── 7_LDSC_Assessment_Files
│   ├── Results_From_Run_1
│   ├── Results_From_Run_2
│   ├── ...
│   └── Results_From_Run_n
├── 8_Model_Comparison_Files
└── 9_Big_Model_Files
    └── Plots
        ├── Euclidean_Distance_Histrograms
        └── Transition_Maxima_Scatter_Plots


LDSC_reference_files
├── PLINK_files
├── Weights
└── GWAS_traits
```

## FilePaths.txt

```text title="FilePaths.txt"
## Data directories

export MAIN_DIR="full/path/to/main/directory"
export DOWNLOAD_DIR="${MAIN_DIR}/0_Downloads"
export RAW_DIR="${MAIN_DIR}/1_RawBamFiles"
export PROCESSED_DIR="${MAIN_DIR}/2_ProcessedBamFiles"
export SUBSAMPLED_DIR="${MAIN_DIR}/3_SubsampledBamFiles"
export BINARY_DIR="${MAIN_DIR}/4_BinarizedBamFiles"
export MODEL_DIR="${MAIN_DIR}/5_ModelFiles"
export OPTIMUM_STATES_DIR="${MAIN_DIR}/6_OptimumNumberOfStates"
export LD_ASSESSMENT_DIR="${MAIN_DIR}/7_LDSCFiles"
export COMPARE_DIR="${MAIN_DIR}/8_ModelComparisonFiles"
export BIG_MODELS_DIR="${MAIN_DIR}/9_BigModelFiles"

## LDSC data directories

export LD_DIR="full/path/to/ldsc/directory"
export LD_PLINK_DIR="${LD_DIR}/PLINK_files"
export LD_WEIGHTS_DIR="${LD_DIR}/weights_files"
export LD_FRQ_DIR="${LD_DIR}/frq_files"
export LD_GWAS_TRAITS_DIR="${LD_DIR}/gwas_traits"
export LD_BASELINE_DIR="${LD_DIR}/baseline"

## LDSC reference file prefixes

export PLINK_PREFIX="prefix.for.plink.files"
export WEIGHTS_PREFIX="prefix.for.weights.files"
export FRQ_PREFIX="prefix.for.frq.files"
export BASELINE_PREFIX="prefix.for.baseline.files"

## Script directories

export SCRIPTS_DIR="full/path/to/this/repository"
export RSCRIPTS_DIR="${SCRIPTS_DIR}/Rscripts"
export LOG_DIR="${SCRIPTS_DIR}/LogFiles"

## Pyega3/conda directories

export CONDA_SHELL="path/to/conda/etc/folder"
export PYEGA_ENVIRONMENT="path/to/pyega3/conda/environment"
export LDSC_ENVIRONMENT="path/to/ldsc/conda/environment"
export CREDENTIALS="path/to/EGA/credentials"

## ChromHMM file locations

export CHROMHMM_MAIN_DIR="/path/to/ChromHMM/main/directory"
export CHROMHMM_CHROM_SIZES="${CHROMHMM_MAIN_DIR}/path/to/chromosome/sizes"

## LDSC software

export LD_SOFTWARE_DIR="path/to/ldsc/main/directory"
```

## config.R

To get a good value for the thresholds in the redundancy parameters section, please consult the [supplementary pipeline](/category/supplementary-pipeline---usage-and-explanation).

```R title="config.R"
## Redundancy parameters

emissions_threshold=VALUE
isolation_threshold=VALUE

## Number of marks used in analysis

number_of_marks=VALUE
```

## LogFileManagement.sh

This script will produce information on the time for the script to complete, reducing repetition of code
\
(Alternatively you can use the `time` command with the scripts, this is not available on our HPC).

```shell title="LogFileManagement.sh"
#!/bin/bash
## ============= ##
##   JOB START   ##
## ============= ##

echo "Job '${SLURM_JOB_NAME}' started at:"
date -u

start_time=$(date +%s)

LOG_FILE_PATH="${LOG_DIR}/$SLURM_JOB_NAME/$USER"
mkdir -p "${LOG_FILE_PATH}"
timestamp=$(date -u +%Y.%m.%d-%H_%M)
export timestamp

## ============= ##
##   FUNCTIONS   ##
## ============= ##

## ====== FUNCTION : finishing_statement() ===========================================
## Description: Give finishing message then exit
## Globals: 
##     start_time
## Locals:
##     end_time
##     time_taken
## Arguments:
##     exit code
## ===================================================================================
finishing_statement(){
    echo "Job finished with exit code $1 at:"
    date -u
    local end_time
    local time_taken
    end_time=$(date +%s)
    time_taken=$((end_time-start_time))
    echo "Job took a total of: ${time_taken} seconds to finish."
    exit "$1"
}
```


## ChromOptimiseConfig.txt

This is the configuration file that enables the user to run all of the files in the [main pipeline](/category/main-pipeline---usage-and-explanation) sequentially. Options are briefly described in comments here, but for a better picture of what to put here we recommend looking at the pipeline explanation. If you do not plan on using certain scripts at all (which is likely the case for 0_EGADownloading.sh for example) you can just remove the options section for those selected scripts.

```text title="ChromOptimiseConfig.txt"
# Which shell script to start from (provide a number from 0 to 6)
export STARTING_SCRIPT=

# This is the list of marks that you intend to use in the analysis
# Please provide this as a white space separated array
# i.e. (mark1 mark2 mark3 mark4 ...)
export LIST_OF_MARKS=

# This is a FULL file path to the FOFN which contains files you want to
# download using the pyega3 client
export FILE_OF_FILE_NAMES=path/to/file

# This is a threshold for the Phred score used in the processing stage
# (which reads to discard due to low base accuracy)
export PRED_SCORE_THRESHOLD=20

# This is the sample size (as a percentage) to use in the subsampling stage
# If your data is small in size, the recommended value is 100
export SAMPLE_SIZE=50

# This is the bin size to use during the binarization stage
# ChromHMM recommends a default of 200
export BIN_SIZE=200

# This is the assembly that your data is alligned to.
export ASSEMBLY=hg19

# This is the number of models you wish to create in the model learning stage
# Read the documentation on 5_batch_CreateIncrementalModel.sh for help
# here
export NUMBER_OF_MODELS=4

# This is the chromosome identifier (1-22, x,y,m) for the chromosome that is
# to be used with the isolation metric (see Pipeline-explanation). Unless
# this matters for your analysis, keep this at 1.
export CHROMOSOME_IDENTIFIER=1

# If you are running the script starting from script 7 (in which case this
# master script is a little overkill), you need to input the number of states
# you deem to be optimal. If starting on script <=6, this is not required
export OPTIMUM_NUMBER_OF_STATES=

# The following are the time you wish to allocate to each job (without this
# jobs have no maximum walltime and can run forever).
# Variables are numbered by their corresponding script number 
# (0 -> EGADownloading.sh)
# Times must be in the format: hh:mm:ss
export MAXTIME_0=96:00:00
export MAXTIME_1=00:30:00
export MAXTIME_2=12:00:00
export MAXTIME_3=12:00:00
export MAXTIME_4=04:00:00
export MAXTIME_5=12:00:00
export MAXTIME_6=00:20:00
export MAXTIME_7=00:10:00


# The following are the array sizes you wish to use for batch jobs (processing
# and model learning scripts). Using a number larger than the number of cores
# on your system is not recommended
export PROCESSING_ARRAY_SIZE=4
export MODEL_LEARNING_ARRAY_SIZE=4
```
