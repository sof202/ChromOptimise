---
sidebar_position: 2
---

# Configuration files setup

You will need to create three configuration files for this pipeline to work:

- FilePaths.txt
- config.R
- LogFileManagement.sh.

These files are used by each of the R and bash scripts to aid in organisation of the scripts and avoid repetition.

**IMPORTANT**: You must ensure that these directories actually exist on your system before running any of the scripts in the pipeline. The scripts do not create their own directories (except sub-directories) as this could lead to unwanted folder positions if the user incorrectly assigns the paths in the config files.

**IMPORTANT**: You must ensure that these files are written with EOL: \n (LF) and not EOL: \r\n (CRLF).

After one creates each of these configuration files, place them in the 'configuration' directory. Then run the `setup` executable from the main directory.

Note: The pipeline was completed with blueprint data in mind, if your data is already downloaded, processed, binarized etc. then the associated lines in the config files might not be required.

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
├── 7_Big_Model_Files
│   └── Plots
│       ├── Euclidean_Distance_Histrograms
│       └── Transition_Maxima_Scatter_Plots
└── 8_Model_Comparison_Files
```

## FilePaths.txt

```text
## Data directories

export MAIN_DIR="path/to/main/directory"
export DOWNLOAD_DIR="${MAIN_DIR}/path/to/downloads"
export RAW_DIR="${MAIN_DIR}/path/to/raw/data"
export PROCESSED_DIR="${MAIN_DIR}/path/to/processed/data"
export SUBSAMPLED_DIR="${MAIN_DIR}/path/to/subsampled/data"
export BINARY_DIR="${MAIN_DIR}/path/to/binary/data"
export MODEL_DIR="${MAIN_DIR}/path/to/chromHMM/models"
export OPTIMUM_STATES_DIR="${MAIN_DIR}/path/to/optimum/states/output"
export COMPARE_DIR="${MAIN_DIR}/path/to/comparison/files"
export BIG_MODELS_DIR="${MAIN_DIR}/path/to/big/models"

## Script directories

export SCRIPTS_DIR="/path/to/this/repository"
export RSCRIPTS_DIR="${SCRIPTS_DIR}/Rscripts"
export SUPPLEMENTARY_DIR="${SCRIPTS_DIR}/supplementary"
export LOG_DIR="${SCRIPTS_DIR}/LogFiles"

## ChromHMM file locations

export CHROMHMM_MAIN_DIR="/path/to/ChromHMM/main/directory"
export CHROMHMM_CHROM_SIZES="${CHROMHMM_MAIN_DIR}/path/to/chromosome/sizes"
```

## config.R

To get a good value for the thresholds in the redundancy parameters section, please consult the [supplementary pipeline](./Supplementary-pipeline-explanation.md).

```R
## Data Directories

main_dir="path/to/main/directory"
model_dir=paste0(main_dir, "path/to/model/files")
optimum_states_dir=paste0(main_dir, "path/to/optimum/states/output")
likelihood_dir=paste0(optimum_states_dir, "Likelihood_Values")
compare_dir=paste0(main_dir, "path/to/comparison/files")
big_models_dir=paste0(main_dir,"path/to/big/model/files")

## Plotting directories

transition_plotting_dir=paste0(big_models_dir,"/path/to/plots")
emission_plotting_dir=paste0(big_models_dir,"/path/to/plots")

## Redundancy parameters

emissions_threshold=VALUE
transitions_threshold=VALUE
isolation_threshold=VALUE

## Number of marks used in analysis

number_of_marks=VALUE
```

## LogFileManagement.sh

This file is essential for managing the temporary log files that are produced when sending the jobs through SLURM workload manager. On top of this, the script will produce information on the time for the script to complete.
\
(Alternatively you can use the `time` command with the scripts, this is not available on our HPC).

```shell
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
## Description: Delete temporary log and error files, give finishing message then exit
## Globals: 
##     SLURM_SUBMIT_DIR
##     SLURM_JOB_ID
##     start_time
## Locals:
##     end_time
##     time_taken
## Arguments:
##     exit code
## ===================================================================================
finishing_statement(){
    rm "${SLURM_SUBMIT_DIR}/temp${SLURM_JOB_ID}.log" 
    rm "${SLURM_SUBMIT_DIR}/temp${SLURM_JOB_ID}.err"
    echo "Job finished with exit code $1 at:"
    date -u
    local end_time
    local time_taken
    end_time=$(date +%s)
    time_taken=$((end_time-start_time))
    echo "Job took a total of: ${time_taken} seconds to finish."
    exit "$1"
}

## ====== FUNCTION : batch_finishing_statement() ======================================
## Description: Delete temporary log and error files, give finishing message then exit
## Globals: 
##     SLURM_SUBMIT_DIR
##     SLURM_ARRAY_JOB_ID
##     SLURM_ARRAY_TASK_ID
##     start_time
## Locals:
##     end_time
##     time_taken
## Arguments:
##     exit code
## ===================================================================================
batch_finishing_statement(){
    rm "${SLURM_SUBMIT_DIR}/temp${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log" 
    rm "${SLURM_SUBMIT_DIR}/temp${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err"
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
