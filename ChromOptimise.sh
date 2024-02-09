#!/bin/bash

## ===========================================================================##
##                                                                            ||
##                                     PREAMBLE                               ||
##                                                                            ||
## ===========================================================================##
## PURPOSE:                                                                   ||
## This is a master script that will run all of the scripts in JobSubmission  ||
## sequentially using a configuration file inputted by the user.              ||
## ===========================================================================##
## AUTHOR: Sam Fletcher                                                       ||
## CONTACT: s.o.fletcher@exeter.ac.uk                                         ||
## CREATED: February 2023                                                     ||
## ===========================================================================##
## PREREQUISITES: Run setup                                                   ||
## ===========================================================================##
## DEPENDENCIES:                                                              ||
## SLURM Workload Manager                                                     ||
## Java, ChromHMM                                                             ||
## R                                                                          ||
## Samtools                                                                   ||
## Python, pyega3, Miniconda/Conda/Anaconda                                   ||
## ===========================================================================##
## INPUTS:                                                                    ||
## $1 -> Full (or relative) file path for configuation file                   ||
## ===========================================================================##
## OUTPUTS:                                                                   ||
## The optimum number of states to use with your data (given user parameters) ||
## ===========================================================================##

## ========= ##
##   SETUP   ##
## ========= ##

# To avoid excessive source commands and specifying full file paths
# we need to change the directory to ChromOptimise/Jobsubmission (so scripts
# can be called using just their file name)

CHROMOPTIMISE_DIR=$(dirname "$0")
JOBSUBMISSION_DIR="${CHROMOPTIMISE_DIR}/JobSubmission"

cd "${JOBSUBMISSION_DIR}" || { echo "An error occured, the JobSubmission 
directory: ${JOBSUBMISSION_DIR} does not exist (please check this)."; exit 1; }

# We need an associative array to store the job ids for each SLURM job
# that are submitted. This is so that we can use these job ids as dependencies.
# Dependencies allow us to queue up job submissions and run certain scripts
# in parallel whenever required

declare -A jobID

## =========================== ##
##   LOAD CONFIGURATION FILE   ##
## =========================== ##

# We just require the configuration directory here instead of the full
# file path as all of the configuration files should be in the same directory,
# reducing repetition
configuration_directory=$1

source "${configuration_directory}/ChromOptimiseConfig.txt" || \
{ echo "The configuration file does not exist in the specified location: \
${configuration_directory}/ChromOptimiseConfig.txt"; exit 1; }

echo "Configuration file used with ChromOptimise master script"
echo ""
cat "${configuration_directory}/ChromOptimiseConfig.txt"
echo ""

## ===================== ##
##   RUN CHROMOPTIMISE   ##
## ===================== ##

# The pipeline can be started from different scripts depending on where the user
# specifies. Sometimes the downloading stage or processing stage (etc.) is not
# required.

## ---0_EGADownloading.sh---------------------------------------------------- ##
if [[ "${STARTING_SCRIPT:=10}" -eq 0 ]]; then
    jobID[Download]=$( \
    sbatch \
    "0_EGADownloading.sh" \
    "${configuration_directory}" \
    "${FILE_OF_FILE_NAMES}" | \
    awk '{print $4}' \
    )

    echo "Submitted 0_EGADownloading.sh under job ID:"
    echo "${jobID[Download]}"
fi
## -------------------------------------------------------------------------- ##


## ---1_MoveFilesToSingleDirectory.sh---------------------------------------- ##
if [[ "${STARTING_SCRIPT}" -eq 1 ]]; then
    # This script can be ran in parallel for each mark as the scripts only 
    # move files corresponding to one mark each
    for mark in "${LIST_OF_MARKS[@]}"; do
        array_index_move="${mark}_move"

        jobID[$array_index_move]=$( \
        sbatch \
        "1_MoveFilesToSingleDirectory.sh" \
        "${configuration_directory}" \
        "${mark}" | \
        awk '{print $4}' \
        )

        echo "Submitted 1_MoveFilesToSingleDirectory.sh for $mark under job ID:"
        echo "${jobID[$array_index_move]}"
    done
   
elif [[ "${STARTING_SCRIPT}" -lt 1 ]]; then
    for mark in "${LIST_OF_MARKS[@]}"; do
        array_index_move="${mark}_move"
        
        # This script depends on the downloading script being finished
        jobID[$array_index_move]=$( \
        sbatch \
        --dependency=afterok:"${jobID[0]}" \
        "1_MoveFilesToSingleDirectory.sh" \
        "${configuration_directory}" \
        "${mark}" | \
        awk '{print $4}' \
        )

        echo "Submitted 1_MoveFilesToSingleDirectory.sh for $mark under job ID:"
        echo "${jobID[$array_index_move]}"
    done
fi
## -------------------------------------------------------------------------- ##


## ---2_batch_ProcessBamFiles.sh--------------------------------------------- ##
if [[ "${STARTING_SCRIPT}" -eq 2 ]]; then
    # This script can be ran in parallel for each mark as the scripts only 
    # process files corresponding to one mark each 
    for mark in "${LIST_OF_MARKS[@]}"; do
        array_index_process="${mark}_process"

        jobID[$array_index_process]=$( \
        sbatch \
        "2_batch_ProcessBamFiles.sh" \
        "${configuration_directory}" \
        "${mark}" \
        "${PRED_SCORE_THRESHOLD}" | \
        awk '{print $4}' \
        )

        echo "Submitted 2_batch_ProcessBamFiles.sh for $mark under job ID:"
        echo "${jobID[$array_index_process]}"
    done
   
elif [[ "${STARTING_SCRIPT}" -lt 2 ]]; then
    for mark in "${LIST_OF_MARKS[@]}"; do
        array_index_move="${mark}_move"
        array_index_process="${mark}_process"
        
        # This script depends on .bam files (for specific mark) being in 
        # the correct place in the file structure
        jobID[$array_index_process]=$( \
        sbatch \
        --dependency=afterok:"${jobID[${array_index_move}]}" \
        "2_batch_ProcessBamFiles.sh" \
        "${configuration_directory}" \
        "${mark}" \
        "${PRED_SCORE_THRESHOLD}" | \
        awk '{print $4}' \
        )

        echo "Submitted 2_batch_ProcessBamFiles.sh for $mark under job ID:"
        echo "${jobID[$array_index_process]}"
    done
fi
## -------------------------------------------------------------------------- ##


## ---3_SubsampleBamFiles.sh------------------------------------------------- ##
if [[ "${STARTING_SCRIPT}" -eq 3 ]]; then
    # This script can be ran in parallel for each mark as the scripts only 
    # merge files corresponding to one mark each
    for mark in "${LIST_OF_MARKS[@]}"; do
        array_index_merge="${mark}_merge"

        jobID[$array_index_merge]=$( \
        sbatch \
        "3_SubsampleBamFiles.sh" \
        "${configuration_directory}" \
        "${mark}" \
        "${SAMPLE_SIZE}" | \
        awk '{print $4}' \
        )
        
        echo "Submitted 3_SubsampleBamFiles.sh for $mark under job ID:"
        echo "${jobID[$array_index_merge]}"
    done
   
elif [[ "${STARTING_SCRIPT}" -lt 3 ]]; then
    for mark in "${LIST_OF_MARKS[@]}"; do
        array_index_process="${mark}_process"
        array_index_merge="${mark}_merge"

        # This script depends on .bam files (for specific mark) being already 
        # processed (and with a specific file name)
        jobID[$array_index_merge]=$( \
        sbatch \
        --dependency=afterok:"${jobID[${array_index_process}]}" \
        "3_SubsampleBamFiles.sh" \
        "${configuration_directory}" \
        "${mark}" \
        "${SAMPLE_SIZE}" | \
        awk '{print $4}' \
        )

        echo "Submitted 3_SubsampleBamFiles.sh for $mark under job ID:"
        echo "${jobID[$array_index_merge]}"
    done
fi
## -------------------------------------------------------------------------- ##


## ---4_BinarizeBamFiles.sh-------------------------------------------------- ##
if [[ "${STARTING_SCRIPT}" -eq 4 ]]; then
    jobID[binarization]=$(
    sbatch \
    "4_BinarizeBamFiles.sh" \
    "${configuration_directory}" \
    "${BIN_SIZE}" \
    "${SAMPLE_SIZE}" \
    "${ASSEMBLY}" \
    )

    echo "Submitted 4_BinarizeBamFiles.sh under job ID:"
    echo "${jobID[binarization]}"

# If the starting point is before this script then:

# This script cannot be ran until all marks have been subsampled.
# However, we do not know ahead of time how many completed jobs to look for
# since the number of marks being processed/subsampled is data-determined.
# Hence we use a checkpoint script that only terminates once the subsampled
# directory has the correct number of files in it (for the sample size given).
elif [[ "${STARTING_SCRIPT}" -lt 4 ]]; then
    # Checkpoint script should 
    jobID[checkpoint]=$(
    sbatch \
    --dependency=afterok:"${jobID[${array_index_process}]}" \
    "$CHROMOPTIMISE_DIR/ChromOptimiseCheckpoints/Subsampling_Checkpoint.sh" \
    "${configuration_directory}" \
    "${SAMPLE_SIZE}" \
    "${LIST_OF_MARKS[@]}" | \
    awk '{print $4}' \
    )

    echo "Submitted Subsampling_Checkpoint.sh under job ID:"
    echo "${jobID[checkpoint]}"

    # ------------------------------------------------------------- #
    
    jobID[binarization]=$(
    sbatch \
    --dependency=afterok:"${jobID[checkpoint]}" \
    "4_BinarizeBamFiles.sh" \
    "${configuration_directory}" \
    "${BIN_SIZE}" \
    "${SAMPLE_SIZE}" \
    "${ASSEMBLY}" \
    )

    echo "Submitted 4_BinarizeBamFiles.sh under job ID:"
    echo "${jobID[binarization]}"
fi
## -------------------------------------------------------------------------- ##


## ---5_batch_CreateIncrementalModels.sh------------------------------------- ##
if [[ "${STARTING_SCRIPT}" -eq 5 ]]; then
    jobID[Model_Learning]=$( \
    sbatch \
    "5_batch_CreateIncrementalModels.sh" \
    "${configuration_directory}" \
    "${NUMBER_OF_MODELS}" \
    "${STATE_INCREMENT}" \
    "${BIN_SIZE}" \
    "${SAMPLE_SIZE}" \
    "${ASSEMBLY}" \
    )

    echo "Submitted 5_batch_CreateIncrementalModels.sh under job ID:"
    echo "${jobID[Model_Learning]}"

elif [[ "${STARTING_SCRIPT}" -lt 5 ]]; then
    jobID[Model_Learning]=$( \
    sbatch \
    --dependency=afterok:"${jobID[binarization]}" \
    "5_batch_CreateIncrementalModels.sh" \
    "${configuration_directory}" \
    "${NUMBER_OF_MODELS}" \
    "${STATE_INCREMENT}" \
    "${BIN_SIZE}" \
    "${SAMPLE_SIZE}" \
    "${ASSEMBLY}" \
    )

    echo "Submitted 5_batch_CreateIncrementalModels.sh under job ID:"
    echo "${jobID[Model_Learning]}"
fi
## -------------------------------------------------------------------------- ##


## ---6_OptimalNumberOfStates.sh--------------------------------------------- ##
if [[ "${STARTING_SCRIPT}" -eq 6 ]]; then
    jobID[Optimal_States]=$( \
    sbatch \
    "6_OptimalNumberOfStates.sh" \
    "${configuration_directory}" \
    )

    echo "Submitted 6_OptimalNumberOfStates.sh under job ID:"
    echo "${jobID[Optimal_States]}"

elif [[ "${STARTING_SCRIPT}" -lt 6 ]]; then
    jobID[Optimal_States]=$( \
    sbatch \
    --dependency=afterok:"${jobID[Model_Learning]}" \
    "6_OptimalNumberOfStates.sh" \
    "${configuration_directory}" \
    )

    echo "Submitted 6_OptimalNumberOfStates.sh under job ID:"
    echo "${jobID[Optimal_States]}"
fi
## -------------------------------------------------------------------------- ##