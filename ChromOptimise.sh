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

# dirname only gets the relative path, we want the full path for robustness.
cd "${CHROMOPTIMISE_DIR}" || exit 1
CHROMOPTIMISE_DIR=$(pwd)

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

## ---1_SubsampleBamFiles.sh------------------------------------------------- ##
if [[ "${STARTING_SCRIPT}" -eq 1 ]]; then
    # This script can be ran in parallel for each mark as the scripts only 
    # merge files corresponding to one mark each
    for mark in "${LIST_OF_MARKS[@]}"; do
        array_index_merge="${mark}_merge"

        jobID[$array_index_merge]=$( \
        sbatch \
        --time="${MAXTIME_1}" \
        "1_SubsampleBamFiles.sh" \
        --config="${configuration_directory}" \
        --mark="${mark}" \
        --samplesize="${SAMPLE_SIZE}" | \
        awk '{print $4}' \
        )
        
        echo "Submitted 1_SubsampleBamFiles.sh for $mark under job ID:"
        echo "${jobID[$array_index_merge]}"
    done
   
elif [[ "${STARTING_SCRIPT}" -lt 1 ]]; then
    for mark in "${LIST_OF_MARKS[@]}"; do
        array_index_process="${mark}_process"
        array_index_merge="${mark}_merge"

        # This script depends on .bam files (for specific mark) being already 
        # processed (and with a specific file name)
        jobID[$array_index_merge]=$( \
        sbatch \
        --time="${MAXTIME_1}" \
        --dependency=afterok:"${jobID[${array_index_process}]}" \
        "1_SubsampleBamFiles.sh" \
        --config="${configuration_directory}" \
        --mark="${mark}" \
        --samplesize="${SAMPLE_SIZE}" | \
        awk '{print $4}' \
        )

        echo "Submitted 1_SubsampleBamFiles.sh for $mark under job ID:"
        echo "${jobID[$array_index_merge]}"
    done
fi
## -------------------------------------------------------------------------- ##


## ---2_BinarizeBamFiles.sh-------------------------------------------------- ##
if [[ "${STARTING_SCRIPT}" -eq 2 ]]; then
    jobID[binarization]=$( \
    sbatch \
    --time="${MAXTIME_2}" \
    "2_BinarizeBamFiles.sh" \
    --config="${configuration_directory}" \
    --binsize="${BIN_SIZE}" \
    --samplesize="${SAMPLE_SIZE}" \
    --assembly="${ASSEMBLY}" | \
    awk '{print $4}' \
    )

    echo "Submitted 2_BinarizeBamFiles.sh under job ID:"
    echo "${jobID[binarization]}"

# If the starting point is before this script then:

# This script cannot be ran until all marks have been subsampled.
# However, we do not know ahead of time how many completed jobs to look for
# since the number of marks being processed/subsampled is data-determined.
# Hence we use a checkpoint script that only terminates once the subsampled
# directory has the correct number of files in it (for the sample size given).
elif [[ "${STARTING_SCRIPT}" -lt 2 ]]; then
    # Checkpoint script should not finish searching until all subsampling has
    # finished. We use the downloading MAXTIME as this is likely to be very
    # high and wait until at least the final mark has been successfully merged
    # to start the search
    array_index_merge="${LIST_OF_MARKS[-1]}_merge"

    jobID[checkpoint]=$( \
    sbatch \
    --time="${MAXTIME_1}" \
    --dependency=afterok:"${jobID[${array_index_merge}]}" \
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
    --time="${MAXTIME_2}" \
    --dependency=afterok:"${jobID[checkpoint]}" \
    "2_BinarizeBamFiles.sh" \
    --config="${configuration_directory}" \
    --binsize="${BIN_SIZE}" \
    --samplesize="${SAMPLE_SIZE}" \
    --assembly="${ASSEMBLY}" | \
    awk '{print $4}' \
    )

    echo "Submitted 2_BinarizeBamFiles.sh under job ID:"
    echo "${jobID[binarization]}"
fi
## -------------------------------------------------------------------------- ##


## ---3_batch_CreateIncrementalModels.sh------------------------------------- ##
if [[ "${STARTING_SCRIPT}" -eq 3 ]]; then
    jobID[Model_Learning]=$( \
    sbatch \
    --time="${MAXTIME_3}" \
    --array=1-"${MODEL_LEARNING_ARRAY_SIZE}" \
    "3_batch_CreateIncrementalModels.sh" \
    --config="${configuration_directory}" \
    --nummodels="${NUMBER_OF_MODELS}" \
    --binsize="${BIN_SIZE}" \
    --samplesize="${SAMPLE_SIZE}" \
    --assembly="${ASSEMBLY}" | \
    awk '{print $4}' \
    )

    echo "Submitted 3_batch_CreateIncrementalModels.sh under job ID:"
    echo "${jobID[Model_Learning]}"

elif [[ "${STARTING_SCRIPT}" -lt 3 ]]; then
    jobID[Model_Learning]=$( \
    sbatch \
    --time="${MAXTIME_3}" \
    --array=1-"${MODEL_LEARNING_ARRAY_SIZE}" \
    --dependency=afterok:"${jobID[binarization]}" \
    "3_batch_CreateIncrementalModels.sh" \
    --config="${configuration_directory}" \
    --nummodels="${NUMBER_OF_MODELS}" \
    --binsize="${BIN_SIZE}" \
    --samplesize="${SAMPLE_SIZE}" \
    --assembly="${ASSEMBLY}" | \
    awk '{print $4}' \
    )

    echo "Submitted 3_batch_CreateIncrementalModels.sh under job ID:"
    echo "${jobID[Model_Learning]}"
fi
## -------------------------------------------------------------------------- ##


## ---4_OptimalNumberOfStates.sh--------------------------------------------- ##
if [[ "${STARTING_SCRIPT}" -eq 4 ]]; then
    jobID[Optimal_States]=$( \
    sbatch \
    --time="${MAXTIME_4}" \
    "4_OptimalNumberOfStates.sh" \
    --config="${configuration_directory}" \
    --chromosome="${CHROMOSOME_IDENTIFIER}" \
    --binsize="${BIN_SIZE}" \
    --samplesize="${SAMPLE_SIZE}" \
    --nummodels="${NUMBER_OF_MODELS}" | \
    awk '{print $4}' \
    )

    echo "Submitted 4_OptimalNumberOfStates.sh under job ID:"
    echo "${jobID[Optimal_States]}"

elif [[ "${STARTING_SCRIPT}" -lt 4 ]]; then
    jobID[Optimal_States]=$( \
    sbatch \
    --time="${MAXTIME_4}" \
    --dependency=afterok:"${jobID[Model_Learning]}" \
    "4_OptimalNumberOfStates.sh" \
    --config="${configuration_directory}" \
    --chromosome="${CHROMOSOME_IDENTIFIER}" \
    --binsize="${BIN_SIZE}" \
    --samplesize="${SAMPLE_SIZE}" \
    --nummodels="${NUMBER_OF_MODELS}" | \
    awk '{print $4}' \
    )

    echo "Submitted 4_OptimalNumberOfStates.sh under job ID:"
    echo "${jobID[Optimal_States]}"
fi
## -------------------------------------------------------------------------- ##

## ---5_ReferenceLDSCore.sh-------------------------------------------------- ##
if [[ "${STARTING_SCRIPT}" -eq 5 ]]; then
    jobID[LDSC]=$( \
    sbatch \
    --time="${MAXTIME_5}" \
    "5_ReferenceLDSCore.sh" \
    --config="${configuration_directory}" \
    --state="${OPTIMUM_NUMBER_OF_STATES}" \
    --binsize="${BIN_SIZE}" \
    --samplesize="${SAMPLE_SIZE}" \
    --nummodels="${NUMBER_OF_MODELS}" | \
    awk '{print $4}' \
    )

    echo "Submitted 5_ReferenceLDSCore.sh under job ID:"
    echo "${jobID[LDSC]}"

elif [[ "${STARTING_SCRIPT}" -lt 5 ]]; then
    jobID[LDSC]=$( \
    sbatch \
    --time="${MAXTIME_5}" \
    --dependency=afterok:"${jobID[Optimal_States]}" \
    "5_ReferenceLDSCore.sh" \
    --config="${configuration_directory}" \
    --binsize="${BIN_SIZE}" \
    --samplesize="${SAMPLE_SIZE}" \
    --nummodels="${NUMBER_OF_MODELS}" | \
    awk '{print $4}' \
    )

    echo "Submitted 5_ReferenceLDSCore.sh under job ID:"
    echo "${jobID[LDSC]}"
fi
## -------------------------------------------------------------------------- ##
