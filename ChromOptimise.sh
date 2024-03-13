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

# The pipeline can be started from different scripts depending on where the user
# specifies. Sometimes the downloading stage or processing stage (etc.) is not
# required.

## ---0_EGADownloading.sh---------------------------------------------------- ##
if [[ "${STARTING_SCRIPT}" -eq 0 ]]; then
    jobID[Download]=$( \
    sbatch \
    --time="${MAXTIME_0}" \
    "0_EGADownloading.sh" \
    --config="${configuration_directory}" \
    --file="${FILE_OF_FILE_NAMES}" | \
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
        --time="${MAXTIME_1}" \
        "1_MoveFilesToSingleDirectory.sh" \
        --config="${configuration_directory}" \
        --mark="${mark}" | \
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
        --time="${MAXTIME_1}" \
        --dependency=afterok:"${jobID[0]}" \
        "1_MoveFilesToSingleDirectory.sh" \
        --config="${configuration_directory}" \
        --mark="${mark}" | \
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
        --time="${MAXTIME_2}" \
        --array=1-"${PROCESSING_ARRAY_SIZE}" \
        "2_batch_ProcessBamFiles.sh" \
        --config="${configuration_directory}" \
        --mark="${mark}" \
        --phred="${PRED_SCORE_THRESHOLD}" | \
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
        --time="${MAXTIME_2}" \
        --array=1-"${PROCESSING_ARRAY_SIZE}" \
        --dependency=afterok:"${jobID[${array_index_move}]}" \
        "2_batch_ProcessBamFiles.sh" \
        --config="${configuration_directory}" \
        --mark="${mark}" \
        --phred="${PRED_SCORE_THRESHOLD}" | \
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
        --time="${MAXTIME_3}" \
        "3_SubsampleBamFiles.sh" \
        --config="${configuration_directory}" \
        --mark="${mark}" \
        --samplesize="${SAMPLE_SIZE}" | \
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
        --time="${MAXTIME_3}" \
        --dependency=afterok:"${jobID[${array_index_process}]}" \
        "3_SubsampleBamFiles.sh" \
        --config="${configuration_directory}" \
        --mark="${mark}" \
        --samplesize="${SAMPLE_SIZE}" | \
        awk '{print $4}' \
        )

        echo "Submitted 3_SubsampleBamFiles.sh for $mark under job ID:"
        echo "${jobID[$array_index_merge]}"
    done
fi
## -------------------------------------------------------------------------- ##


## ---4_BinarizeBamFiles.sh-------------------------------------------------- ##
if [[ "${STARTING_SCRIPT}" -eq 4 ]]; then
    jobID[binarization]=$( \
    sbatch \
    --time="${MAXTIME_4}" \
    "4_BinarizeBamFiles.sh" \
    --config="${configuration_directory}" \
    --binsize="${BIN_SIZE}" \
    --samplesize="${SAMPLE_SIZE}" \
    --assembly="${ASSEMBLY}" | \
    awk '{print $4}' \
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
    # Checkpoint script should not finish searching until all subsampling has
    # finished. We use the downloading MAXTIME as this is likely to be very
    # high and wait until at least the final mark has been successfully merged
    # to start the search
    array_index_merge="${LIST_OF_MARKS[-1]}_merge"

    jobID[checkpoint]=$( \
    sbatch \
    --time="${MAXTIME_0}" \
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
    --time="${MAXTIME_4}" \
    --dependency=afterok:"${jobID[checkpoint]}" \
    "4_BinarizeBamFiles.sh" \
    --config="${configuration_directory}" \
    --binsize="${BIN_SIZE}" \
    --samplesize="${SAMPLE_SIZE}" \
    --assembly="${ASSEMBLY}" | \
    awk '{print $4}' \
    )

    echo "Submitted 4_BinarizeBamFiles.sh under job ID:"
    echo "${jobID[binarization]}"
fi
## -------------------------------------------------------------------------- ##


## ---5_batch_CreateIncrementalModels.sh------------------------------------- ##
if [[ "${STARTING_SCRIPT}" -eq 5 ]]; then
    jobID[Model_Learning]=$( \
    sbatch \
    --time="${MAXTIME_5}" \
    --array=1-"${MODEL_LEARNING_ARRAY_SIZE}" \
    "5_batch_CreateIncrementalModels.sh" \
    --config="${configuration_directory}" \
    --nummodels="${NUMBER_OF_MODELS}" \
    --binsize="${BIN_SIZE}" \
    --samplesize="${SAMPLE_SIZE}" \
    --assembly="${ASSEMBLY}" | \
    awk '{print $4}' \
    )

    echo "Submitted 5_batch_CreateIncrementalModels.sh under job ID:"
    echo "${jobID[Model_Learning]}"

elif [[ "${STARTING_SCRIPT}" -lt 5 ]]; then
    jobID[Model_Learning]=$( \
    sbatch \
    --time="${MAXTIME_5}" \
    --array=1-"${MODEL_LEARNING_ARRAY_SIZE}" \
    --dependency=afterok:"${jobID[binarization]}" \
    "5_batch_CreateIncrementalModels.sh" \
    --config="${configuration_directory}" \
    --nummodels="${NUMBER_OF_MODELS}" \
    --binsize="${BIN_SIZE}" \
    --samplesize="${SAMPLE_SIZE}" \
    --assembly="${ASSEMBLY}" | \
    awk '{print $4}' \
    )

    echo "Submitted 5_batch_CreateIncrementalModels.sh under job ID:"
    echo "${jobID[Model_Learning]}"
fi
## -------------------------------------------------------------------------- ##


## ---6_OptimalNumberOfStates.sh--------------------------------------------- ##
if [[ "${STARTING_SCRIPT}" -eq 6 ]]; then
    jobID[Optimal_States]=$( \
    sbatch \
    --time="${MAXTIME_6}" \
    "6_OptimalNumberOfStates.sh" \
    --config="${configuration_directory}" \
    --chromosome="${CHROMOSOME_IDENTIFIER}" \
    --binsize="${BIN_SIZE}" \
    --samplesize="${SAMPLE_SIZE}" \
    --nummodels="${NUMBER_OF_MODELS}" | \
    awk '{print $4}' \
    )

    echo "Submitted 6_OptimalNumberOfStates.sh under job ID:"
    echo "${jobID[Optimal_States]}"

elif [[ "${STARTING_SCRIPT}" -lt 6 ]]; then
    jobID[Optimal_States]=$( \
    sbatch \
    --time="${MAXTIME_6}" \
    --dependency=afterok:"${jobID[Model_Learning]}" \
    "6_OptimalNumberOfStates.sh" \
    --config="${configuration_directory}" \
    --chromosome="${CHROMOSOME_IDENTIFIER}" \
    --binsize="${BIN_SIZE}" \
    --samplesize="${SAMPLE_SIZE}" \
    --nummodels="${NUMBER_OF_MODELS}" | \
    awk '{print $4}' \
    )

    echo "Submitted 6_OptimalNumberOfStates.sh under job ID:"
    echo "${jobID[Optimal_States]}"
fi
## -------------------------------------------------------------------------- ##

## ---7_RunLDSC.sh----------------------------------------------------------- ##
if [[ "${STARTING_SCRIPT}" -eq 7 ]]; then
    jobID[LDSC]=$( \
    sbatch \
    --time="${MAXTIME_7}" \
    "7_RunLDSC.sh" \
    --config="${configuration_directory}" \
    --state="${OPTIMUM_NUMBER_OF_STATES}" \
    --gwas="${GWAS_PATTERN}" \
    --binsize="${BIN_SIZE}" \
    --samplesize="${SAMPLE_SIZE}" \
    --nummodels="${NUMBER_OF_MODELS}" | \
    awk '{print $4}' \
    )

    echo "Submitted 7_RunLDSC.sh under job ID:"
    echo "${jobID[LDSC]}"

elif [[ "${STARTING_SCRIPT}" -lt 7 ]]; then
    jobID[LDSC]=$( \
    sbatch \
    --time="${MAXTIME_7}" \
    --dependency=afterok:"${jobID[Optimal_States]}" \
    "7_RunLDSC.sh" \
    --config="${configuration_directory}" \
    --gwas="${GWAS_PATTERN}" \
    --binsize="${BIN_SIZE}" \
    --samplesize="${SAMPLE_SIZE}" \
    --nummodels="${NUMBER_OF_MODELS}" | \
    awk '{print $4}' \
    )

    echo "Submitted 7_RunLDSC.sh under job ID:"
    echo "${jobID[LDSC]}"
fi
## -------------------------------------------------------------------------- ##
