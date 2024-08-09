#!/bin/bash

## ========= ##
##   SETUP   ##
## ========= ##

usage() {
cat << EOF
================================================================================
ChromOptimise.sh
================================================================================
Purpose: Runs the ChromOptimise pipeline using the specified configuration files
Author: Sam Fletcher
Contact: s.o.fletcher@exeter.ac.uk
Inputs:
\$1 -> Full/relative file path for configuation file directory
================================================================================
EOF
    exit 0
}

if [[ $# -ne 1 ]]; then usage; fi

configuration_directory=$1

source "${configuration_directory}/Config.txt" || \
{ echo "The configuration file does not exist in the specified location: \
${configuration_directory}/Config.txt"; exit 1; }

# We need an associative array to store the job ids for each SLURM job
# that are submitted. This is so that we can use these job ids as dependencies.
# Dependencies allow us to queue up job submissions and run certain scripts
# in parallel whenever required

declare -A jobID

## ===================== ##
##   RUN CHROMOPTIMISE   ##
## ===================== ##

## ---1_BinarizeBamFiles.sh-------------------------------------------------- ##
if [[ "${STARTING_SCRIPT}" -eq 1 ]]; then
    jobID[binarization]=$( \
    sbatch \
    --time="${MAXTIME_1}" \
    "${JOBSUBMISSION_DIR}/1_BinarizeFiles.sh" \
    "${configuration_directory}" | \
    awk '{print $4}' \
    )

    echo "Submitted 1_BinarizeFiles.sh under job ID:"
    echo "${jobID[binarization]}"
fi
## -------------------------------------------------------------------------- ##


## ---2_batch_CreateIncrementalModels.sh------------------------------------- ##
if [[ "${STARTING_SCRIPT}" -eq 2 ]]; then
    jobID[Model_Learning]=$( \
    sbatch \
    --time="${MAXTIME_2}" \
    --array=1-"${MODEL_LEARNING_ARRAY_SIZE}" \
    --mem=$((MODEL_LEARNING_ARRAY_SIZE * 4))G \
    "${JOBSUBMISSION_DIR}/2_batch_CreateIncrementalModels.sh" \
    "${configuration_directory}" | \
    awk '{print $4}' \
    )

    echo "Submitted 2_batch_CreateIncrementalModels.sh under job ID:"
    echo "${jobID[Model_Learning]}"

elif [[ "${STARTING_SCRIPT}" -lt 2 ]]; then
    jobID[Model_Learning]=$( \
    sbatch \
    --time="${MAXTIME_2}" \
    --array=1-"${MODEL_LEARNING_ARRAY_SIZE}" \
    --mem=$((MODEL_LEARNING_ARRAY_SIZE * 4))G \
    --dependency=afterok:"${jobID[binarization]}" \
    "${JOBSUBMISSION_DIR}/2_batch_CreateIncrementalModels.sh" \
    "${configuration_directory}" | \
    awk '{print $4}' \
    )

    echo "Submitted 2_batch_CreateIncrementalModels.sh under job ID:"
    echo "${jobID[Model_Learning]}"
fi
## -------------------------------------------------------------------------- ##


## ---3_OptimalNumberOfStates.sh--------------------------------------------- ##
if [[ "${STARTING_SCRIPT}" -eq 3 ]]; then
    jobID[Optimal_States]=$( \
    sbatch \
    --time="${MAXTIME_3}" \
    "${JOBSUBMISSION_DIR}/3_OptimalNumberOfStates.sh" \
    "${configuration_directory}" | \
    awk '{print $4}' \
    )

    echo "Submitted 3_OptimalNumberOfStates.sh under job ID:"
    echo "${jobID[Optimal_States]}"

elif [[ "${STARTING_SCRIPT}" -lt 3 ]]; then
    jobID[Optimal_States]=$( \
    sbatch \
    --time="${MAXTIME_3}" \
    --dependency=afterok:"${jobID[Model_Learning]}" \
    "${JOBSUBMISSION_DIR}/3_OptimalNumberOfStates.sh" \
    "${configuration_directory}" | \
    awk '{print $4}' \
    )

    echo "Submitted 3_OptimalNumberOfStates.sh under job ID:"
    echo "${jobID[Optimal_States]}"
fi
## -------------------------------------------------------------------------- ##

## ---4_ReferenceLDSCore.sh-------------------------------------------------- ##
if [[ "${STARTING_SCRIPT}" -eq 4 ]]; then
    jobID[LDSC]=$( \
    sbatch \
    --time="${MAXTIME_4}" \
    "${JOBSUBMISSION_DIR}/4_ReferenceLDSCore.sh" \
    "${configuration_directory}" | \
    awk '{print $4}' \
    )

    echo "Submitted 4_ReferenceLDSCore.sh under job ID:"
    echo "${jobID[LDSC]}"

elif [[ "${STARTING_SCRIPT}" -lt 4 ]]; then
    jobID[LDSC]=$( \
    sbatch \
    --time="${MAXTIME_4}" \
    --dependency=afterok:"${jobID[Optimal_States]}" \
    "${JOBSUBMISSION_DIR}/4_ReferenceLDSCore.sh" \
    "${configuration_directory}" | \
    awk '{print $4}' \
    )

    echo "Submitted 4_ReferenceLDSCore.sh under job ID:"
    echo "${jobID[LDSC]}"
fi
## -------------------------------------------------------------------------- ##

## ---5_PartitionedHeritability.sh------------------------------------------- ##
if [[ "${STARTING_SCRIPT}" -eq 5 ]]; then
    jobID[Heritability]=$( \
    sbatch \
    --time="${MAXTIME_5}" \
    "${JOBSUBMISSION_DIR}/5_PartitionedHeritability.sh" \
    "${configuration_directory}" | \
    awk '{print $4}' \
    )

    echo "Submitted 5_PartitionedHeritability.sh under job ID:"
    echo "${jobID[Heritability]}"

elif [[ "${STARTING_SCRIPT}" -lt 5 ]]; then
    jobID[Heritability]=$( \
    sbatch \
    --time="${MAXTIME_5}" \
    --dependency=afterok:"${jobID[LDSC]}" \
    "${JOBSUBMISSION_DIR}/5_PartitionedHeritability.sh" \
    "${configuration_directory}" | \
    awk '{print $4}' \
    )

    echo "Submitted 5_PartitionedHeritability.sh under job ID:"
    echo "${jobID[Heritability]}"
fi
## -------------------------------------------------------------------------- ##
