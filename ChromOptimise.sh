#!/bin/bash

## ========= ##
##   SETUP   ##
## ========= ##
#
cat << START_MESSAGE
Config file used:
$(cat "${configuration_directory}/Config.txt")

R config used:
$(cat "${configuration_directory}/Config.R")

START_MESSAGE

# We need an associative array to store the job ids for each SLURM job
# that are submitted. This is so that we can use these job ids as dependencies.
# Dependencies allow us to queue up job submissions and run certain scripts
# in parallel whenever required

declare -A jobID

## =========================== ##
##   LOAD CONFIGURATION FILE   ##
## =========================== ##

configuration_directory=$1

source "${configuration_directory}/Config.txt" || \
{ echo "The configuration file does not exist in the specified location: \
${configuration_directory}/Config.txt"; exit 1; }

source "${configuration_directory}/Config.R" || \
{ echo "The configuration file does not exist in the specified location: \
${configuration_directory}/Config.R"; exit 1; }


## ===================== ##
##   RUN CHROMOPTIMISE   ##
## ===================== ##

## ---1_BinarizeBamFiles.sh-------------------------------------------------- ##
if [[ "${STARTING_SCRIPT}" -eq 1 ]]; then
    jobID[binarization]=$( \
    sbatch \
    --time="${MAXTIME_2}" \
    "${JOBSUBMISSION_DIR}/1_BinarizeBamFiles.sh" \
    --config="${configuration_directory}" \
    --binsize="${BIN_SIZE}" \
    --samplesize="${SAMPLE_SIZE}" \
    --assembly="${ASSEMBLY}" | \
    awk '{print $4}' \
    )

    echo "Submitted 1_BinarizeBamFiles.sh under job ID:"
    echo "${jobID[binarization]}"
fi
## -------------------------------------------------------------------------- ##


## ---2_batch_CreateIncrementalModels.sh------------------------------------- ##
if [[ "${STARTING_SCRIPT}" -eq 2 ]]; then
    jobID[Model_Learning]=$( \
    sbatch \
    --time="${MAXTIME_2}" \
    --array=1-"${MODEL_LEARNING_ARRAY_SIZE}" \
    "${JOBSUBMISSION_DIR}/2_batch_CreateIncrementalModels.sh" \
    --config="${configuration_directory}" \
    --nummodels="${NUMBER_OF_MODELS}" \
    --binsize="${BIN_SIZE}" \
    --samplesize="${SAMPLE_SIZE}" \
    --assembly="${ASSEMBLY}" | \
    awk '{print $4}' \
    )

    echo "Submitted 2_batch_CreateIncrementalModels.sh under job ID:"
    echo "${jobID[Model_Learning]}"

elif [[ "${STARTING_SCRIPT}" -lt 2 ]]; then
    jobID[Model_Learning]=$( \
    sbatch \
    --time="${MAXTIME_2}" \
    --array=1-"${MODEL_LEARNING_ARRAY_SIZE}" \
    --dependency=afterok:"${jobID[binarization]}" \
    "${JOBSUBMISSION_DIR}/2_batch_CreateIncrementalModels.sh" \
    --config="${configuration_directory}" \
    --nummodels="${NUMBER_OF_MODELS}" \
    --binsize="${BIN_SIZE}" \
    --samplesize="${SAMPLE_SIZE}" \
    --assembly="${ASSEMBLY}" | \
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
    --config="${configuration_directory}" \
    --chromosome="${CHROMOSOME_IDENTIFIER}" \
    --binsize="${BIN_SIZE}" \
    --samplesize="${SAMPLE_SIZE}" \
    --nummodels="${NUMBER_OF_MODELS}" | \
    awk '{print $4}' \
    )

    echo "Submitted 3_OptimalNumberOfStates.sh under job ID:"
    echo "${jobID[Optimal_States]}"

elif [[ "${STARTING_SCRIPT}" -lt 3 ]]; then
    jobID[Optimal_States]=$( \
    sbatch \
    --time="${MAXTIME_4}" \
    --dependency=afterok:"${jobID[Model_Learning]}" \
    "${JOBSUBMISSION_DIR}/3_OptimalNumberOfStates.sh" \
    --config="${configuration_directory}" \
    --chromosome="${CHROMOSOME_IDENTIFIER}" \
    --binsize="${BIN_SIZE}" \
    --samplesize="${SAMPLE_SIZE}" \
    --nummodels="${NUMBER_OF_MODELS}" | \
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
    --config="${configuration_directory}" \
    --state="${OPTIMUM_NUMBER_OF_STATES}" \
    --binsize="${BIN_SIZE}" \
    --samplesize="${SAMPLE_SIZE}" \
    --nummodels="${NUMBER_OF_MODELS}" | \
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
    --config="${configuration_directory}" \
    --binsize="${BIN_SIZE}" \
    --samplesize="${SAMPLE_SIZE}" \
    --nummodels="${NUMBER_OF_MODELS}" | \
    awk '{print $4}' \
    )

    echo "Submitted 4_ReferenceLDSCore.sh under job ID:"
    echo "${jobID[LDSC]}"
fi
## -------------------------------------------------------------------------- ##
