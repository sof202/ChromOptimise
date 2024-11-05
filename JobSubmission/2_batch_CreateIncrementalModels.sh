#!/bin/bash
# Export all environment variables to the batch job
#SBATCH --export=ALL
# Submit to the mrc queue for faster queue times
#SBATCH -p mrcq
# Consult information/Processing_Times.md for expected time
# Forward backwards algorithm has time complexity of O(N^2T), where T is the 
# number of genomic bins and N is the number of states
#SBATCH --time=01:00:00 
#SBATCH -A Research_Project-MRC190311 
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=16
# Make sure array is not higher than the number of models being learned
# as this ends up with all models being processed by max index array
#SBATCH --array=1-4
# Memory consumption doesn't appear to be dependent on number of states
# Consult information/Memory_Profiling.md for expected memory usage
#SBATCH --mem=5G
#SBATCH --mail-type=END 
#SBATCH --output=temp%A_%a.log
#SBATCH --error=temp%A_%a.err
#SBATCH --job-name=2_Model_Learning

usage() {
cat <<EOF
================================================================================
$(basename "$0")
================================================================================
Purpose: Uses ChromHMM's LearnModel command to generate several models
with increasing numbers of states.
Author: Sam Fletcher
Contact: s.o.fletcher@exeter.ac.uk
Dependencies: Java, ChromHMM, awk, gzip
Inputs:
\$1 -> Full/relative file path for configuation file directory
================================================================================
EOF
    exit 0
}

if [[ $# -eq 0 ]]; then usage; fi

## ============ ##
##    SET UP    ##
## ============ ##

configuration_directory=$1

source "${configuration_directory}/Config.txt" || \
{ echo "The configuration file does not exist in the specified location: \
${configuration_directory}"; exit 1; }

source "${WRAPPER_SCRIPT}" || exit 1

## =============== ##
##    VARIABLES    ##
## =============== ##

## ====== DEFAULTS =============================================================
if ! [[ "${NUMBER_OF_MODELS}" =~ ^[0-9]+$ ]]; then
    NUMBER_OF_MODELS=4
    echo "Invalid value of number of models to learn given."
    echo "Using the default value of: ${NUMBER_OF_MODELS} instead."
fi

if ! [[ "${MAX_ITERATIONS}" =~ ^[0-9]+$ ]]; then
    MAX_ITERATIONS=200
    echo "Invalid max iterations given, using the default value of" \
    "${MAX_ITERATIONS} instead."
fi

if [[ -z "${ASSEMBLY}" ]]; then
    ASSEMBLY=hg19
    echo "No assembly was given, using the default value of" \
    "${ASSEMBLY} instead."
fi
# ==============================================================================

## =============================== ##
##   CLEAN UP AND ERROR CATCHING   ##
## =============================== ##

full_binary_path="${BINARY_DIR}/BinSize_${BIN_SIZE}"

if [[ ! -d ${full_binary_path} ]]; then
   >&2 echo "ERROR: Binary directory for ${BIN_SIZE} is empty." \
   "Ensure that binary files exist before running this script."
   finishing_statement 1
fi

output_directory="${MODEL_DIR}/BinSize_${BIN_SIZE}_models_${NUMBER_OF_MODELS}"

if [[ "${SLURM_ARRAY_TASK_ID}" -eq 1 ]]; then
    rm -rf "${output_directory}"
    mkdir -p \
        "${output_directory}/Likelihood_Values" \
        "${output_directory}/STATEBYLINE"
fi

# Ensures all array tasks have correct output directory structure before model
# learning (avoiding stale file handles)
sleep 5

## ========================== ##
##   PARALLELISATION SET UP   ##
## ========================== ##

number_of_models_per_array=$((NUMBER_OF_MODELS / SLURM_ARRAY_TASK_COUNT))
remainder=$((NUMBER_OF_MODELS % SLURM_ARRAY_TASK_COUNT))

# By default, we assume that the user wants to learn every model from 2 upwards
states_increment=1

# If the number of models to be generated isn't a multiple of the size of the
# array, the array with the smallest id will learn the left over smaller models
# (spreading the larger models evenly across the other array elements).
if [[ "${SLURM_ARRAY_TASK_ID}" -eq 1 ]]; then
    starting_number_of_states=2
    ending_number_of_states=$(( \
    2 + states_increment*(remainder+number_of_models_per_array-1) \
    ))
else
    starting_number_of_states=$(( \
    (((SLURM_ARRAY_TASK_ID-1)*number_of_models_per_array) + remainder)*\
    states_increment + 2 \
    )) 
    ending_number_of_states=$(( \
    (((SLURM_ARRAY_TASK_ID)*number_of_models_per_array) + remainder -1 )*\
    states_increment + 2 \
    ))
fi

## ===================== ##
##    MODEL GENERATION   ##
## ===================== ##

source "${CONDA_SHELL}" || { echo "Could not find conda shell at:
${CONDA_SHELL}"; exit 1; }
conda deactivate
conda activate ChromOptimise-R-java

sequence=$(\
seq "$starting_number_of_states" "$states_increment" "$ending_number_of_states"\
)

echo "Learning models using a bin size of ${BIN_SIZE}..." 

cd "${output_directory}/Likelihood_Values" || finishing_statement 1

# Fixes problems with X11 servers (when local PC shuts down)
unset DISPLAY
for numstates in ${sequence}; do
    echo "Learning model with: ${numstates} states..."

    # -noautoopen used so html files are not opened after model learning finshes.
    # -printstatebyline used to get the state assignment for isolation metrics
    java \
        -mx4G \
        -jar "${CHROMHMM_MAIN_DIR}/ChromHMM.jar" \
        LearnModel \
        -noautoopen \
        -printstatebyline \
        -b "${BIN_SIZE}" \
        -r "${MAX_ITERATIONS}" \
        "${full_binary_path}" \
        "${output_directory}" \
        "${numstates}" \
        "${ASSEMBLY}" > \
        "ChromHMM_output_numstates_${numstates}.txt"

    echo "Writing estimated log likelihood to: likelihoods.txt..."
    echo -n "Estimated Log Likelihood for ${numstates} states: " >> \
    "likelihoods.txt"

    # grep selects terminal logs that are not associated with writing to files.
    grep "  " "ChromHMM_output_numstates_${numstates}.txt" | \
    tail -1 | \
    awk '{print $2}' >> \
    "likelihoods.txt"

    # Instead of storing chromHMM's log file in a separate location, it is
    # easier to just store the log in the existing log file
    grep "  " "ChromHMM_output_numstates_${numstates}.txt" >> \
    "${LOG_FILE_PATH}/${SLURM_ARRAY_JOB_ID}~${SLURM_ARRAY_TASK_ID}~${timestamp}.log"

    rm "ChromHMM_output_numstates_${numstates}.txt"
done

conda deactivate

## ============ ##
##   CLEAN UP   ##
## ============ ##

# html files are not required for subsequent analysis
rm "${output_directory}"/*.html

finishing_statement 0
