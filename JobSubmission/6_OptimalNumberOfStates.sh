#!/bin/bash
# Export all environment variables to the batch job
#SBATCH --export=ALL 
# Submit to the mrc queue for faster queue times
#SBATCH -p mrcq 
# Script usually takes less than one minute
#SBATCH --time=00:10:00 
#SBATCH -A Research_Project-MRC190311 
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=16
# Memory consumption is very low in testing
# Consult information/Memory_Profiling.md for expected memory usage
#SBATCH --mem=1G 
# Send an email after the job is done
#SBATCH --mail-type=END 
# Temporary log file, later to be removed
#SBATCH --output=temp%j.log
# Temporary error file, later to be removed
#SBATCH --error=temp%j.err
#SBATCH --job-name=Optimal_States

## =================================================================================##
##                                                                                  ||
##                                     PREAMBLE                                     ||
##                                                                                  ||
## =================================================================================##
## PURPOSE:                                                                         ||
## Determines the optimum number of states by searching for redundant states in     ||
## the model files (starting with most complex). Redundant states are states that   ||
## satisfy the following criteria:                                                  ||
##  (i) The state's emissions parameter vector is close to another state's under    ||
##      the Euclidean distance metric,                                              ||
## (ii) The state's transition parameter vector (towards the state) has a low       ||
##      maximum value.                                                              ||
## If a model has redundant states it is rejected in favour of a simpler model.     ||
## This then repeats, iterating across smaller and smaller models until no          ||
## no redundant states are found.                                                   ||
##                                                                                  ||
## The script also creates a plot of the estimated log likelihood against the       ||
## number of states in each model for human sense checking.                         ||
##                                                                                  ||
## Note: If the largest model has no redundant states, the optimum model size may   ||
##       be larger than the largest model that was trained.                         ||
## =================================================================================##
## AUTHOR: Sam Fletcher                                                             ||
## CONTACT: s.o.fletcher@exeter.ac.uk                                               ||
## CREATED: December 2023                                                           ||
## =================================================================================##
## PREREQUISITES: Run: 5_batch_CreateIncrementalModels.sh                           ||
## =================================================================================##
## DEPENDENCIES: R                                                                  ||
## =================================================================================##
## INPUTS:                                                                          ||
## NONE                                                                             ||
## =================================================================================##
## OUTPUTS:                                                                         ||
## File containing why models with too many states were rejected                    ||
## The optimum number of states to use with the data                                ||
## Plot between estimated log likelihood and number of states                       ||
## =================================================================================##


## ======================== ##
##    HELP FUNCTIONALITY    ##
## ======================== ##

if [[ "$1" == "--help" || "$1" == "-h" ]]; then
    echo "======================================================================"
    echo "Purpose: Determines the optimum number of states to use with the data."
    echo "Author: Sam Fletcher"
    echo "Contact: s.o.fletcher@exeter.ac.uk"
    echo "Dependencies: R"
    echo "Inputs:"
    echo "NONE"
    echo "======================================================================"
    exit 0
fi

## ============ ##
##    SET UP    ##
## ============ ##

echo "Job '${SLURM_JOB_NAME}' started at:"
date -u

start_time=$(date +%s)

# Activate config.txt to access all file paths
# CHANGE THIS TO YOUR OWN CONFIG FILE
source "/lustre/projects/Research_Project-MRC190311\
/scripts/integrative/ChromHMM_OptimumStates/config/config.txt"

# Rename the output and error files to have format:
# [job id]~[date]-[time]
# This requires a hard link as you cannot rename log files
# whilst running the script without a wrapper function
LOG_FILE_PATH="${LOG_DIR}/$SLURM_JOB_NAME/$USER"
mkdir -p "${LOG_FILE_PATH}"
timestamp=$(date -u +%Y.%m.%d-%H:%M)

ln "${SLURM_SUBMIT_DIR}/temp${SLURM_JOB_ID}.log" \
"${LOG_FILE_PATH}/${SLURM_JOB_ID}~$timestamp.log"
ln "${SLURM_SUBMIT_DIR}/temp${SLURM_JOB_ID}.err" \
"${LOG_FILE_PATH}/${SLURM_JOB_ID}~$timestamp.err"


## ============================= ##
##    VARIABLES AND FUNCTIONS    ##
## ============================= ##

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

# Set bin/sample size by searching through the model directory
cd "${MODEL_DIR}" || \
{ >&2 echo "ERROR: \${MODEL_DIR} - ${MODEL_DIR} doesn't exist, \
make sure config.txt is pointing to the correct directory."; finishing_statement 1; }

bin_size=$(find . -type f -name "*.txt" | head -1 | cut -d "_" -f 3)
sample_size=$(find . -type f -name "*.txt" | head -1 | cut -d "_" -f 5)

## =================== ##
##   FILE MANAGEMENT   ##
## =================== ##

if [[ -z "$(ls -A)" ]]; then
    { >&2 echo -e "ERROR: No files found in \${MODEL_DIR} - ${MODEL_DIR}.\n\
    Please run 5_CreateIncrementalModels.sh before this script."
    finishing_statement 1; }
fi

mkdir -p "${OPTIMUM_STATES_DIR}/temp"
cd "${OPTIMUM_STATES_DIR}/temp" || finishing_statement 1
rm -f ./*

cd "${MODEL_DIR}" || \
{ >&2 echo "ERROR: \${MODEL_DIR} - ${MODEL_DIR} doesn't exist, \
make sure config.txt is pointing to the correct directory."; finishing_statement 1; }

emission_text_files=$(find . -type f -name "Emissions*.txt")
for file in $emission_text_files; do
    cp "$file" "${OPTIMUM_STATES_DIR}/temp"
done
transition_text_files=$(find . -type f -name "Transitions*.txt")
for file in $transition_text_files; do
    cp "$file" "${OPTIMUM_STATES_DIR}/temp"
done

## =============== ##
##    MAIN LOOP    ##
## =============== ##

module purge
module load R/4.2.1-foss-2022a

cd "${RSCRIPTS_DIR}" || \
{ >&2 echo "\${RSCRIPTS_DIR} - ${RSCRIPTS_DIR} doesn't exist, \
make sure config.txt is pointing to the correct directory"; finishing_statement 1; }

max_model_number=$(find "${OPTIMUM_STATES_DIR}/temp" -type f -name "*.txt" | \
grep -oP "\d+(?=.txt)"| \
sort -g | \
tail -1) 

output_directory="${OPTIMUM_STATES_DIR}\
/BinSize_${bin_size}_SampleSize_${sample_size}_MaxModelSize_${max_model_number}"

mkdir -p "${output_directory}"
rm -f "${output_directory}"/*

while [[ ${max_model_number} -gt 2 ]]; do
    max_model_number=$(find "${OPTIMUM_STATES_DIR}/temp" -type f -name "*.txt" | \
    grep -oP "\d+(?=.txt)"| \
    sort -g | \
    tail -1) 

    echo "Running RedundantStateChecker.R for: ${max_model_number}..."

    Rscript RedundantStateChecker.R "${max_model_number}" "${bin_size}" \
    "${sample_size}" "${output_directory}"

    redundant_states=$(tail -1 \
    "${output_directory}/Redundant_States_Modelsize_${max_model_number}.txt")

    if [[ "$redundant_states" == "NONE" ]]; then
        echo "Model with ${max_model_number} states has no redundant states." >> \
        "${output_directory}/OptimumNumberOfStates.txt"
        break
    else
        rm -f "${OPTIMUM_STATES_DIR}"/temp/*"${max_model_number}".txt
        echo -n "Model with ${max_model_number} states has redundant states: " >> \
        "${output_directory}/OptimumNumberOfStates.txt"
        echo "${redundant_states}" >> "${output_directory}/OptimumNumberOfStates.txt"
    fi
done

rm -r "${OPTIMUM_STATES_DIR}/temp"

## ========================= ##
##   OPTIMUM STATES CHECK    ##
## ========================= ##

# If the largest model learned has no redundant states, this doesn't necessarily imply
# that it has the optimum number of states, perhaps a more complex model does. This
# section checks for this scenario.

if [[ $(wc -l "${output_directory}/OptimumNumberOfStates.txt") -eq 1 ]]; then
    {
    echo -e "${max_model_number} states may not be the optimum number of states.\n\ 
    Try increasing the size of the most complex model or increasing \
    the thresholds in the config.r file." 
    } >> "${output_directory}/OptimumNumberOfStates.txt"
else
    echo "Optimum number of states for the data is: ${max_model_number}" >> \
    "${output_directory}/OptimumNumberOfStates.txt"
fi

## ============== ##
##    PLOTTING    ##
## ============== ##

# Plots the estimated log likelihood against the number of states across all models
echo "Plotting the estimated log likelihoods of learned models against one another..."
Rscript PlotLikelihoods.R "${bin_size}" "${sample_size}" "${output_directory}"

finishing_statement 0

