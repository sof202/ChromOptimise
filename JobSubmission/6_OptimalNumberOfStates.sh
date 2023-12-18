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
# Script uses very little memory consumption
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
## satisfy the following critereon:                                                 ||
##  (i) The state's emissions parameter vector is close to another state's under    ||
##      the Euclidean distance metric,                                              ||
## (ii) The state's transition parameter vector (towards the state) has a low       ||
##      maximum value.                                                              ||
## If a model has redundant states it is rejected in favour of a simpler model.     ||
## This then repeats, iterating across smaller and smaller models until no          ||
## no redundant states are found.                                                   ||
##                                                                                  ||
## The script also creates a plot of the log likelihood against the number of       ||
## states in each model for human sense checking.                                   ||
##                                                                                  ||
## Note: If the largest model has no redundant states, the optimum model size may   ||
##       be larger than the largest model that was trained                          ||
## =================================================================================##
## AUTHOR: Sam Fletcher s.o.fletcher@exeter.ac.uk                                   ||
## CREATED: December 2023                                                           ||
## =================================================================================##
## PREREQUISITES: Run: 5_batch_CreateIncrementalModels.sh                           ||
## =================================================================================##
## DEPENDENCIES: R                                                                  ||
## =================================================================================##
## INPUTS:                                                                          ||
## $3 -> Bin size, WARNING: Use the same bin size as was used in                    ||
##       4_BinarizeBamFiles.sh                                                      ||
## $4 -> Sample Size, WARNING: Use the same sample size as was used in              ||
##       3_SubsampleBamFiles.sh                                                     ||
## =================================================================================##
## OUTPUTS:                                                                         ||
## File containing why models with too many states were rejected                    ||
## The optimum number of states to use with the model                               ||
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
    echo "\$1 -> Bin size, WARNING: Use the same bin size as was used in"
    echo "4_BinarizeBamFiles.sh"
    echo "\$2 -> Sample size, WARNING: Use the same sample size as was used in"
    echo "3_SubsampleBamFiles.sh"
    echo "======================================================================"
    exit 0
fi

## ============ ##
##    SET UP    ##
## ============ ##

echo "Job '$SLURM_JOB_NAME' started at:"
date -u

start_time=$(date +%s)

# Activate config.txt to access all file paths
# CHANGE THIS TO YOUR OWN CONFIG FILE
echo "Loading config file..."
source "/lustre/projects/Research_Project-MRC190311\
/scripts/integrative/blueprint/config/config.txt"

# Rename the output and error files to have format:
# [job id]~[date]-[time]
# This requires a hard link as you cannot rename log files
# whilst running the script without a wrapper function
LOG_FILE_PATH="${LOG_DIR}/$USER/$SLURM_JOB_NAME/"
mkdir -p "${LOG_FILE_PATH}"
timestamp=$(date -u +%Y.%m.%d-%H:%M)

ln "${SLURM_SUBMIT_DIR}/temp${SLURM_JOB_ID}.log" \
"${LOG_FILE_PATH}/${SLURM_JOB_ID}~$timestamp.log"
ln "${SLURM_SUBMIT_DIR}/temp${SLURM_JOB_ID}.err" \
"${LOG_FILE_PATH}/${SLURM_JOB_ID}~$timestamp.err"


## ========================= ##
##    VARIABLE ASSIGNMENT    ##
## ========================= ##

BIN_SIZE=$1
SAMPLE_SIZE=$2

cd "${MODEL_DIR}" || { echo "Model directory doesn't exist, \
make sure config.txt is pointing to the correct directory"; exit 1; }
if [ -z "$(ls -A)" ]; then
    echo "No files found in the model directory."
    echo "Please run 5_CreateIncrementalModels.sh before this script."
    echo "Aborting..."

    rm "${SLURM_SUBMIT_DIR}/temp${SLURM_JOB_ID}.log"
    rm "${SLURM_SUBMIT_DIR}/temp${SLURM_JOB_ID}.err"

    exit 1
fi

# 'Intelligently' set defaults by searching through the model directory
if [ -z "$BIN_SIZE" ]; then
    BIN_SIZE=$(find . -type f -name "Emissions*.txt" | head -1 | cut -d "_" -f 3)
    echo "No bin size was given, assuming a default value of ${BIN_SIZE}..."
fi

if [ -z "$SAMPLE_SIZE" ]; then
    SAMPLE_SIZE=$(find . -type f -name "Emissions*.txt" | head -1 | cut -d "_" -f 5)
    echo "No sample size was given, assuming a default value of ${SAMPLE_SIZE}..."
fi

## =================== ##
##   FILE MANAGEMENT   ##
## =================== ##

mkdir -p "${OPTIMUM_STATES_DIR}/temp"
cd "${OPTIMUM_STATES_DIR}/temp" || exit 1
rm -f ./*

cd "${MODEL_DIR}" || { echo "Model directory doesn't exist, \
make sure config.txt is pointing to the correct directory"; exit 1; }
Emission_Text_Files=$(find . -type f -name "Emissions*.txt")
for file in $Emission_Text_Files; do
    cp "$file" "${OPTIMUM_STATES_DIR}/temp"
done
Transition_Text_Files=$(find . -type f -name "Transitions*.txt")
for file in $Transition_Text_Files; do
    cp "$file" "${OPTIMUM_STATES_DIR}/temp"
done

## =============== ##
##    MAIN LOOP    ##
## =============== ##

module purge
module load R/4.2.1-foss-2022a
cd "${RSCRIPTS_DIR}" || { echo "Rscripts directory doesn't exist, \
make sure config.txt is pointing to the correct directory"; exit 1; }

Max_Model_Number=$(find "${OPTIMUM_STATES_DIR}/temp" -type f -name "*.txt" | \
grep -oP "\d+(?=.txt)"| \
sort -g | \
tail -1) 

Output_Directory="${OPTIMUM_STATES_DIR}\
/BinSize_${BIN_SIZE}_SampleSize_${SAMPLE_SIZE}_MaxModelSize_${Max_Model_Number}"

mkdir -p "${Output_Directory}"
rm -f "${Output_Directory}"/*

while [[ $Max_Model_Number -gt 2 ]]; do
    Max_Model_Number=$(find "${OPTIMUM_STATES_DIR}/temp" -type f -name "*.txt" | \
    grep -oP "\d+(?=.txt)"| \
    sort -g | \
    tail -1) 

    Rscript RedundantStateChecker.R "${Max_Model_Number}" "${BIN_SIZE}" \
    "${SAMPLE_SIZE}" "${Output_Directory}"

    Redundant_States=$(tail -1 "${Output_Directory}\
    /Redundant_States_Modelsize_${Max_Model_Number}.txt")

    if [[ "$Redundant_States" == "NONE" ]]; then
        echo "Model with ${Max_Model_Number} states has no redundant states." >> \
        "${Output_Directory}\
        /OptimumNumberOfStates.BinSize.${BIN_SIZE}.SampleSize.${SAMPLE_SIZE}.txt"
        break
    else
        rm -f "${OPTIMUM_STATES_DIR}"/temp/*"${Max_Model_Number}".txt
        echo -n "Model with ${Max_Model_Number} states has redundant states: "
        echo "${Redundant_States}" >> "${Output_Directory}\
        /OptimumNumberOfStates.BinSize.${BIN_SIZE}.SampleSize.${SAMPLE_SIZE}.txt"
    fi
done

rm -r "${OPTIMUM_STATES_DIR}/temp"


## ============== ##
##    PLOTTING    ##
## ============== ##

# Plots the estimated log likelihood against the number of states across all models 
Rscript PlotLikelihoods.R "${BIN_SIZE}" "${SAMPLE_SIZE}" "${Output_Directory}"

## ======================= ##
##   LOG FILE MANAGEMENT   ##
## ======================= ##

echo "Job completed at:"
date -u
end_time=$(date +%s)
time_taken=$((end_time-start_time))
echo "Job took a total of: ${time_taken} seconds to complete"

rm "${SLURM_SUBMIT_DIR}/temp${SLURM_JOB_ID}.log"
rm "${SLURM_SUBMIT_DIR}/temp${SLURM_JOB_ID}.err"

