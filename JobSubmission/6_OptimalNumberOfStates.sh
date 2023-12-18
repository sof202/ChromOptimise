#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -p mrcq # submit to the mrc queue for faster queue times
#SBATCH --time=01:00:00 
#SBATCH -A Research_Project-MRC190311 
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=16 
#SBATCH --mem=10G # specify bytes memory to reserve
#SBATCH --mail-type=END # Send an email after the job is done
# Temporary log file, later to be removed
#SBATCH --output=temp%j.log
# Temporary error file, later to be removed
#SBATCH --error=temp%j.err
#SBATCH --job-name=Optimal_States

## -------------------------------------------------------------------------------------------- ##
##                                                                                              ##
##                                            PREAMBLE                                          ##
##                                                                                              ##
## -------------------------------------------------------------------------------------------- ##
##                                            PURPOSE                                           ##
## The models have at this point been learned and now is the final stage. The optimal number of ##
##        states needs to be determined. The way this is accomplished is by searching for       ##
## 'redundant states'. By using the compare models command in ChromHMM, we can see how similar  ## 
## certain states are amongst models. If we compare the model with the maximum number of states ##
##                    against the other models we will see one of two scenarios:                ##
##                                                                                              ##
## 1) Multiple states will be explained by a single state in a less complex model => redundant  ## 
##                               states exist in the more complex model                         ##
##   2) Every state in the more complex model is described by AT MOST one state in each less    ##
##    complex model => The model has the optimum number of states, less complex models have     ##
##                        lower complexity but don't capture the whole system                   ##
##                                                                                              ##
##                 Note that, by the pigeonhole principle, the above is a dichotomy.            ##
##                                                                                              ##
## In addition to the above, if a model with less states has a larger estimated log likelihood  ##
## than a more complex model, the more complex model will be rejected (punishing complexity to  ##
##                         avoid overfitting due to needless complexity).                       ## 
## -------------------------------------------------------------------------------------------- ##
##                        AUTHOR: Sam Fletcher s.o.fletcher@exeter.ac.uk                        ##
##                                     CREATED: November 2023                                   ##
## -------------------------------------------------------------------------------------------- ##
##                                         PREREQUISITES                                        ##
##                            Run: 5_batch_CreateIncrementalModels.sh                           ##
##                                    Run: 6_CompareModels.sh                                   ##
## -------------------------------------------------------------------------------------------- ##
##                                          DEPENDENCIES                                        ##
##                                               R                                              ##
## -------------------------------------------------------------------------------------------- ##
##                                            INPUTS                                            ##
##            $1 -> Bin Size (make sure to be consistent in choice between scripts)             ##
##          $2 -> Sample Size (make sure to be consistent in choice between scripts)            ##
## -------------------------------------------------------------------------------------------- ##
##                                            OUTPUTS                                           ##
##                  File containing why models with too many states were rejected               ##
##                      The optimum number of states to use with the model                      ##
##       A line graph displaying the relationship between estimated log likelihood and          ##
##                                       number of states                                       ##
## -------------------------------------------------------------------------------------------- ##

## ------------------------ ##
##    HELP FUNCTIONALITY    ##
## ------------------------ ##

if [[ "$1" == "--help" || "$1" == "-h" ]]; then
    echo "==========================================================================================="
    echo "Purpose: Determines the optimum number of states to use with the data"
    echo "Author: Sam Fletcher"
    echo "Contact: s.o.fletcher@exeter.ac.uk"
    echo "Dependencies: R"
    echo "Inputs:"
    echo "\$1 -> Bin size, WARNING: Use the same bin size as was used in 4_BinarizeBamFiles.sh"
    echo "\$2 -> Sample size, WARNING: Use the same sample size as was used in 3_SubsampleBamFiles.sh"
    echo "==========================================================================================="
    exit 0
fi

## ------------ ##
##    SET UP    ##
## ------------ ##

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


# Variables
BIN_SIZE=$1
SAMPLE_SIZE=$2

cd "${MODEL_DIR}" || { echo "Model directory doesn't exist, make sure config.txt is pointing to the correct directory"; exit 1; }
if [ -z "$(ls -A)" ]; then
    echo "No files found in the model directory, please run 5_CreateIncrementalModels.sh before this script."
    echo "Aborting..."

    # Remove temporary log files
    rm "${SLURM_SUBMIT_DIR}/temp${SLURM_JOB_ID}.log"
    rm "${SLURM_SUBMIT_DIR}/temp${SLURM_JOB_ID}.err"

    exit 1
fi

# Set defaults if no arguments are given by the user 'intelligently' by searching through the model directory
if [ -z "$BIN_SIZE" ]; then
    BIN_SIZE=$(find . -type f -name "Emissions*.txt" | head -1 | cut -d "_" -f 3)
    echo "No bin size was given, assuming a default value of ${BIN_SIZE}..."
fi
if [ -z "$SAMPLE_SIZE" ]; then
    SAMPLE_SIZE=$(find . -type f -name "Emissions*.txt" | head -1 | cut -d "_" -f 5)
    echo "No sample size was given, assuming a default value of ${SAMPLE_SIZE}..."
fi

## ------------------- ##
##   FILE MANAGEMENT   ##
## ------------------- ##

# Remove any contents of temporary folder
mkdir -p "${OPTIMUM_STATES_DIR}/temp"
cd "${OPTIMUM_STATES_DIR}/temp" || exit 1
rm -f ./*

cd "${MODEL_DIR}" || { echo "Model directory doesn't exist, make sure config.txt is pointing to the correct directory"; exit 1; }
Emission_Text_Files=$(find . -type f -name "Emissions*.txt")
for file in $Emission_Text_Files; do
    cp "$file" "${OPTIMUM_STATES_DIR}/temp"
done
Transition_Text_Files=$(find . -type f -name "Transitions*.txt")
for file in $Transition_Text_Files; do
    cp "$file" "${OPTIMUM_STATES_DIR}/temp"
done

# Comparison files [CURRENTLY NOT IN USE]
# cd ${COMPARE_DIR}
# Comparison_Text_Files=$(find . -type f -name "Comparison*")
# for file in $Comparison_Text_Files; do
#     cp $file "${OPTIMUM_STATES_DIR}/temp"
# done


## --------------- ##
##    MAIN LOOP    ##
## --------------- ##


module purge
module load R/4.2.1-foss-2022a
cd "${RSCRIPTS_DIR}" || { echo "Rscripts directory doesn't exist, make sure config.txt is pointing to the correct directory"; exit 1; }

Max_Model_Number=$(find "${OPTIMUM_STATES_DIR}/temp" -type f -name "*.txt" | grep -oP "\d+(?=.txt)"| sort -g | tail -1) # gets the number of states in each model, sorts them and takes the largest value

Output_Directory="${OPTIMUM_STATES_DIR}/BinSize_${BIN_SIZE}_SampleSize_${SAMPLE_SIZE}_MaxModelSize_${Max_Model_Number}"
mkdir -p "${Output_Directory}"
rm -f "${Output_Directory}"/*

while [[ $Max_Model_Number -gt 2 ]]; do
    Max_Model_Number=$(find "${OPTIMUM_STATES_DIR}/temp" -type f -name "*.txt" | grep -oP "\d+(?=.txt)"| sort -g | tail -1) 

    Rscript RedundantStateChecker.R "${Max_Model_Number}" "${BIN_SIZE}" "${SAMPLE_SIZE}" "${Output_Directory}"
    Redundant_States=$(tail -1 "${Output_Directory}/Redundant_States_Modelsize_${Max_Model_Number}.txt")

    if [[ "$Redundant_States" == "NONE" ]]; then
        echo "Model with ${Max_Model_Number} states has no redundant states." >> "${Output_Directory}/OptimumNumberOfStates.BinSize.${BIN_SIZE}.SampleSize.${SAMPLE_SIZE}.txt"
        break
    else
        rm -f "${OPTIMUM_STATES_DIR}"/temp/*"${Max_Model_Number}".txt
        echo "Model with ${Max_Model_Number} states has redundant states: ${Redundant_States}" >> "${Output_Directory}/OptimumNumberOfStates.BinSize.${BIN_SIZE}.SampleSize.${SAMPLE_SIZE}.txt"
    fi
done

rm -r "${OPTIMUM_STATES_DIR}/temp"


## -------------- ##
##    PLOTTING    ##
## -------------- ##

# Plots the estimated log likelihood against the number of states across all models 
Rscript PlotLikelihoods.R "${BIN_SIZE}" "${SAMPLE_SIZE}" "${Output_Directory}"

## ----------------------- ##
##   LOG FILE MANAGEMENT   ##
## ----------------------- ##

# Finishing message
echo "Job completed at:"
date -u
end_time=$(date +%s)
time_taken=$((end_time-start_time))
echo "Job took a total of: ${time_taken} seconds to complete"

# Removing temporary log files
rm "${SLURM_SUBMIT_DIR}/temp${SLURM_JOB_ID}.log"
rm "${SLURM_SUBMIT_DIR}/temp${SLURM_JOB_ID}.err"

