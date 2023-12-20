#!/bin/bash
# Export all environment variables to the batch job
#SBATCH --export=ALL 
# Submit to the mrc queue for faster queue times
#SBATCH -p mrcq 
# This can be a lengthy process for particularly large models
# A model with 80 states took 26 hours off of 1MB of binary data.
# See BigModelLearningTimes.pdf for an idea of maximum wall time
#SBATCH --time=40:00:00 
#SBATCH -A Research_Project-MRC190311
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
# Possible large memory consumption for big models, difficult to give good estimates,  
# see 5_batch_CreateIncrementalModels.sh for a idea with less states
#SBATCH --mem=32G 
#SBATCH --mail-type=END # Send an email after the job is done
# Temporary log file, later to be removed
#SBATCH --output=temp%j.log
# Temporary error file, later to be removed
#SBATCH --error=temp%j.err
#SBATCH --job-name=Big_Model

## =================================================================================##
##                                                                                  ||
##                                     PREAMBLE                                     ||
##                                                                                  ||
## =================================================================================##
## PURPOSE:                                                                         ||
## Generate a very large model using the binarized data generated in                ||
## 4_BinarizeBamFiles.sh with a random initialisation of parameters                 ||
## (with a set seed).                                                               ||
## =================================================================================##
## AUTHOR: Sam Fletcher s.o.fletcher@exeter.ac.uk                                   ||
## CREATED: November 2023                                                           ||
## =================================================================================##
## PREREQUISITES: Run 4_BinarizeBamFiles.sh                                         ||
## =================================================================================##
## DEPENDENCIES: Java, ChromHMM                                                     ||
## =================================================================================##
## INPUTS:                                                                          ||
## $1 -> Size of model                                                              ||
## $2 -> Random Seed                                                                ||
## =================================================================================##
## OUTPUTS:                                                                         ||
## The emission parameter matrix of the model (.png,.txt,.svg)                      ||
## The transition parameter matrix of the model (.png,.txt,.svg)                    ||
## Full model file                                                                  ||
## The overlap and fold enrichment files for existing genomic annotations           ||
## The log files produced by ChromHMM.jar                                           ||
## The maximum achieved estimated log likelihood achieved by the model              ||
## =================================================================================##

## ======================== ##
##    HELP FUNCTIONALITY    ##
## ======================== ##

if [[ "$1" == "--help" || "$1" == "-h" ]]; then
    echo "================================================================"
    echo "Purpose: Generates a model unrestricted by ChromHMM LearnModel's"
    echo "default size limits by using random initialisation."
    echo "Author: Sam Fletcher"
    echo "Contact: s.o.fletcher@exeter.ac.uk"
    echo "Dependencies: Java, ChromHMM"
    echo "Inputs:"
    echo "\$1 -> Size of model"
    echo "\$2 -> Random seed"
    echo "================================================================"
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
# ModelSize-[model size]~[job id]~[date]-[time]
# This requires a hard link as you cannot rename log files
# whilst running the script without a wrapper function
LOG_FILE_PATH="${LOG_DIR}/$SLURM_JOB_NAME/$USER"
mkdir -p "${LOG_FILE_PATH}"
timestamp=$(date -u +%Y.%m.%d-%H:%M)

ln "${SLURM_SUBMIT_DIR}/temp${SLURM_JOB_ID}.log" \
"${LOG_FILE_PATH}/ModelSize-$1~${SLURM_JOB_ID}~$timestamp.log"
ln "${SLURM_SUBMIT_DIR}/temp${SLURM_JOB_ID}.err" \
"${LOG_FILE_PATH}/ModelSize-$1~${SLURM_JOB_ID}~$timestamp.err"

cd "${BINARY_DIR}" || { echo "Binary directory doesn't exist, \
make sure config.txt is pointing to the correct directory"; exit 1; }
if [ -z "$(ls -A)" ]; then
    echo "4_BinarizedFiles is empty."
    echo "Ensure that 4_BinarizeBamFiles.sh has been ran before this script."
    echo "Aborting..."

    rm "${SLURM_SUBMIT_DIR}/temp${SLURM_JOB_ID}.log" 
    rm "${SLURM_SUBMIT_DIR}/temp${SLURM_JOB_ID}.err" 

    exit 1
fi

## ========================= ##
##    VARIABLE ASSIGNMENT    ##
## ========================= ##

MODEL_SIZE=$1
SEED=$2

if [ -z "$MODEL_SIZE" ]; then
    echo "No model size was given by the user, using default value of 20."
    MODEL_SIZE=20
fi

if [ -z "$SEED" ]; then
    echo "No random seed was given by the user, using defualt value of 1."
    SEED=1
fi

echo "Learning a model with ${MODEL_SIZE} states and with random seed: ${SEED}."

## ========== ##
##    MAIN    ##
## ========== ##

module purge
module load Java

cd "${BIG_MODELS_DIR}" || { echo "Big models directory doesn't exist, \
make sure config.txt is pointing to the correct directory"; exit 1; }

java -mx30G \
-jar "${CHROMHMM_MAIN_DIR}/ChromHMM.jar" LearnModel \
-noautoopen \
-nobed \
-init random \
-s "${SEED}" \
"${BINARY_DIR}" "${BLUEPRINT_MAIN_DIR}/Big_Model_Files" "${MODEL_SIZE}" hg19 > \
"ChromHMM.Output.ModelSize.${MODEL_SIZE}.txt"

echo "Model Learning finished"
echo "Writing estimated log likelihood to file"
# grep removes the terminal logs associated with writing to files.
# The tail and awk locate the final estimated log likelihood
grep "       " "ChromHMM.Output.ModelSize.${MODEL_SIZE}.txt" | \
tail -1 | \
awk '{print $2}' >> \
"likelihood.ModelSize.${MODEL_SIZE}.txt" 

## ======================= ##
##   LOG FILE MANAGEMENT   ##
## ======================= ##

echo "Job completed at:"
date -u
end_time=$(date +%s)
time_taken=$((end_time - start_time))
echo "Job took a total of: ${time_taken} seconds to complete."

rm "${SLURM_SUBMIT_DIR}/temp${SLURM_JOB_ID}.log" 
rm "${SLURM_SUBMIT_DIR}/temp${SLURM_JOB_ID}.err" 
