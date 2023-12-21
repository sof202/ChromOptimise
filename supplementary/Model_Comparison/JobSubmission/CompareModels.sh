#!/bin/bash
# Export all environment variables to the batch job
#SBATCH --export=ALL
# Submit to the mrc queue for faster queue times
#SBATCH -p mrcq
# Running this script on model files (2-8 states) took 10 seconds
#SBATCH --time=00:01:00 
#SBATCH -A Research_Project-MRC190311 
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=16
# When comparing a relatively small number of models (7 models) the peak
# heap memory consumption was 44 KB.
#SBATCH --mem=1G
# Send an email after the job is done
#SBATCH --mail-type=END 
# Temporary log file, later to be removed
#SBATCH --output=temp%j.log
# Temporary error file, later to be removed
#SBATCH --error=temp%j.err
#SBATCH --job-name=Model_Comparing

## =================================================================================##
##                                                                                  ||
##                                     PREAMBLE                                     ||
##                                                                                  ||
## =================================================================================##
## PURPOSE:                                                                         ||
## Compares models produced by 5_batch_CreateIncrementalModels.sh using ChromHMM's  ||
## CompareModels command. The base model used for comparing is the most complex     ||                       
## model. The emission file for the most complex model is then deleted and the      ||
## process is repeated for the next most complex model.                             ||
## This continues until all emission files have been deleted.                       ||
## =================================================================================##
## AUTHOR: Sam Fletcher s.o.fletcher@exeter.ac.uk                                   ||
## CREATED: November 2023                                                           ||
## =================================================================================##
## PREREQUISITES: Run: 5_batch_CreateIncrementalModels.sh                           ||
## =================================================================================##
## DEPENDENCIES: Java, ChromHMM                                                     ||
## =================================================================================##
## INPUTS:                                                                          ||
## NONE                                                                             ||
## =================================================================================##
## OUTPUTS:                                                                         ||
## Model comparison files in (.txt,.svg,.png)                                       ||
## =================================================================================##

## ======================== ##
##    HELP FUNCTIONALITY    ##
## ======================== ##

if [[ "$1" == "--help" || "$1" == "-h" ]]; then
    echo "==========================================================================="
    echo "Purpose: Uses ChromHMM's CompareModels to compare model files sequentially."
    echo "Author: Sam Fletcher"
    echo "Contact: s.o.fletcher@exeter.ac.uk"
    echo "Dependencies: Java, ChromHMM"
    echo "Inputs:"
    echo "NONE"
    echo "==========================================================================="
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
/scripts/integrative/blueprint/config/config.txt"

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

## ============= ##
##   FUNCTIONS   ##
## ============= ##

## ====== FUNCTION : delete_logs() ========================
## Delete temporary log and error files then exit
## Globals: 
##   SLURM_SUBMIT_DIR
##   SLURM_JOB_ID
## Arguments:
##   exit code
## ========================================================
delete_logs(){
    rm "${SLURM_SUBMIT_DIR}/temp${SLURM_JOB_ID}.log" 
    rm "${SLURM_SUBMIT_DIR}/temp${SLURM_JOB_ID}.err"
    exit "$1"
}

## ===================== ##
##    FILE MANAGEMENT    ##
## ===================== ##

cd "${MODEL_DIR}" || \
{ >&2 echo "ERROR: \${MODEL_DIR} - ${MODEL_DIR} doesn't exist, \
make sure config.txt is pointing to the correct directory."; delete_logs 1; }

if [ -z "$(ls -A)" ]; then
    { echo -e "ERROR: \${MODEL_DIR} - ${MODEL_DIR} is empty.\n\
    Ensure that 5_CreateIncrementalModels.sh has been ran before this script."; }

    delete_logs 1
fi

cd "${COMPARE_DIR}" || { echo "ERROR: \${COMPARE_DIR} - ${COMPARE_DIR} doesn't exist, \
make sure config.txt is pointing to the correct directory"; delete_logs 1; }
mkdir -p temp
cd temp || delete_logs 1
rm ./*

cd "${MODEL_DIR}" || \
{ >&2 echo "ERROR: \${MODEL_DIR} - ${MODEL_DIR} doesn't exist, \
make sure config.txt is pointing to the correct directory."; delete_logs 1; }

emission_text_files=$(find . -type f -name "Emission*.txt")

echo "Copying all found emission files to a temporary directory..."
for file in $emission_text_files; do
    echo "Copying ${file}..."

    # ChromHMM's CompareModels requires files to start with 'emissions'
    # This is case sensitive, hence we need to convert to lower case.
    new_file_name=$(echo "$file" | tr '[:upper:]' '[:lower:]')
    cp "$file" "${COMPARE_DIR}/temp/${new_file_name}"
done

## -------------------- ##
##   COMPARING MODELS   ##
## -------------------- ##

# Steps:
# 1) Sort the files by the number of states
# 2) Create a comparison file using the most complex model as a base
# 3) Delete the most complex model

module purge
module load Java

for file in $emission_text_files; do
    # 1)
    cd "${COMPARE_DIR}" || delete_logs 1
    most_complex_model_number=$(find ./temp -type f -name "emissions*.txt" | \
    grep -oP "\d+(?=.txt)"| \
    sort -g | \
    tail -1)

    most_complex_model_file=$(find ./temp -type f -name "emissions*.txt" | \
    grep "${most_complex_model_number}.txt")

    echo "${most_complex_model_file}"

    # 2)
    echo -e "Comparing model with ${most_complex_model_file} states to the \ 
    less complex models..."

    java -mx1G \
    -jar "${CHROMHMM_MAIN_DIR}/ChromHMM.jar" CompareModels \
    "${most_complex_model_file}" temp/ \
    "Comparison_To_${most_complex_model_number}_states" 

    # 3)
    rm "${most_complex_model_file}"
done

rm -rf temp

## ======================= ##
##   LOG FILE MANAGEMENT   ##
## ======================= ##

echo "Job completed at:"
date -u
end_time=$(date +%s)
time_taken=$((end_time-start_time))
echo "Job took a total of: ${time_taken} seconds to complete"

delete_logs 0



