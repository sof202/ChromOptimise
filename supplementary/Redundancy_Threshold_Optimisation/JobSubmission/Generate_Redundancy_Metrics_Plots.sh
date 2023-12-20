#!/bin/bash
# Export all environment variables to the batch job
#SBATCH --export=ALL 
# Submit to the mrc queue for faster queue times
#SBATCH -p mrcq 
# This script uses two Rscripts with low computational complexity
#SBATCH --time=00:10:00 
#SBATCH -A Research_Project-MRC190311
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
# Very low memory usage is expected from generating the plots
#SBATCH --mem=1G 
#SBATCH --mail-type=END # Send an email after the job is done
# Temporary log file, later to be removed
#SBATCH --output=temp%j.log
# Temporary error file, later to be removed
#SBATCH --error=temp%j.err
#SBATCH --job-name=Redundancy_Threshold_Plotting

## =================================================================================##
##                                                                                  ||
##                                     PREAMBLE                                     ||
##                                                                                  ||
## =================================================================================##
## PURPOSE:                                                                         ||
## Generate plots for the euclidean distances between the emission parameters (for  ||
## pairs of states) and the maximum transition probability towards each state in    ||
## the selected model.                                                              ||
## =================================================================================##
## AUTHOR: Sam Fletcher s.o.fletcher@exeter.ac.uk                                   ||
## CREATED: December 2023                                                           ||
## =================================================================================##
## PREREQUISITES: Run Generate_Big_Model.sh or ChromHMM's LearnModel Command        ||
## =================================================================================##
## DEPENDENCIES: R                                                                  ||
## =================================================================================##
## INPUTS:                                                                          ||
## $1 -> Size of model                                                              ||
## $2 -> Random Seed                                                                ||
## $3 -> Path to directory containing model files                                   ||
## =================================================================================##
## OUTPUTS:                                                                         ||
## Histogram plot of Euclidean distances between emission parameters of pairs       ||
## of states. A threshold suggestion value is also given.                           ||
## A scatter plot of the maximum transition probabilitites towards each state.      ||
## =================================================================================##

## ======================== ##
##    HELP FUNCTIONALITY    ##
## ======================== ##

if [[ "$1" == "--help" || "$1" == "-h" ]]; then
    echo "======================================================================="
    echo "Purpose: Generates plots to aid in thresholds used in config.r which"
    echo "are used in 6_OptimumNumberOfStates.sh in determining redundant states."
    echo "Author: Sam Fletcher"
    echo "Contact: s.o.fletcher@exeter.ac.uk"
    echo "Dependencies: R"
    echo "Inputs:"
    echo "\$1 -> Size of model"
    echo "\$2 -> Random seed"
    echo "\$3 -> Path to directory containing the model files"
    echo "======================================================================="
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

## ========================= ##
##    VARIABLE ASSIGNMENT    ##
## ========================= ##

MODEL_SIZE=$1
SEED=$2
MODEL_FILE_DIR=$3

if [ -z "$MODEL_FILE_DIR" ]; then
    echo "Model file directory was given using the default of ${BIG_MODELS_DIR}"
    MODEL_FILE_DIR="${BIG_MODELS_DIR}"
fi

if [ -z "$MODEL_SIZE" ]; then
    echo "No model size was given by the user, using default value of 20."
    MODEL_SIZE=20
fi

if [ -z "$SEED" ]; then
    echo "No random seed was given by the user, using default value of 1."
    SEED=1
fi

## ================== ##
##   FILE EXISTANCE   ##
## ================== ##

cd "${MODEL_FILE_DIR}" ||  { echo "Directory given doesn't exist, \
ensure that the directory exists or that config.txt is pointing \
to the correct directory if default path was used."; exit 1; }

if [[ -z $(find . -type f -name emissions*) ]]; then
    echo "No model files were found in the directory given."
    echo "Aborting..."

    rm "${SLURM_SUBMIT_DIR}/temp${SLURM_JOB_ID}.log" 
    rm "${SLURM_SUBMIT_DIR}/temp${SLURM_JOB_ID}.err" 

    exit 1
fi

## ======== ##
##   MAIN   ##
## ======== ##

module purge
module load R/4.2.1-foss-2022a

cd "${SUPPLEMENTARY_DIR}/Redundancy_Threshold_Optimisation/Rscripts" || exit 1

Rscript HistogramPlotsForSimilarityMetrics_Emissions.r \
"${MODEL_SIZE}" "${SEED}" "${MODEL_FILE_DIR}"

Rscript HistogramPlotsForSimilarityMetrics_Transistions.r \
"${MODEL_SIZE}" "${SEED}" "${MODEL_FILE_DIR}"


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
