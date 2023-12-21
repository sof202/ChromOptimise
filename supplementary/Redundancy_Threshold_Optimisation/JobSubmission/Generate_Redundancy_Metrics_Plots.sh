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
## $1 -> Size of model (default: 20)                                                ||
## $2 -> Random seed (default: 1)"                                                  ||
## $3 -> Path to directory containing the model files                               ||
##       (default: \${BIG_MODELS_DIR} in config.txt)                                ||
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
    echo "\$1 -> Size of model (default: 20)"
    echo "\$2 -> Random seed (default: 1)"
    echo "\$3 -> Path to directory containing the model files"
    echo "       (default: \${BIG_MODELS_DIR} in config.txt)"
    echo "======================================================================="
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
scripts/integrative/ChromHMM_OptimumStates/config/config.txt"

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

## ============================= ##
##    VARIABLES AND FUNCTIONS    ##
## ============================= ##

model_size=$1
seed=$2
model_file_dir=$3

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


if [[ -z "$model_file_dir" ]]; then
    model_file_dir="${BIG_MODELS_DIR}"
    echo "Model file directory was not given, using the default of: ${BIG_MODELS_DIR}"
fi

if [[ -z "$model_size" ]]; then
    model_size=20
    echo "No model size was given by the user, using default value of: ${model_size}."
fi

if [[ -z "$seed" ]]; then
    seed=1
    echo "No random seed was given by the user, using defualt value of: ${seed}." 
fi

## ================== ##
##   FILE EXISTANCE   ##
## ================== ##

cd "${model_file_dir}" ||  { >&2 echo "ERROR: ${model_file_dir} doesn't exist, \
ensure that the directory exists before running this script."; finishing_statement 1; }

if [[ -z $(find . -type f -name "emissions*") ]]; then
    { >&2 echo -e "ERROR: No model files were found in ${model_file_dir}.\n\
    Ensure that you have ran Generate_Big_Model.sh or ChromHMM's \
    LearnModel Command before using this script." 
    }

    finishing_statement 1
fi

## ======== ##
##   MAIN   ##
## ======== ##

module purge
module load R/4.2.1-foss-2022a

cd "${SUPPLEMENTARY_DIR}/Redundancy_Threshold_Optimisation/Rscripts" ||  \
{ >&2 echo "ERROR: ${SUPPLEMENTARY_DIR}/Redundancy_Threshold_Optimisation/Rscripts \
doesn't exist, make sure config.txt is pointing to the correct directory"; \
finishing_statement 1; }

Rscript HistogramPlotForEuclideanDistances.R \
"${model_size}" "${seed}" "${model_file_dir}"

Rscript ScatterPlotForTransitionMaxima.R \
"${model_size}" "${seed}" "${model_file_dir}"

finishing_statement 0
