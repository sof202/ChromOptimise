#!/bin/bash
# Export all environment variables to the batch job
#SBATCH --export=ALL 
# Submit to the mrc queue for faster queue times
#SBATCH -p mrcq 
# Script usually takes less than three minutes
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
## AUTHOR: Sam Fletcher                                                             ||
## CONTACT: s.o.fletcher@exeter.ac.uk                                               ||
## CREATED: December 2023                                                           ||
## =================================================================================##
## PREREQUISITES: Run Generate_Big_Model.sh or ChromHMM's LearnModel command        ||
## =================================================================================##
## DEPENDENCIES: R                                                                  ||
## =================================================================================##
## INPUTS:                                                                          ||
## $1 -> Size of model (default: 20)                                                ||
## $2 -> Random seed (default: 1)                                                   ||
## $3 -> Path to directory containing the model files                               ||
##       (default: \${BIG_MODELS_DIR} in FilePaths.txt)                             ||
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
    echo "Purpose: Generates plots to aid in thresholds used in config.R which"
    echo "are used in 6_OptimumNumberOfStates.sh in determining redundant states."
    echo "Author: Sam Fletcher"
    echo "Contact: s.o.fletcher@exeter.ac.uk"
    echo "Dependencies: R"
    echo "Inputs:"
    echo "\$1 -> Size of model (default: 20)"
    echo "\$2 -> Random seed (default: 1)"
    echo "\$3 -> Path to directory containing the model files"
    echo "       (default: \${BIG_MODELS_DIR} in FilePaths.txt)"
    echo "======================================================================="
    exit 0
fi

## ============ ##
##    SET UP    ##
## ============ ##

# CHANGE THESE TO YOUR OWN CONFIG FILES
source "/lustre/projects/Research_Project-MRC190311/scripts/integrative\
/ChromOptimise/configuration/FilePaths.txt"
source "/lustre/projects/Research_Project-MRC190311/scripts/integrative\
/ChromOptimise/configuration/LogFileManagement.sh"

# Output and error files renamed to:
# ModelSize-[model size]~[job id]~[date]-[time]

ln "${SLURM_SUBMIT_DIR}/temp${SLURM_JOB_ID}.log" \
"${LOG_FILE_PATH}/ModelSize-$1~${SLURM_JOB_ID}~${timestamp:=}.log"
ln "${SLURM_SUBMIT_DIR}/temp${SLURM_JOB_ID}.err" \
"${LOG_FILE_PATH}/ModelSize-$1~${SLURM_JOB_ID}~$timestamp.err"

## =============== ##
##    VARIABLES    ##
## =============== ##

model_size=$1
seed=$2
model_file_dir=$3

## ====== DEFAULTS ====================================================================
if [[ -z "$model_file_dir" ]]; then
    model_file_dir="${BIG_MODELS_DIR}"
    echo "Model file directory was not given, using the default of: ${BIG_MODELS_DIR}"
fi

if ! [[ "${model_size}" =~ ^[0-9]+$ ]]; then
    model_size=20
    echo "Model size given is invalid, using default value of: ${model_size}."
fi

if ! [[ "$seed" =~ ^[0-9]+$ ]]; then
    seed=1
    echo "Random seed given is invalid, using defualt value of: ${seed}." 
fi
# =====================================================================================

## ================== ##
##   FILE EXISTANCE   ##
## ================== ##

cd "${model_file_dir}" ||  { >&2 echo "ERROR: ${model_file_dir} doesn't exist, \
ensure that the directory exists before running this script."; finishing_statement 1; }

if [[ -z $(find . -type f -name "emissions*") ]]; then
    { >&2 echo -e "ERROR: No model files were found in ${model_file_dir}.\n"\
    "Ensure that you have ran Generate_Big_Model.sh or ChromHMM's "\
    "LearnModel command before using this script."; finishing_statement 1; }  
fi

## ======== ##
##   MAIN   ##
## ======== ##

module purge
module load R/4.2.1-foss-2022a

cd "${SUPPLEMENTARY_DIR}/Redundancy_Threshold_Optimisation/Rscripts" ||  \
{ >&2 echo "ERROR: ${SUPPLEMENTARY_DIR}/Redundancy_Threshold_Optimisation/Rscripts \
doesn't exist, make sure [\${SUPPLEMENTARY_DIR} - ${SUPPLEMENTARY_DIR}] \
in config.txt is pointing to the correct directory"; finishing_statement 1; }

Rscript HistogramPlotForEuclideanDistances.R \
"${model_size}" "${seed}" "${model_file_dir}"

Rscript ScatterPlotForTransitionMaxima.R \
"${model_size}" "${seed}" "${model_file_dir}"

finishing_statement 0
