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

usage() {
cat <<EOF
===========================================================================
Execute_Redundancy_Metrics
===========================================================================
Purpose: Generates plots to aid in thresholds used in config.R which
are used in 6_OptimumNumberOfStates.sh in determining redundant states.
Author: Sam Fletcher
Contact: s.o.fletcher@exeter.ac.uk
Dependencies: R
Inputs:
\$1 -> Full/relative file path for configuation file directory
\$2 -> Size of model to learn
===========================================================================
EOF
    exit 0
}

if [[ $# -eq 0 ]]; then usage; fi

## ============ ##
##    SET UP    ##
## ============ ##

configuration_directory=$1
model_size=$2

source "${configuration_directory}/Config.txt" || \
{ echo "The configuration file does not exist in the specified location: \
${configuration_directory}"; exit 1; }

source "${WRAPPER_SCRIPT}" || exit 1

## =============== ##
##    VARIABLES    ##
## =============== ##

## ====== DEFAULTS =============================================================
seed=1

if ! [[ "${model_size}" =~ ^[0-9]+$ ]]; then
    model_size=20
    echo "Model size given is invalid, using default value of: ${model_size}."
fi

if [[ -z "${CHROMOSOME_IDENTIFIER}" ]]; then
    CHROMOSOME_IDENTIFIER=1
    echo "No chromosome identifier was given, using the default value of:" \
    "${CHROMOSOME_IDENTIFIER} instead."
fi
# ==============================================================================

## ================== ##
##   FILE EXISTANCE   ##
## ================== ##

# The R scripts below will fail to execute if there are no model files
# in the directory given. This just reduces the number of error messages
# and is more descriptive than what R will output.
if [[ -z $(find "${BIG_MODELS_DIR}" -type f -name "emissions*") ]]; then
    { >&2 echo -e "ERROR: No model files were found in ${BIG_MODELS_DIR}.\n"\
    "Ensure that you have ran Generate_Big_Model.sh or ChromHMM's "\
    "LearnModel command before using this script."; finishing_statement 1; }  
fi

## ======== ##
##   MAIN   ##
## ======== ##

module purge
module load R/4.2.1-foss-2022a

cd "${RSCRIPTS_DIR}" || \
{ >&2 echo "ERROR: make sure [\${RSCRIPTS_DIR} - ${RSCRIPTS_DIR}] \
in FilePaths.txt is pointing to the correct directory"; finishing_statement 1; }

mkdir -p \
    "${BIG_MODELS_DIR}/Isolation_scores" \
    "${BIG_MODELS_DIR}/Euclidean_distances" \
    "${BIG_MODELS_DIR}/Flanking_states"

# State assignments are named:
# CellType_BinSize_ModelSize_Seed_Chromosome_statebyline.txt

# We only look at one chromosome as the decision of how to handle the
# following case is rather arbitrary (see wiki): 
# "A state is not assigned on one chromosome but has dense assignment
# on another"
state_assignment_file=$( \
find "${BIG_MODELS_DIR}" \
-name "*${model_size}_${seed}_chr${CHROMOSOME_IDENTIFIER}_*" \
)

emissions_file=$(\
find "${BIG_MODELS_DIR}" \
-name "emissions*${model_size}_${seed}.txt*" \
)

transitions_file=$(\
find "${BIG_MODELS_DIR}" \
-name "transitions*${model_size}_${seed}.txt*" \
)

if [[ -z "$state_assignment_file" ]]; then
    { >&2 echo "ERROR: No state assignment file found for chromosome:" \
    "${CHROMOSOME_IDENTIFIER}, please check ${BIG_MODELS_DIR}/STATEBYLINE for" \
    "the existence of this state assignment file"; finishing_statement 1; }
fi

echo "Running SimilarEmissions.R for: ${model_size} states..."

Rscript SimilarEmissions.R \
    "${emissions_file}" \
    "${BIG_MODELS_DIR}/Euclidean_distances" \
    TRUE


echo "Running FlankingStates.R for: ${model_size} states..."

Rscript FlankingStates.R \
    "${transitions_file}" \
    "${model_file_dir}/Flanking_states"


# IsolationScores.R is ran with a sample size of 100% 
# (all data is considered) because the slow down is not that significant
echo "Running IsolationScores.R for: ${model_size} states..."

Rscript IsolationScores.R \
    "${state_assignment_file}" \
    "${model_file_dir}/Isolation_scores" \
    "${model_size}" \
    100 \
    TRUE


finishing_statement 0
