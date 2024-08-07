#!/bin/bash
# Export all environment variables to the batch job
#SBATCH --export=ALL 
# Submit to the mrc queue for faster queue times
#SBATCH -p mrcq 
# This can be a lengthy process for particularly large models
# Consult information/Processing_Times.md for expected time
#SBATCH --time=40:00:00 
#SBATCH -A Research_Project-MRC190311
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
# Possible large memory consumption for big models, difficult to give good estimates
#SBATCH --mem=32G 
#SBATCH --mail-type=END # Send an email after the job is done
# Temporary log file, later to be removed
#SBATCH --output=temp%j.log
# Temporary error file, later to be removed
#SBATCH --error=temp%j.err
#SBATCH --job-name=Big_Model

usage() {
cat <<EOF
===========================================================================
Generate_Big_Model
===========================================================================
Purpose: Generates a model unrestricted by ChromHMM LearnModel's
default size limits by using random initialisation.
Author: Sam Fletcher
Contact: s.o.fletcher@exeter.ac.uk
Dependencies: Java, ChromHMM
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

## ============== ##
##    VARIABLES   ##
## ============== ##

## ====== DEFAULTS =============================================================
if ! [[ "${model_size}" =~ ^[0-9]+$ ]]; then
    model_size=20
    echo "Model size given is invalid, using default value of: ${model_size}."
fi

# For reproducability
seed=1

if [[ -z "${ASSEMBLY}" ]]; then
    ASSEMBLY=hg19
    echo "No assembly was given, using the default value of" \
    "${ASSEMBLY} instead."
fi
# ==============================================================================

## ================== ##
##   FILE EXISTANCE   ##
## ================== ##

full_binary_path="${BINARY_DIR}/BinSize_${BIN_SIZE}"

cd "${full_binary_path}" || \
{ >&2 echo -e "ERROR: Binary directory for bin/sample size is empty.\n" \
"Ensure that 4_BinarizeBamFiles.sh has been ran before this script."
finishing_statement 1; }

## ========== ##
##    MAIN    ##
## ========== ##

echo "Learning a model with ${model_size} states and random seed: ${seed}..."

module purge
module load Java

cd "${BIG_MODELS_DIR}" || \
{ >&2 echo "ERROR: [\${BIG_MODELS_DIR} - ${BIG_MODELS_DIR}] \
doesn't exist, make sure Config.txt is pointing to the correct directory."
finishing_statement 1; }

java \
    -mx30G \
    -jar "${CHROMHMM_MAIN_DIR}/ChromHMM.jar" LearnModel \
    -noautoopen \
    -nobed \
    -printstatebyline \
    -init random \
    -s "${seed}" \
    "${full_binary_path}" \
    "${BIG_MODELS_DIR}" \
    "${model_size}" \
    "${ASSEMBLY}" > \
    "ChromHMM.Output.ModelSize.${model_size}.txt"

# grep removes the terminal logs associated with writing to files.
# The tail and awk locate the final estimated log likelihood
grep "       " "ChromHMM.Output.ModelSize.${model_size}.txt" | \
    tail -1 | \
    awk '{print $2}' >> \
    "likelihood.ModelSize.${model_size}.txt" 

grep "       " "ChromHMM.Output.ModelSize.${model_size}.txt" >> \
    "${LOG_FILE_PATH}/ModelSize-$1~${SLURM_JOB_ID}~${timestamp}.log"

rm "ChromHMM.Output.ModelSize.${model_size}.txt"

finishing_statement 0
