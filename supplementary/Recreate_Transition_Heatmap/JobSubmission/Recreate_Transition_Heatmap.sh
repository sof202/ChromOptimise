#!/bin/bash
#SBATCH --export=ALL
#SBATCH -p mrcq 
#SBATCH --time=00:05:00
#SBATCH -A Research_Project-MRC190311 
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=16
#SBATCH --mem=1G
#SBATCH --mail-type=END 
#SBATCH --output=Recreate_Heatmap%j.log
#SBATCH --error=Recreate_Heatmap%j.err
#SBATCH --job-name=Recreate_Heatmap


usage() {
cat << EOF
================================================================================
$(basename "$0")
================================================================================
Purpose: Creates a heatmap for transition parameters without the main diagonal
Author: Sam Fletcher
Contact: s.o.fletcher@exeter.ac.uk
Inputs:
\$1 -> configuration directory
\$2 -> path to transitions file
\$3 -> output directory
================================================================================
EOF
    exit 0
}

configuration_directory=$1
transitions_file=$2
output_directory=$3

source "${configuration_directory}/Config.txt" || \
{ echo "The configuration file does not exist in the specified location: \
${configuration_directory}"; exit 1; }

source "${WRAPPER_SCRIPT}" || exit 1

if [[ $# -ne 3 ]]; then usage; fi

## ======== ##
##   MAIN   ##
## ======== ##

source "${CONDA_SHELL}" || { echo "Could not find conda shell at:
${CONDA_SHELL}"; exit 1; }
conda deactivate
conda activate ChromOptimise-R-java

Rscript \
    "${RSCRIPTS_DIR}/RecreateTransitionMatrix.R" \
    "${REPO_DIR}" \
    "${transitions_file}" \
    "${output_directory}"

conda deactivate

finishing_statement 0
