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
Recreate_Transition_Heatmap.sh
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

module purge
module load R/4.2.1-foss-2022a

Rscript \
    "${RSCRIPTS_DIR}/RecreateTransitionMatrix.R" \
    "${transitions_file}" \
    "${output_directory}"


finishing_statement 0
