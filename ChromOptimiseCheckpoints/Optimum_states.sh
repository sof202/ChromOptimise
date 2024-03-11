#!/bin/bash
# Export all environment variables to the batch job
#SBATCH --export=ALL
# Submit to the mrc queue for faster queue times
#SBATCH -p mrcq 
# Script should take seconds (just reads a file)
#SBATCH --time=00:01:00 
#SBATCH -A Research_Project-MRC190311 
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=16 
# Predicted that memory consumption will rise massively when merging lots of files
#SBATCH --mem=1G 
# Send an email after the job is done
#SBATCH --mail-type=END 
# Temporary log file, later to be removed
#SBATCH --output=temp%j.log
# Temporary error file, later to be removed
#SBATCH --error=temp%j.err
#SBATCH --job-name=Subsampling_Checkpoint

## ===========================================================================##
##                                                                            ||
##                                  PREAMBLE                                  ||
##                                                                            ||
## ===========================================================================##
## PURPOSE:                                                                   ||
## This script simply checks the output of the optimum number of states script||
## and stores it in a temporary file for the master script to utilise         ||
## ===========================================================================##
## AUTHOR: Sam Fletcher                                                       ||
## CONTACT: s.o.fletcher@exeter.ac.uk                                         ||
## CREATED: March 2023                                                        ||
## ===========================================================================##
## PREREQUISITES: Run 6_OptimalNumberOfStates.sh                              ||
## ===========================================================================##
## DEPENDENCIES: NONE                                                         ||
## ===========================================================================##
## INPUTS:                                                                    ||
## $1 -> Full (or relative) file path for configuation file directory         ||
## $2 -> Number of models made                                                ||
## $3 -> State increment used                                                 ||
## $4 -> Bin size                                                             ||
## $5 -> Sample Size                                                          ||
## ===========================================================================##
## OUTPUTS:                                                                   ||
## NONE                                                                       ||
## ===========================================================================##

## ============ ##
##    SET UP    ##
## ============ ##

# Configuration files are required for file paths and log file management
configuration_directory=$1

source "${configuration_directory}/FilePaths.txt" || \
{ echo "The configuration file does not exist in the specified location: \
${configuration_directory}"; exit 1; }



# Log and error files are not required for this script as all relevant
# information is displayed in the masterscript
rm "${SLURM_SUBMIT_DIR}/temp${SLURM_JOB_ID}.log" 
rm "${SLURM_SUBMIT_DIR}/temp${SLURM_JOB_ID}.err"

## ============= ##
##   VARIABLES   ##
## ============= ##

number_of_models=$2
states_increment=$3
bin_size=$4
sample_size=$5

## ======== ##
##   MAIN   ##
## ======== ##

max_model_size=$(( 2 + (number_of_models - 1) * states_increment ))

input_file="${OPTIMUM_STATES_DIR}\
/BinSize_${bin_size}_SampleSize_${sample_size}_MaxModelSize_${max_model_size}\
/OptimumNumberOfStates.txt"

optimum_number_of_states=$(tail -1 "${input_file}" | cut -d: -f2 | tr -d ' ')

echo "${optimum_number_of_states}" > "$(dirname "$0")/../Optimum_state.txt"
