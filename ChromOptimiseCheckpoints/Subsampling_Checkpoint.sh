#!/bin/bash
# Export all environment variables to the batch job
#SBATCH --export=ALL
# Submit to the mrc queue for faster queue times
#SBATCH -p mrcq 
# Script should have the same time alotted as 3_SubsampleBamFiles.sh
#SBATCH --time=48:00:00 
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

## =================================================================================##
##                                                                                  ||
##                                     PREAMBLE                                     ||
##                                                                                  ||
## =================================================================================##
## PURPOSE:                                                                         ||
## This checks if all marks that have been selected have been selected. The script  ||
## is required for ChromOptimise.sh (master script). The script serves as a flag    ||
## (or checkpoint) for when the binarization process should begin.                  ||
## =================================================================================##
## AUTHOR: Sam Fletcher                                                             ||
## CONTACT: s.o.fletcher@exeter.ac.uk                                               ||
## CREATED: February 2023                                                           ||
## =================================================================================##
## PREREQUISITES: Run 3_SubsampleBamFiles.sh                                        ||
## =================================================================================##
## DEPENDENCIES: NONE                                                               ||
## =================================================================================##
## INPUTS:                                                                          ||
## $1 -> Full (or relative) file path for configuation file directory               ||
## $2 -> Sample size                                                                ||
## remaining -> List of epigenetic marks that are being subsampled                  ||
## =================================================================================##
## OUTPUTS:                                                                         ||
## NONE                                                                             ||
## =================================================================================##

## ============ ##
##    SET UP    ##
## ============ ##

# Configuration files are required for file paths and log file management
configuration_directory=$1; shift

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

sample_size=$1; shift
list_of_marks=( "$@" )

number_of_marks=${#list_of_marks[@]}

cd "${SUBSAMPLED_DIR}" || \
{ >&2 echo "ERROR: [\${SUBSAMPLED_DIR} - ${SUBSAMPLED_DIR} doesn't exist,] \
make sure FilePaths.txt is pointing to the correct directory"; finishing_statement 1; }

## ================================= ##
##   REGULAR EXPRESSION PRODUCTION   ##
## ================================= ##

# We want a regular expression that will capture all files that satisfy
# "Subsampled.{samplesize}.{markname}.bam" for all of the marks that are in 
# list_of_marks

# To make such a regular expression we need to convert our array of marks
# into a pipe-separated list of marks
# i.e. (1 2 3 4 5) --> 1|2|3|4|5
pipe_separated_marks=""


echo "${list_of_marks[2]}"

for mark in "${list_of_marks[@]}"; do
    pipe_separated_marks="${pipe_separated_marks}${mark}|"
done

# The above loops results in a pipe character at the end of the list of marks
# We don't want this as leaving it in will result in the regular expression
# matching a file with name: Subsampled.{samplesize}..bam
pipe_separated_marks="${pipe_separated_marks%|}"



regex_expected_file_names="Subsampled\.${sample_size}\.(${pipe_separated_marks})*"

## ============= ##
##   FUNCTIONS   ##
## ============= ##

## ====== FUNCTION : count_subsampled_files() =========================================
## Description: Find the number of files that satisfy:
## "Subsampled.samplesize.markname.bam"
## Globals: 
##     regular_expression
##     count_of_subsampled_files
## ===================================================================================
count_subsampled_files(){
    count_subsampled_files=$(find . | grep -cE "$regex_expected_file_names")
    echo "$count_subsampled_files"
}


## ============= ##
##   MAIN LOOP   ##
## ============= ##

count_subsampled_files

while [[ "${count_subsampled_files}" -ne "${number_of_marks}" ]]; do
    # We do not need to check for file updates every fraction of a second
    # every minute is good enough
    sleep 60
    count_subsampled_files
done