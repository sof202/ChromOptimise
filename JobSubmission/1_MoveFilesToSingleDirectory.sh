#!/bin/bash
# Export all environment variables to the batch job
#SBATCH --export=ALL
# Submit to the mrc queue for faster queue times
#SBATCH -p mrcq
# Job is unlikely to take a long time unless moving lots of very large files
#SBATCH --time=01:00:00
#SBATCH -A Research_Project-MRC190311 
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=16
# Moving files on linux takes very little memory
#SBATCH --mem=1G
# Send an email after the job is done
#SBATCH --mail-type=END
# Temporary log file, later to be removed
#SBATCH --output=temp%j.log
# Temporary error file, later to be removed
#SBATCH --error=temp%j.err
#SBATCH --job-name=Moving_Files

## =================================================================================##
##                                                                                  ||
##                                     PREAMBLE                                     ||
##                                                                                  ||
## =================================================================================##
## PURPOSE:                                                                         ||
## Moves .bam files that include the epigenetic mark into a single folder.          ||
## =================================================================================##
## AUTHOR: Sam Fletcher s.o.fletcher@exeter.ac.uk                                   ||
## CREATED: November 2023                                                           ||
## =================================================================================##
## PREREQUISITES: Downloaded files from EGA                                         ||
## =================================================================================##
## DEPENDENCIES: NONE                                                               ||
## =================================================================================##
## INPUTS:                                                                          ||
## $1 -> Epigenetic Mark to process                                                 ||
## =================================================================================##
## OUTPUTS:                                                                         ||
## NONE                                                                             ||
## =================================================================================##

## ======================== ##
##    HELP FUNCTIONALITY    ##
## ======================== ##

if [[ "$1" == "--help" || "$1" == "-h" ]]; then
    echo "=========================================================="
    echo "Purpose: Moves .bam files that include the epigenetic mark"
    echo "into a single folder."
    echo "Author: Sam Fletcher"
    echo "Contact: s.o.fletcher@exeter.ac.uk"
    echo "Dependencies: NONE"
    echo "Inputs:"
    echo "\$1 -> Epigenetic mark name"
    echo "=========================================================="
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
# [epigenetic mark name]~[job id]~[date]-[time]
# This requires a hard link as you cannot rename log files
# whilst running the script without a wrapper function
LOG_FILE_PATH="${LOG_DIR}/$SLURM_JOB_NAME/$USER"
mkdir -p "${LOG_FILE_PATH}"
timestamp=$(date -u +%Y.%m.%d-%H:%M)

ln "${SLURM_SUBMIT_DIR}/temp${SLURM_JOB_ID}.log" \
"${LOG_FILE_PATH}/$1~${SLURM_JOB_ID}~$timestamp.log"
ln "${SLURM_SUBMIT_DIR}/temp${SLURM_JOB_ID}.err" \
"${LOG_FILE_PATH}/$1~${SLURM_JOB_ID}~$timestamp.err"

## ============================= ##
##    VARIABLES AND FUNCTIONS    ##
## ============================= ##

mark_name=$1

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

## ================== ##
##    MOVING FILES    ##
## ================== ##

cd "${MAIN_DIR}" || \
{ >&2 echo "ERROR: \${MAIN_DIR} - ${MAIN_DIR} doesn't exist, \
make sure config.txt is pointing to the correct directory"; finishing_statement 1; }

list_of_files_with_mark_name=\
$(find . -type f -name "*${mark_name}*.bam")

# If the above list is empty then the epigenetic mark must not exist.
number_of_files_to_move=\
$(find . -type f -name "*${mark_name}*.bam" | wc -l)
echo "Number of .bam files to be moved is: ${number_of_files_to_move}"

# Exit here if there are no files 
# There is no point in generating a folder if the mark doesn't exist
if [[ "${number_of_files_to_move}" -eq 0 ]]; then
    { >&2 echo -e "ERROR: No files with epigenetic mark: ${mark_name} were found.\n\
    Please input a epigenetic mark name that exists (note that this is \
    case sensitive)."; }

    finishing_statement 1
fi

mkdir -p "${RAW_DIR}/${mark_name}"

echo "Moving .bam files to ${RAW_DIR}/${mark_name}..."

for file in ${list_of_files_with_mark_name}; do
    echo "Moving file ${file}..."
    mv "${file}" "${RAW_DIR}/${mark_name}"
done

finishing_statement 0