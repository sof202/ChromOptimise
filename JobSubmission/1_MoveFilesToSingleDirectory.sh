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
# If the downloaded files are very large, you may want to increase this value 
#SBATCH --mem=5G
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
## $1 -> Blueprint Epigenetic Mark to process                                       ||
## =================================================================================##
## OUTPUTS:                                                                         ||
## NONE                                                                             ||
## =================================================================================##

## ======================== ##
##    HELP FUNCTIONALITY    ##
## ======================== ##

if [[ "$1" == "--help" || "$1" == "-h" ]]; then
    echo "==============================================================="
    echo -n "Purpose: Moves .bam files that include the epigenetic mark."
    echo "into a single folder."
    echo "Author: Sam Fletcher"
    echo "Contact: s.o.fletcher@exeter.ac.uk"
    echo "Dependencies: NONE"
    echo "Inputs:"
    echo "\$1 -> Epigenetic mark name"
    echo "==============================================================="
    exit 0
fi

## ============ ##
##    SET UP    ##
## ============ ##

# Print start date/time
echo "Job '$SLURM_JOB_NAME' started at:"
date -u

# Get the start time for the program
start_time=$(date +%s)

# Activate config.txt to access all file paths
# CHANGE THIS TO YOUR OWN CONFIG 
echo "Loading config file..."
source "/lustre/projects/Research_Project-MRC190311\
/scripts/integrative/blueprint/config/config.txt"

# Rename the output and error files to have format:
# [epigenetic mark name]~[job id]~[date]-[time]
# This requires a hard link as you cannot rename log files
# whilst running the script without a wrapper function
LOG_FILE_PATH="${LOG_DIR}/$USER/$SLURM_JOB_NAME/"
timestamp=$(date -u +%Y.%m.%d-%H:%M)
ln "${SLURM_SUBMIT_DIR}/temp${SLURM_JOB_ID}.log" \
"${LOG_FILE_PATH}/$1~${SLURM_JOB_ID}~$timestamp.log"
ln "${SLURM_SUBMIT_DIR}/temp${SLURM_JOB_ID}.err" \
"${LOG_FILE_PATH}/$1~${SLURM_JOB_ID}~$timestamp.err"

## ========================= ##
##    VARIABLE ASSIGNMENT    ##
## ========================= ##

BLUEPRINT_MARK_NAME=$1

## ================== ##
##    MOVING FILES    ##
## ================== ##

cd "${BLUEPRINT_MAIN_DIR}" || \
{ echo "Main directory doesn't exist, \
make sure config.txt is pointing to the correct directory"; exit 1; }

# Find all of the bam files with the mark name given by the user
List_Of_Files_With_Mark_Name=\
$(find . -type f -name "*${BLUEPRINT_MARK_NAME}*.bam")

# If the above list is empty then the epigenetic mark must not exist.
Number_Of_Files_To_Move=\
$(find . -type f -name "*${BLUEPRINT_MARK_NAME}*.bam" | wc -l)
echo "Number of .bam files to be moved is: ${Number_Of_Files_To_Move}"

# Exit here if there are no files 
# There is no point in generating a folder if the mark doesn't exist
if [ "${Number_Of_Files_To_Move}" -eq 0 ]; then
    echo "No files were found."
    echo "Please input a epigenetic mark name that exists (case sensitive)." 
    echo "Aborting..."

    # Remove temporary log files
    rm "${SLURM_SUBMIT_DIR}/temp${SLURM_JOB_ID}.log"
    rm "${SLURM_SUBMIT_DIR}/temp${SLURM_JOB_ID}.err"

    exit 1
fi

mkdir -p "${RAW_DIR}/${BLUEPRINT_MARK_NAME}"

echo "Moving bam files to ${RAW_DIR}/${BLUEPRINT_MARK_NAME}"

for file in ${List_Of_Files_With_Mark_Name}; do
    echo "Moving file ${file}..."
    mv "${file}" "${RAW_DIR}/${BLUEPRINT_MARK_NAME}"
done


## ----------------------- ##
##   LOG FILE MANAGEMENT   ##
## ----------------------- ##

#Finish message and time
echo "Job completed at:"
date -u
end_time=$(date +%s)
time_taken=$((end_time-start_time))
echo "Job took a total of: ${time_taken} seconds to complete"

# Remove temporary log files
rm "${SLURM_SUBMIT_DIR}/temp${SLURM_JOB_ID}.log"
rm "${SLURM_SUBMIT_DIR}/temp${SLURM_JOB_ID}.err"