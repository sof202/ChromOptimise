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
# Moving files on linux takes very little memory unless moving between drives
#SBATCH --mem=1G
# Send an email after the job is done
#SBATCH --mail-type=END
# Temporary log file, later to be removed
#SBATCH --output=temp%j.log
# Temporary error file, later to be removed
#SBATCH --error=temp%j.err
#SBATCH --job-name=1_Moving_Files

## =================================================================================##
##                                                                                  ||
##                                     PREAMBLE                                     ||
##                                                                                  ||
## =================================================================================##
## PURPOSE:                                                                         ||
## Moves .bam files that include the inputted epigenetic mark into a single folder. ||
## =================================================================================##
## AUTHOR: Sam Fletcher                                                             ||
## CONTACT: s.o.fletcher@exeter.ac.uk                                               ||
## CREATED: November 2023                                                           ||
## =================================================================================##
## PREREQUISITES: Downloaded files from EGA                                         ||
## =================================================================================##
## DEPENDENCIES: NONE                                                               ||
## =================================================================================##
## INPUTS:                                                                          ||
## $1 -> Full (or relative) file path for configuation file directory               ||
## $2 -> Epigenetic Mark to process                                                 ||
## =================================================================================##
## OUTPUTS:                                                                         ||
## NONE                                                                             ||
## =================================================================================##

## ======================== ##
##    HELP FUNCTIONALITY    ##
## ======================== ##

if [[ "$1" == "--help" || "$1" == "-h" ]]; then
    echo "===================================================================="
    echo "Purpose: Moves .bam files that include the epigenetic mark."
    echo "into a single folder."
    echo "Author: Sam Fletcher"
    echo "Contact: s.o.fletcher@exeter.ac.uk"
    echo "Dependencies: NONE"
    echo "Inputs:"
    echo "\$1 -> Full (or relative) file path for configuation file directory"
    echo "\$2 -> Epigenetic mark name"
    echo "===================================================================="
    exit 0
fi

## ============ ##
##    SET UP    ##
## ============ ##

# Configuration files are required for file paths and log file management
configuration_directory=$1

source "${configuration_directory}/FilePaths.txt" || \
{ echo "The configuration file does not exist in the specified location: \
${configuration_directory}"; exit 1; }

# If a configuration file is changed during analysis, it is hard to tell
# what configuration was used for a specific run through, below accounts for 
# this
echo "Configuration file used with this script: \
${configuration_directory}/FilePaths.txt"
echo ""
cat "${configuration_directory}/FilePaths.txt"
echo ""

source "${configuration_directory}/LogFileManagement.sh" || \
{ echo "The log file management script does not exist in the specified \
location: ${configuration_directory}"; exit 1; }



# Output and error files renamed to:
# [epigenetic mark name]~[job id]~[date]-[time]

mv "${SLURM_SUBMIT_DIR}/temp${SLURM_JOB_ID}.log" \
"${LOG_FILE_PATH}/$2~${SLURM_JOB_ID}~${timestamp:=}.log"
mv "${SLURM_SUBMIT_DIR}/temp${SLURM_JOB_ID}.err" \
"${LOG_FILE_PATH}/$2~${SLURM_JOB_ID}~$timestamp.err"

## =============== ##
##    VARIABLES    ##
## =============== ##

mark_name=$2

## ================== ##
##    MOVING FILES    ##
## ================== ##

cd "${MAIN_DIR}" || \
{ >&2 echo "ERROR: [\${MAIN_DIR} - ${MAIN_DIR}] doesn't exist, \
make sure FilePaths.txt is pointing to the correct directory"; finishing_statement 1; }

list_of_files_with_mark_name=\
$(find . -type f -name "*${mark_name}*.bam")

# If the above list is empty then the epigenetic mark must not exist.
if [[ -z ${list_of_files_with_mark_name} ]]; then
    { >&2 echo -e "ERROR: No files with epigenetic mark: ${mark_name} were found.\n"\
    "Please input a epigenetic mark name that exists (note that this is "\
    "case sensitive)."; finishing_statement 1; }
fi

mkdir -p "${RAW_DIR}/${mark_name}"

echo "Moving .bam files to ${RAW_DIR}/${mark_name}..."

for file in ${list_of_files_with_mark_name}; do
  mv "${file}" "${RAW_DIR}/${mark_name}"
done

finishing_statement 0