#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -p mrcq # submit to the mrc queue for faster queue times
#SBATCH --time=01:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC190311 
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=16 
#SBATCH --mem=10G # specify bytes memory to reserve
#SBATCH --mail-type=END # Send an email after the job is done
#SBATCH --output=/lustre/projects/Research_Project-MRC190311/scripts/integrative/blueprint/LogFiles/Moving_Files/temp.o #Put output file in log files with temporary name
#SBATCH --error=/lustre/projects/Research_Project-MRC190311/scripts/integrative/blueprint/LogFiles/Moving_Files/temp.e #Put error file in log files with temporary name
#SBATCH --job-name=Moving_Files

## -------------------------------------------------------------------------------------------- ##
##                                                                                              ##
##                                            PREAMBLE                                          ##
##                                                                                              ##
## -------------------------------------------------------------------------------------------- ##
##                                            PURPOSE                                           ##
##           Moves .bam files that include the epigenetic mark into a single folder.            ##
## -------------------------------------------------------------------------------------------- ##
##                        AUTHOR: Sam Fletcher s.o.fletcher@exeter.ac.uk                        ##
##                                     CREATED: November 2023                                   ##
## -------------------------------------------------------------------------------------------- ##
##                                         PREREQUISITES                                        ##
##                                   Downloaded files from EGA                                  ##
## -------------------------------------------------------------------------------------------- ##
##                                          DEPENDENCIES                                        ##
##                                              NONE                                            ##
## -------------------------------------------------------------------------------------------- ##
##                                             INPUTS                                           ##
##                           $1 -> Blueprint Epigenetic Mark to process                         ##
## -------------------------------------------------------------------------------------------- ##
##                                            OUTPUTS                                           ##
##                                              NONE                                            ##
## -------------------------------------------------------------------------------------------- ##

## ------------------------ ##
##    HELP FUNCTIONALITY    ##
## ------------------------ ##

if [[ "$1" == "--help" || "$1" == "-h" ]]; then
    echo "==========================================================================================="
    echo "Purpose: Moves .bam files that include the epigenetic mark into a single folder."
    echo "Author: Sam Fletcher"
    echo "Contact: s.o.fletcher@exeter.ac.uk"
    echo "Dependencies: NONE"
    echo "Inputs:"
    echo "\$1 -> Epigenetic mark name"
    echo "==========================================================================================="
    exit 0
fi

## ------------ ##
##    SET UP    ##
## ------------ ##

# Rename the output and error files to have format: [epigenetic mark name]~[job id]~([array id])~[date]-[time]
# This requires a hard link as you cannot rename log files whilst running the script without a wrapper function
LOG_FILE_PATH=/lustre/projects/Research_Project-MRC190311/scripts/integrative/blueprint/LogFiles/$SLURM_JOB_NAME/
cd "${LOG_FILE_PATH}" || exit 1
timestamp=$(date -u +%Y.%m.%d-%H:%M)
ln temp.e "$1~${SLURM_JOB_ID}~$timestamp.e"
ln temp.o "$1~${SLURM_JOB_ID}~$timestamp.o"

# Print start date/time
echo "Job '$SLURM_JOB_NAME' started at:"
date -u

# Get the start time for the program
start_time=$(date +%s)

# Activate config.txt to access all file paths
echo "Loading config file..."
source /lustre/projects/Research_Project-MRC190311/scripts/integrative/blueprint/config/config.txt


# Open the main blueprint directory
BLUEPRINT_MARK_NAME=$1
echo "Changing directory to: ${BLUEPRINT_MAIN_DIR}"
cd "${BLUEPRINT_MAIN_DIR}" || { echo "Main directory doesn't exist, make sure config.txt is pointing to the correct directory"; exit 1; }


## ------------------ ##
##    MOVING FILES    ##
## ------------------ ##


# Find all of the bam files with the mark name given by the user
List_Of_Files_With_Mark_Name=$(find . -type f -name "*${BLUEPRINT_MARK_NAME}*.bam")

# If the above list is empty then the epigenetic mark must not exist. Print this as an error
Number_Of_Files_To_Move=$(find . -type f -name "*${BLUEPRINT_MARK_NAME}*.bam" | wc -l)
echo "Number of .bam files to be moved is: ${Number_Of_Files_To_Move}"

# Exit here if there are no files (There is no point in generating a folder if the mark doesn't exist)
if [ "${Number_Of_Files_To_Move}" -eq 0 ]; then
    echo "No files were found, please input a epigenetic mark name that exists (case sensitive)." 
    echo "Aborting..."

    # Remove temporary log files
    cd "${LOG_FILE_PATH}" || exit 1
    rm temp.e
    rm temp.o

    exit 1
fi

mkdir -p "${RAW_DIR}/${BLUEPRINT_MARK_NAME}"

# Use a for loop to move all of the found files to the directory that holds each individual mark
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
cd "${LOG_FILE_PATH}" || exit
rm temp.e
rm temp.o