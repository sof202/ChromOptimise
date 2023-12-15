#!/bin/bash
# export all enviroment variables to the batch job
#SBATCH --export=ALL
# submit to the mrcq for faster queue times
#SBATCH -p mrcq
# Downloading can take a long time if files are large/abundant
#SBATCH --time=150:00:00
#SBATCH -A Research_Project-MRC190311
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
# Send email at job completion
#SBATCH --mail-type=END
# Temporary log file, later to be removed
#SBATCH --output=temp%j.log
# Temporary error file, later to be removed
#SBATCH --error=temp%j.err
#SBATCH --job-name=Download_Files

## =================================================================================##
##                                                                                  ||
##                                     PREAMBLE                                     ||
##                                                                                  ||
## =================================================================================##
## PURPOSE:                                                                         ||
## Downloads files from EGA.                                                        ||
## =================================================================================##
## AUTHOR: Sam Fletcher s.o.fletcher@exeter.ac.uk                                   ||
## CREATED: November 2023                                                           ||
## =================================================================================##
## PREREQUISITES:                                                                   ||
## Create a conda environement that has pyega installed in it                       ||
## Change source file in |MAIN| to be your personal etc/profile.d/conda.sh file     ||
## Create a .json file containing your EGA login credentials                        ||
## =================================================================================##
## DEPENDENCIES:                                                                    ||
## Miniconda/Conda/Anaconda                                                         ||
## Python                                                                           ||
## Pyega3                                                                           ||
## =================================================================================##
## INPUTS:                                                                          ||
## $1 -> File of Names of files/directories to download                             ||
## =================================================================================##
## OUTPUTS:                                                                         ||
## NONE                                                                             ||
## =================================================================================##

## ======================== ##
##    HELP FUNCTIONALITY    ##
## ======================== ##

if [[ "$1" == "--help" || "$1" == "-h" ]]; then
    echo "======================================================================"
    echo "Purpose: Downloads files from EGA from a list of file/directory names."
    echo "Author: Sam Fletcher"
    echo "Contact: s.o.fletcher@exeter.ac.uk"
    echo "Dependencies: Miniconda/Conda/Anaconda, EGA account details, Python"
    echo "Inputs:"
    echo "\$1 -> File of names of files/directories to download from EGA"
    echo "======================================================================"
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
# CHANGE THIS TO YOUR OWN CONFIG FILE
echo "Loading config file..."
source "/lustre/projects/Research_Project-MRC190311\
/scripts/integrative/blueprint/config/config.txt"

# Rename the output and error files to have format:
# [job id]~[date]-[time]
# This requires a hard link as you cannot rename log files
# whilst running the script without a wrapper function
LOG_FILE_PATH="${LOG_DIR}/$USER/$SLURM_JOB_NAME/"
timestamp=$(date -u +%Y.%m.%d-%H:%M)
ln "${SLURM_SUBMIT_DIR}/temp${SLURM_JOB_ID}.log" \
"${LOG_FILE_PATH}/${SLURM_JOB_ID}~$timestamp.log"
ln "${SLURM_SUBMIT_DIR}/temp${SLURM_JOB_ID}.err" \
"${LOG_FILE_PATH}/${SLURM_JOB_ID}~$timestamp.err"

## ========================= ##
##    VARIABLE ASSIGNMENT    ##
## ========================= ##

Text_File_Containing_Inodes=$1

## ========== ##
##    MAIN    ##
## ========== ##

module purge
module load Miniconda3

# Conda environments will not be activated until one uses `conda init bash`
# However, running this will result in a new shell being created.
# This means one cannot have their environment activatable and activate it
# Using the conda shell script is a known work around for this.
source /lustre/home/sof202/miniconda3/etc/profile.d/conda.sh 

# CHANGE THIS TO YOUR CONDA ENVIRONMENT NAME
conda activate pyega

# Currently configured to download into a specific directory
# change this if you want a different one.
while IFS= read -r line; do
    # mkdir -p "${DOWNLOAD_DIR}"
    # pyega3 -c 5 -cf ~/pyegaDownloading/egaConfig.json fetch \
    # "$line" --output-dir "${DOWNLOAD_DIR}" 

    mkdir -p \
    "/lustre/projects/Research_Project-MRC190311/blueprint/EGAD00001002670/$line"
    # CHANGE egaConfig.json TO FILE WITH EGA LOGIN CREDENTIALS
    pyega3 -c 5 -cf ~/pyegaDownloading/egaConfig.json fetch "$line" --output-dir \
    "/lustre/projects/Research_Project-MRC190311/blueprint/EGAD00001002670/" 
done < "${Text_File_Containing_Inodes}"


## ======================= ##
##   LOG FILE MANAGEMENT   ##
## ======================= ##


#Finish message and time
echo "Job completed at:"
date -u
end_time=$(date +%s)
time_taken=$((end_time-start_time))
echo "Job took a total of: ${time_taken} seconds to complete."

# Remove temporary log files
rm "${SLURM_SUBMIT_DIR}/temp${SLURM_JOB_ID}.log"
rm "${SLURM_SUBMIT_DIR}/temp${SLURM_JOB_ID}.err"