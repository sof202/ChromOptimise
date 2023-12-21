#!/bin/bash
# Export all enviroment variables to the batch job
#SBATCH --export=ALL
# Submit to the mrcq for faster queue times
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
    echo "======================================================================="
    echo "Purpose: Downloads files from EGA using a list of file/directory names."
    echo "Author: Sam Fletcher"
    echo "Contact: s.o.fletcher@exeter.ac.uk"
    echo "Dependencies: Miniconda/Conda/Anaconda, EGA login credentials, Python"
    echo "Inputs:"
    echo "\$1 -> File of file names to download from EGA."
    echo "======================================================================="
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
/scripts/integrative/blueprint/config/config.txt"

# Rename the output and error files to have format:
# [file name]~[job id]~[date]-[time]
# This requires a hard link as you cannot rename log files
# whilst running the script without a wrapper function
LOG_FILE_PATH="${LOG_DIR}/$SLURM_JOB_NAME/$USER"
mkdir -p "${LOG_FILE_PATH}"
timestamp=$(date -u +%Y.%m.%d-%H:%M)
ln "${SLURM_SUBMIT_DIR}/temp${SLURM_JOB_ID}.log" \
"${LOG_FILE_PATH}/File-$1~${SLURM_JOB_ID}~$timestamp.log"
ln "${SLURM_SUBMIT_DIR}/temp${SLURM_JOB_ID}.err" \
"${LOG_FILE_PATH}/File-$1~${SLURM_JOB_ID}~$timestamp.err"

## =============== ##
##    VARIABLES    ##
## =============== ##

text_file_containing_inodes=$1

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
conda activate /lustre/home/sof202/miniconda3/envs/pyega

mkdir -p "${DOWNLOAD_DIR}"
# Read each line and handle the last line without a newline character
while IFS= read -r line || [[ -n "$line" ]]; do
    # CHANGE egaConfig.json TO FILE WITH EGA LOGIN CREDENTIALS
    pyega3 -c 5 -cf ~/Tools/pyegaDownloading/egaConfig.json fetch \
    "$line" --output-dir "${DOWNLOAD_DIR}"
done < "${text_file_containing_inodes}"


## ======================= ##
##   LOG FILE MANAGEMENT   ##
## ======================= ##

echo "Job completed at:"
date -u
end_time=$(date +%s)
time_taken=$((end_time-start_time))
echo "Job took a total of: ${time_taken} seconds to complete."

rm "${SLURM_SUBMIT_DIR}/pyega3_output.log"
rm "${SLURM_SUBMIT_DIR}/temp${SLURM_JOB_ID}.log"
rm "${SLURM_SUBMIT_DIR}/temp${SLURM_JOB_ID}.err"