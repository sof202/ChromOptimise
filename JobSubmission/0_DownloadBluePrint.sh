#!/bin/bash
#SBATCH --export=ALL #export all enviroment variables to the batch job
#SBATCH -p mrcq #submit to the mrcq for faster timing
#SBATCH --time=150:00:00 ##maximum wall time for the job
#SBATCH -A Research_Project-MRC190311 #research project to submit under
#SBATCH --nodes=1 #specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --output=temp%j.log #Put output file in log files with temporary name
#SBATCH --error=temp%j.err #Put error file in log files with temporary name
#SBATCH --job-name=Download_Files

## -------------------------------------------------------------------------------------------- ##
##                                                                                              ##
##                                            PREAMBLE                                          ##
##                                                                                              ##
## -------------------------------------------------------------------------------------------- ##
##                                            PURPOSE                                           ##
##                                    Downloads files from EGA.                                 ##
## -------------------------------------------------------------------------------------------- ##
##                        AUTHOR: Sam Fletcher s.o.fletcher@exeter.ac.uk                        ##
##                                     CREATED: November 2023                                   ##
## -------------------------------------------------------------------------------------------- ##
##                                         PREREQUISITES                                        ##
##                 Created a conda environement that has pyega installed in it                  ##
##              Changed souce file in |SET UP| to be your personal conda.sh file                ##
##                  Created .json file containing your EGA login credentials                    ##
## -------------------------------------------------------------------------------------------- ##
##                                          DEPENDENCIES                                        ##
##                                   Miniconda/Conda/Anaconda                                   ##
## -------------------------------------------------------------------------------------------- ##
##                                             INPUTS                                           ##
##                      $1 -> File of Names of files/directories to download                    ##
## -------------------------------------------------------------------------------------------- ##
##                                            OUTPUTS                                           ##
##                      Files that were chosen to be downloaded from EGA                        ##
## -------------------------------------------------------------------------------------------- ##

## ------------------------ ##
##    HELP FUNCTIONALITY    ##
## ------------------------ ##

if [[ "$1" == "--help" || "$1" == "-h" ]]; then
    echo "======================================================================"
    echo "Purpose: Downloads files from EGA from a list of file/directory names."
    echo "Author: Sam Fletcher"
    echo "Contact: s.o.fletcher@exeter.ac.uk"
    echo "Dependencies: Miniconda/Conda/Anaconda, EGA account details"
    echo "Inputs:"
    echo "\$1 -> File of names of files/directories to download from EGA"
    echo "======================================================================"
    exit 0
fi

## ------------ ##
##    SET UP    ##
## ------------ ##

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
# [job id]~[date]-[time]
# This requires a hard link as you cannot rename log files
# whilst running the script without a wrapper function
LOG_FILE_PATH="${LOG_DIR}/$USER/$SLURM_JOB_NAME/"
timestamp=$(date -u +%Y.%m.%d-%H:%M)
ln "${SLURM_SUBMIT_DIR}/temp${SLURM_JOB_ID}.log" \
"${LOG_FILE_PATH}/${SLURM_JOB_ID}~$timestamp.log"
ln "${SLURM_SUBMIT_DIR}/temp${SLURM_JOB_ID}.err" \
"${LOG_FILE_PATH}/${SLURM_JOB_ID}~$timestamp.err"


## ---------- ##
##    MAIN    ##
## ---------- ##


#set up conda environment
module purge
module load Miniconda3
source /lustre/home/sof202/miniconda3/etc/profile.d/conda.sh 
# Conda environments will not be avtivated until one uses `conda init bash`
# However, running this will result in a new shell being created.
# This means one cannot have their environment activatable and activate it
# Using the conda shell script is a known work around for this
conda activate pyega
# Name of conda environment may be different depending on system

Text_File_Containing_Inodes=$1

# Looping over each line in the Text file that contains the files/directories to be downloaded...
# We use pyega3 to download these files that are being read in.

# Currently configured to download into a specific directory, change this if you want a different one.
while IFS= read -r line; do
    # mkdir -p "${DOWNLOAD_DIR}"
    #pyega3 -c 5 -cf ~/pyegaDownloading/egaConfig.json fetch "$line" --output-dir "${DOWNLOAD_DIR}" 

    mkdir -p "/lustre/projects/Research_Project-MRC190311/blueprint/EGAD00001002670/$line"
    pyega3 -c 5 -cf ~/pyegaDownloading/egaConfig.json fetch "$line" --output-dir "/lustre/projects/Research_Project-MRC190311/blueprint/EGAD00001002670/" 
    #Location of config file (containing login credentials) will likely differ for different users
done < "${Text_File_Containing_Inodes}"


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