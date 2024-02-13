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
#SBATCH --job-name=0_Download_Files

## =================================================================================##
##                                                                                  ||
##                                     PREAMBLE                                     ||
##                                                                                  ||
## =================================================================================##
## PURPOSE:                                                                         ||
## Downloads files from EGA using pyega3.                                           ||
## =================================================================================##
## AUTHOR: Sam Fletcher                                                             ||
## CONTACT: s.o.fletcher@exeter.ac.uk                                               ||
## CREATED: November 2023                                                           ||
## =================================================================================##
## PREREQUISITES:                                                                   ||
## Create a conda environement that has pyega3 installed in it                      ||
## Create a .json file containing your EGA login credentials                        ||
## =================================================================================##
## DEPENDENCIES:                                                                    ||
## Miniconda/Conda/Anaconda                                                         ||
## Pyega3 conda environment                                                         ||
## =================================================================================##
## INPUTS:                                                                          ||
## -c|--config -> Full/relative file path for configuation file directory           ||
## -f|--file -> File of file names to download from EGA.                            ||
## =================================================================================##
## OUTPUTS:                                                                         ||
## NONE                                                                             ||
## =================================================================================##

## ===================== ##
##   ARGUMENT PARSING    ##
## ===================== ##

usage() {
cat <<EOF
=======================================================================
0_EGADownloading
=======================================================================
Purpose: Downloads files from EGA using a list of file/directory names.
Author: Sam Fletcher
Contact: s.o.fletcher@exeter.ac.uk
Dependencies: Conda, EGA login credentials, Pyega3 conda environment
Inputs:
-c|--config -> Full/relative file path for configuation file directory
-f|--file -> File of file names to download from EGA.
=======================================================================
EOF
    exit 0
}

if [[ $# -eq 0 ]]; then
    usage
fi

while [[ $# -gt 0 ]]; do
    case $1 in
        -c|--config)
        if [[ "$2" != -* && -n "$2" ]]; then 
            configuration_directory="$2"; shift 2
        else 
            shift 1
        fi
        ;;
        -f|--file)
        if [[ "$2" != -* && -n "$2" ]]; then 
            text_file_containing_inodes="$2"; shift 2
        else 
            shift 1
        fi
        ;;
        *)
        usage
    esac
done

## ============ ##
##    SET UP    ##
## ============ ##

echo "Job '${SLURM_JOB_NAME}' started at:"
date -u

start_time=$(date +%s)

source "${configuration_directory}/FilePaths.txt" || \
{ echo "The configuration file does not exist in the specified location: \
${configuration_directory}"; exit 1; }


LOG_FILE_PATH="${LOG_DIR}/$SLURM_JOB_NAME/$USER"
mkdir -p "${LOG_FILE_PATH}"
timestamp=$(date -u +%Y.%m.%d-%H_%M)

# Output and error files renamed to:
# [file name]~[job id]~[date]-[time]

mv "${SLURM_SUBMIT_DIR}/temp${SLURM_JOB_ID}.log" \
"${LOG_FILE_PATH}/${SLURM_JOB_ID}~$timestamp.log"
mv "${SLURM_SUBMIT_DIR}/temp${SLURM_JOB_ID}.err" \
"${LOG_FILE_PATH}/${SLURM_JOB_ID}~$timestamp.err"

## ========== ##
##    MAIN    ##
## ========== ##

module purge
module load Miniconda3

# Conda environments will not be activated until one uses `conda init bash`
# However, running this will result in a new shell being created.
# This means one cannot have their environment activatable and activate it
# Using the conda shell script in the [conda]/etc folder is a work around for this.
source "${CONDA_SHELL}/profile.d/conda.sh" || \
{ echo "profile.d/conda.sh does not exist in specified location: \
[\${CONDA_SHELL} - ${CONDA_SHELL}]"; exit 1; }

conda activate "${PYEGA_ENVIRONMENT}" || \
{ echo "conda environment does not exist in specified location: \
[\${PYEGA_ENVIRONMENT} - ${PYEGA_ENVIRONMENT}]"; exit 1; }

# Read each line of text file
# [[ -n "$line" ]] handles the last line that has no newline character
while IFS= read -r line || [[ -n "$line" ]]; do
    # CHANGE "egaConfig.json" TO FILE WITH EGA LOGIN CREDENTIALS
    # -c 5 -> Failed downloads are retried 5 times before moving on to next file
    pyega3 -c 5 -cf ~/Tools/pyegaDownloading/egaConfig.json fetch \
    "$line" --output-dir "${DOWNLOAD_DIR}"
done < "${text_file_containing_inodes}"


## ======================= ##
##   FINISHING STATEMENT   ##
## ======================= ##

rm "${SLURM_SUBMIT_DIR}/pyega3_output.log"
echo "Job finished with exit code 0 at:"
date -u
end_time=$(date +%s)
time_taken=$((end_time-start_time))
echo "Job took a total of: ${time_taken} seconds to finish."
