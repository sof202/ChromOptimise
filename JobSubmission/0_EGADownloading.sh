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

## ===========================================================================##
##                                                                            ||
##                                  PREAMBLE                                  ||
##                                                                            ||
## ===========================================================================##
## PURPOSE:                                                                   ||
## Downloads files from EGA using pyega3.                                     ||
## ===========================================================================##
## AUTHOR: Sam Fletcher                                                       ||
## CONTACT: s.o.fletcher@exeter.ac.uk                                         ||
## CREATED: November 2023                                                     ||
## ===========================================================================##
## PREREQUISITES:                                                             ||
## Create a conda environement that has pyega3 installed in it                ||
## Create a .json file containing your EGA login credentials                  ||
## ===========================================================================##
## DEPENDENCIES:                                                              ||
## Miniconda/Conda/Anaconda                                                   ||
## Pyega3 conda environment                                                   ||
## ===========================================================================##
## INPUTS:                                                                    ||
## -c|--config= -> Full/relative file path for configuation file directory    ||
## -f|--file=   -> File of file names to download from EGA.                   ||
## ===========================================================================##
## OUTPUTS:                                                                   ||
## NONE                                                                       ||
## ===========================================================================##

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
-c|--config= -> Full/relative file path for configuation file directory
-f|--file=   -> File of file names to download from EGA.
=======================================================================
EOF
    exit 0
}

needs_argurment() {
    # Required check in case user uses -a -b or -b -a (no argument given).
    if [[ -z "$OPTARG" || "${OPTARG:0:1}" == - ]]; then usage; fi
}

while getopts f:c:-: OPT; do
  # Adds support for long options by reformulating OPT and OPTARG
  # This assumes that long options are in the form: "--long=option"
  if [ "$OPT" = "-" ]; then
    OPT="${OPTARG%%=*}"
    OPTARG="${OPTARG#"$OPT"}"
    OPTARG="${OPTARG#=}"
  fi
  case "$OPT" in
    c | config )  needs_argurment; configuration_directory="$OPTARG" ;;
    f | file )    needs_argurment; text_file_containing_inodes="$OPTARG" ;;
    \? )          usage ;;  # Illegal short options are caught by getopts
    * )           usage ;;  # bad long option
  esac
done
shift $((OPTIND-1))

## ============ ##
##    SET UP    ##
## ============ ##

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
# [file name]~[job id]~[date]-[time]

mv "${SLURM_SUBMIT_DIR}/temp${SLURM_JOB_ID}.log" \
"${LOG_FILE_PATH}/${SLURM_JOB_ID}~${timestamp:=}.log"
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
# Using the conda shell script in [conda]/etc is a work around for this.
source "${CONDA_SHELL}/profile.d/conda.sh" || \
{ echo "profile.d/conda.sh does not exist in specified location: \
[\${CONDA_SHELL} - ${CONDA_SHELL}]"; finishing_statement 1; }

conda activate "${PYEGA_ENVIRONMENT}" || \
{ echo "conda environment does not exist in specified location: \
[\${PYEGA_ENVIRONMENT} - ${PYEGA_ENVIRONMENT}]"; finishing_statement 1; }

# Read each line of text file
# [[ -n "$line" ]] handles the last line that has no newline character
while IFS= read -r line || [[ -n "$line" ]]; do
    # -c 5 -> Failed downloads are retried 5 times before moving on to next file
    pyega3 -c 5 -cf "${CREDENTIALS}" fetch \
    "$line" --output-dir "${DOWNLOAD_DIR}"
done < "${text_file_containing_inodes}"

rm "${SLURM_SUBMIT_DIR}/pyega3_output.log"
finishing_statement 1
