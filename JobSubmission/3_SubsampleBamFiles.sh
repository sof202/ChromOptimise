#!/bin/bash
# Export all environment variables to the batch job
#SBATCH --export=ALL
# Submit to the mrc queue for faster queue times
#SBATCH -p mrcq
# Change this depending on the number of files that are to be merged. 
# In previous tests, merging 2 files of 1-2 Gigabytes each took ~7 minutes.
#SBATCH --time=01:00:00 
#SBATCH -A Research_Project-MRC190311 
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=16 
# Lots of memory is required for this script as the merging process 
# opens up every file at the same time piece by piece
#SBATCH --mem=100G 
# Send an email after the job is done
#SBATCH --mail-type=END 
# Temporary log file, later to be removed
#SBATCH --output=temp%j.log
# Temporary error file, later to be removed
#SBATCH --error=temp%j.err
#SBATCH --job-name=Merging_and_Subsampling

## =================================================================================##
##                                                                                  ||
##                                     PREAMBLE                                     ||
##                                                                                  ||
## =================================================================================##
## PURPOSE:                                                                         ||
## Obtain a sample of the blueprint bam files. However, the bam files have varying  ||
## sizes due to the number of reads. Sampling the files produced in                 ||
## 2_ProcessBamFiles.sh will lead to samples with the same number of files that     ||
## contain a different number of reads. This removes the reproducability of the     ||
## proceedure. To get around this, this script merges all of the processed          ||
## blueprint .bam files and subsequently samples this larger file randomly.         ||
## To save space, the merged file is deleted.                                       ||
## =================================================================================##
## AUTHOR: Sam Fletcher s.o.fletcher@exeter.ac.uk                                   ||
## CREATED: November 2023                                                           ||
## =================================================================================##
## PREREQUISITES: Run 2_batch_ProcessBamFiles.sh                                    ||
## =================================================================================##
## DEPENDENCIES: Samtools                                                           ||
## =================================================================================##
## INPUTS:                                                                          ||
## $1 -> Epigenetic mark to process                                                 ||
## $2 -> Sample size as a percentage                                                ||
## =================================================================================##
## OUTPUTS:                                                                         ||
## Subsampled .bam files                                                            ||
## =================================================================================##

## ======================== ##
##    HELP FUNCTIONALITY    ##
## ======================== ##

if [[ "$1" == "--help" || "$1" == "-h" ]]; then
    echo "==================================================="
    echo "Purpose: Merges and subsamples processed .bam files"
    echo "present in specified folder"
    echo "Author: Sam Fletcher"
    echo "Contact: s.o.fletcher@exeter.ac.uk"
    echo "Dependencies: Samtools"
    echo "Inputs:"
    echo "\$1 -> Name of epigenetic mark"
    echo "\$2 -> Sample size as a percentage"
    echo "==================================================="
    exit 0
fi

## ============ ##
##    SET UP    ##
## ============ ##

echo "Job '$SLURM_JOB_NAME' started at:"
date -u

start_time=$(date +%s)

# Activate config.txt to access all file paths
# CHANGE THIS TO YOUR OWN CONFIG FILE
echo "Loading config file..."
source "/lustre/projects/Research_Project-MRC190311\
/scripts/integrative/blueprint/config/config.txt"

# Rename the output and error files to have format:
# [epigenetic mark name]~[Sample size]~[job id]~[date]-[time]
# This requires a hard link as you cannot rename log files
# whilst running the script without a wrapper function
LOG_FILE_PATH="${LOG_DIR}/$SLURM_JOB_NAME/$USER"
mkdir -p "${LOG_FILE_PATH}"
timestamp=$(date -u +%Y.%m.%d-%H:%M)

ln "${SLURM_SUBMIT_DIR}/temp${SLURM_JOB_ID}.log" \
"${LOG_FILE_PATH}/$1~$2~${SLURM_JOB_ID}~$timestamp.log"
ln "${SLURM_SUBMIT_DIR}/temp${SLURM_JOB_ID}.err" \
"${LOG_FILE_PATH}/$1~$2~${SLURM_JOB_ID}~$timestamp.err"

## ========================= ##
##    VARIABLE ASSIGNMENT    ##
## ========================= ##

BLUEPRINT_MARK_NAME=$1
SAMPLE_SIZE=$2
BLUEPRINT_PROCESSED_FILE_PATH="${PROCESSED_DIR}/${BLUEPRINT_MARK_NAME}"

if [ -z "${BLUEPRINT_MARK_NAME}" ]; then
    echo "No Blueprint epigenetic mark name given."
    echo "Ensure first argument is the name of the epigenetic mark."
    echo "Aborting..."

    rm "${SLURM_SUBMIT_DIR}/temp${SLURM_JOB_ID}.log"
    rm "${SLURM_SUBMIT_DIR}/temp${SLURM_JOB_ID}.err"
    exit 1
fi

if [ -z "${SAMPLE_SIZE}" ]; then
    SAMPLE_SIZE=50
    echo "No sample size was given, using default value of 50 percent."
fi
echo -n "Subsampling processed .bam files using sample size of: "
echo "${SAMPLE_SIZE} percent for epigenetic mark: ${BLUEPRINT_MARK_NAME}"

## ========================= ##
##   MERGING OF .BAM FILES   ##
## ========================= ##

if [ -d "${BLUEPRINT_PROCESSED_FILE_PATH}" ]; then
    echo "Changing directory to: ${BLUEPRINT_PROCESSED_FILE_PATH}" 
    cd "${BLUEPRINT_PROCESSED_FILE_PATH}" || exit 1
else
    echo "Directory does not exist yet." 
    echo -n "Make sure you typed the epigenetic mark correctly "
    echo "and that you have ran 2_ProcessBamFiles.sh first."
    echo "Aborting..."

    rm "${SLURM_SUBMIT_DIR}/temp${SLURM_JOB_ID}.log"
    rm "${SLURM_SUBMIT_DIR}/temp${SLURM_JOB_ID}.err"
    exit 1
fi

echo "Finding suitable .bam files to merge..."
find . -type f -name "*.sorted.filtered.noDuplicates.bam" \
> List_Of_Bam_Files_To_Merge.txt

module purge
module load SAMtools
output_file_path="${SUBSAMPLED_DIR}/FullMerged.${BLUEPRINT_MARK_NAME}.bam"

echo "Merging..."
samtools merge -b List_Of_Bam_Files_To_Merge.txt "${output_file_path}"

## ===================================== ##
##    SUBSAMPLING OF MERGED .BAM FILE    ##
## ===================================== ##

cd "${SUBSAMPLED_DIR}" || { echo "Subsampled directory doesn't exist, \
make sure config.txt is pointing to the correct directory"; exit 1; }

sample_size_decimal=$(echo "scale=2; $SAMPLE_SIZE /100" | bc)
echo "Subsampling..."
# Ensure headers are kept in subsampled file to avoid errors later in pipeline
samtools view -H "${output_file_path}" \
> "Subsampled.${SAMPLE_SIZE}.${BLUEPRINT_MARK_NAME}.bam"
samtools view -s "${sample_size_decimal}" "${output_file_path}" \
>> "Subsampled.${SAMPLE_SIZE}.${BLUEPRINT_MARK_NAME}.bam"

rm "FullMerged.${BLUEPRINT_MARK_NAME}.bam"
cd "${BLUEPRINT_PROCESSED_FILE_PATH}" || exit 1
rm List_Of_Bam_Files_To_Merge.txt

## ======================= ##
##   LOG FILE MANAGEMENT   ##
## ======================= ##

echo "Job completed at:"
date -u
end_time=$(date +%s)
time_taken=$((end_time-start_time))
echo "Job took a total of: ${time_taken} seconds to complete"

rm "${SLURM_SUBMIT_DIR}/temp${SLURM_JOB_ID}.log"
rm "${SLURM_SUBMIT_DIR}/temp${SLURM_JOB_ID}.err"