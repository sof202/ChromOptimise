#!/bin/bash
# Export all environment variables to the batch job
#SBATCH --export=ALL
# Submit to the mrc queue for faster queue times
#SBATCH -p mrcq 
# Previous tests:(50%) 3 and 4 GB -> 4 mins. 5 and 7 GB -> 7 mins, 6 and 7 GB -> 7 mins
#SBATCH --time=01:00:00 
#SBATCH -A Research_Project-MRC190311 
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=16 
# Predicted that memory consumption will rise massively when merging lots of files
# For 2 files, the merging only uses 5 MB. Though there is a large number of temporary
# allocations
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
## Obtain a sample of the bam files. However, the bam files have varying sizes due  ||
## to the number of reads. Sampling the files produced in 2_ProcessBamFiles.sh will ||
## lead to samples with the same number of files that contain a different number    ||
## of reads. This removes the reproducability of the proceedure. To get around      ||
## this, this script merges all of the processed .bam files and subsequently        ||
## samples this larger file randomly. To save space, the merged file is deleted.    ||
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
## $2 -> Sample size as a percentage (default : 50)                                 ||
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
    echo "\$2 -> Sample size as a percentage (default : 50)"
    echo "==================================================="
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

## ============================= ##
##    VARIABLES AND FUNCTIONS    ##
## ============================= ##

mark_name=$1
sample_size=$2
PROCESSED_FILE_PATH="${PROCESSED_DIR}/${mark_name}"

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

## ====== DEFAULTS ====================================================================
if [[ -z "${mark_name}" ]]; then
    { >&2 echo -e "ERROR: No epigenetic mark name given.\n\
    Ensure that the first argument is the name of a processed epigenetic mark." ;}

    finishing_statement 1
fi

if [[ -z "${sample_size}" || "${sample_size}" =~ ^[^0-9]+$ ]]; then
    sample_size=50
    echo "Invalid sample size was given, using default value of: ${sample_size}%."
fi

if [[ "${sample_size}" -gt 100 || "${sample_size}" -le 0 ]]; then
    echo "Sample size must be greater than 0 and less than or equal to 100."
    sample_size=50
    echo "Using a default value of: ${sample_size}% instead."
fi
## ====================================================================================

## ========================= ##
##   MERGING OF .BAM FILES   ##
## ========================= ##

echo -n "Merging processed .bam files for epigenetic mark: ${mark_name}"

cd "${PROCESSED_FILE_PATH}" || \
{ >&2 echo "ERROR: \${PROCESSED_DIR}/\${mark_name} - ${PROCESSED_DIR}/${mark_name} \
doesn't exist, make sure that you typed the epigenetic mark correctly and that \
config.txt is pointing to the correct directory"; finishing_statement 1; }

find . -type f -name "*.sorted.filtered.noDuplicates.bam" \
> List_Of_Bam_Files_To_Merge.txt

module purge
module load SAMtools
output_file_path="${SUBSAMPLED_DIR}/FullMerged.${mark_name}.bam"

echo "Merging the following files:"
cat List_Of_Bam_Files_To_Merge.txt

samtools merge -b List_Of_Bam_Files_To_Merge.txt "${output_file_path}"

## ===================================== ##
##    SUBSAMPLING OF MERGED .BAM FILE    ##
## ===================================== ##

cd "${SUBSAMPLED_DIR}" || \
{ >&2 echo "ERROR: \${SUBSAMPLED_DIR} - ${SUBSAMPLED_DIR} doesn't exist, \
make sure config.txt is pointing to the correct directory"; finishing_statement 1; }

sample_size_decimal=$(echo "scale=2; $sample_size /100" | bc)
echo "Subsampling merged .bam file with sample size ${sample_size}%..."
# Ensure headers are kept in subsampled file to avoid errors later in pipeline
samtools view -H "${output_file_path}" \
> "Subsampled.${sample_size}.${mark_name}.bam"
samtools view -s "${sample_size_decimal}" "${output_file_path}" \
>> "Subsampled.${sample_size}.${mark_name}.bam"

rm "FullMerged.${mark_name}.bam"
cd "${PROCESSED_FILE_PATH}" || finishing_statement 1
rm List_Of_Bam_Files_To_Merge.txt

finishing_statement 0