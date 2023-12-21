#!/bin/bash
# Export all environment variables to the batch job
#SBATCH --export=ALL
# Submit to the mrc queue for faster queue times
#SBATCH -p mrcq
# Tests have thus far shown a linear relationship between file size and time
# Current approximation is: [time (mins)] = 0.1 + 9*[file size (GB)]
#SBATCH --time=01:00:00 
#SBATCH -A Research_Project-MRC190311 
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=16
# Make sure array is not higher than the number of files being processed
# as this ends up with all files being processed by max index array
#SBATCH --array=1-4
# Peak memory consumption appears to come from samtools sort
# Previous tests show 1GB peak heap memory consumption with 1.5GB files
#SBATCH --mem=10G 
# Send an email after the job is done
#SBATCH --mail-type=END
# Temporary log file, later to be removed
#SBATCH --output=temp%A_%a.log
# Temporary error file, later to be removed
#SBATCH --error=temp%A_%a.err
#SBATCH --job-name=Processing

## =================================================================================##
##                                                                                  ||
##                                     PREAMBLE                                     ||
##                                                                                  ||
## =================================================================================##
## PURPOSE:                                                                         ||
## 1_RawBamFiles contains .bam files for each epigentic mark, however these have    ||
## yet to be processed, they still contain duplicates, multimapped fragments and    ||
## low quality reads. This script processes the bam files.                          ||
## =================================================================================##
## AUTHOR: Sam Fletcher s.o.fletcher@exeter.ac.uk                                   ||
## CREATED: November 2023                                                           ||
## =================================================================================##
## PREREQUISITES: All .bam files for a specific epigenetic mark must be in 1 folder ||
## =================================================================================##
## DEPENDENCIES: Samtools                                                           ||
## =================================================================================##
## INPUTS:                                                                          ||
## $1 -> Epigenetic mark to process                                                 ||
## $2 -> Phred score threshold value (default: 20)                                  ||
## =================================================================================##
## OUTPUTS:                                                                         ||
## Processed .bam files                                                             ||
## Index files for raw .bam files and processed .bam files                          ||
## Per chromosome stats and general stats for raw and processed .bam files          ||
## =================================================================================##

## ======================== ##
##    HELP FUNCTIONALITY    ##
## ======================== ##

if [[ "$1" == "--help" || "$1" == "-h" ]]; then
    echo "=============================================================="
    echo "Purpose: Processes .bam files by removing duplicates,"
    echo "filtering out poor quality reads and sorting."
    echo "Author: Sam Fletcher"
    echo "Contact: s.o.fletcher@exeter.ac.uk"
    echo "Dependencies: Samtools"
    echo "Inputs:"
    echo "\$1 -> Name of epigenetic mark"
    echo "\$2 -> Phred score threshold value (default: 20)"
    echo "Optional:"
    echo "Specify --array in sbatch options, to set a custom array size."
    echo "=============================================================="
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
# [epigenetic mark name]~[job id]~[array id]~[date]-[time]
# This requires a hard link as you cannot rename log files
# whilst running the script without a wrapper function
LOG_FILE_PATH="${LOG_DIR}/$SLURM_JOB_NAME/$USER"
mkdir -p "${LOG_FILE_PATH}"
timestamp=$(date -u +%Y.%m.%d-%H:%M)

ln "${SLURM_SUBMIT_DIR}/temp${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log" \
"${LOG_FILE_PATH}/$1~${SLURM_ARRAY_JOB_ID}~${SLURM_ARRAY_TASK_ID}~$timestamp.log"
ln "${SLURM_SUBMIT_DIR}/temp${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err" \
"${LOG_FILE_PATH}/$1~${SLURM_ARRAY_JOB_ID}~${SLURM_ARRAY_TASK_ID}~$timestamp.err"

## ============================= ##
##    VARIABLES AND FUNCTIONS    ##
## ============================= ##

blueprint_mark_name=$1
minimum_tolerated_phred_score=$2
BLUEPRINT_FULL_FILE_PATH="${RAW_DIR}/${blueprint_mark_name}"
BLUEPRINT_PROCESSED_FULL_FILE_PATH="${PROCESSED_DIR}/${blueprint_mark_name}"

## ====== FUNCTION : delete_logs() ========================
## Delete temporary log and error files then exit
## Globals: 
##   SLURM_SUBMIT_DIR
##   SLURM_JOB_ID
## Arguments:
##   exit code
## ========================================================
delete_logs(){
    rm "${SLURM_SUBMIT_DIR}/temp${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log" 
    rm "${SLURM_SUBMIT_DIR}/temp${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err"
    exit "$1"
}

if [[ -z "${minimum_tolerated_phred_score}" ]]; then
    minimum_tolerated_phred_score=20
    echo "No Phred score threshold was given, \
    using the default value of ${minimum_tolerated_phred_score}."
elif [[ "${minimum_tolerated_phred_score}" =~ ^[^0-9]+$ ]]; then
    minimum_tolerated_phred_score=20
    echo "Phred score threshold given is invalid (non-integer), \
    using the default value of ${minimum_tolerated_phred_score}."
fi

## ===================== ##
##    FILE MANAGEMENT    ##
## ===================== ##

cd "${BLUEPRINT_FULL_FILE_PATH}" || \
{ >&2 echo "ERROR: \${BLUEPRINT_FULL_FILE_PATH} - ${BLUEPRINT_FULL_FILE_PATH} \
doesn't exist, make sure you typed the epigenetic mark correctly and that you \
have ran 1_MoveFilesToSingleDirectory.sh first."; delete_logs 1; }

# Get base name of the files in 3 steps
# 1) Find all of the .bam files in the mark directory
# 2) Remove the "./" at the start of these file paths
# 3) Remove the .bam found at the end of these file names using substition
list_of_files=$(find . -type f -name "*.bam" | cut -d "/" -f 2 | sed 's/.bam//')
mkdir -p "${BLUEPRINT_PROCESSED_FULL_FILE_PATH}"


## =============================== ##
##    PARALLEL PROCESSING LOGIC    ##
## =============================== ##

# Split the directory into chunks which are determined by the array id.
total_number_of_files=$(echo "${list_of_files}" | wc -w)

# In the event that the number of files is not a multiple of the array size some
# files won't be processed if each array element processes the same number of files.
# The remaining files are processed in the highest indexed array using logic below.

number_of_files_for_each_array=$((total_number_of_files / SLURM_ARRAY_TASK_COUNT))
start_file_index=$((SLURM_ARRAY_TASK_ID * number_of_files_for_each_array))

remainder=$((total_number_of_files % SLURM_ARRAY_TASK_COUNT)) 
left_over_files=$((remainder + number_of_files_for_each_array))


if [ "${SLURM_ARRAY_TASK_ID}" -eq "${SLURM_ARRAY_TASK_COUNT}" ]; then
    files_to_process=$(find . -type f -name "*.bam" | \
    cut -d "/" -f 2 | \
    sed 's/.bam//' | \
    tail -$left_over_files)
else
    files_to_process=$(find . -type f -name "*.bam" | \
    cut -d "/" -f 2 | \
    sed 's/.bam//' | \
    head -$start_file_index | \
    tail -$number_of_files_for_each_array )
fi

## ====================== ##
##    PROCESSING STAGE    ##
## ====================== ##

echo -n "Processing the following files using Phred score threshold of "
echo "${minimum_tolerated_phred_score}:"
echo "${files_to_process}"

module purge
module load SAMtools

# The processing is characterised into 4 stages:
# 1) Create an index file, an index stats file and a stats file for the original files
# 2) Sort the .bam files, remove reads with a phred score that is below:
#     $minimum_tolerated_phred_score
#     [Note that blueprint files have already processed to 
#      remove reads with phred score below 15]
# 3) Delete intermediate files
# 4) Create an index file, an index stats file and a stats file for the processed files

for file in ${files_to_process}; do
    cd "${BLUEPRINT_FULL_FILE_PATH}" || delete_logs 1
    # 1)
    samtools index "${file}.bam" 
    samtools idxstats "${file}.bam" > "${file}.PerChromosomeStats.txt"
    samtools stats "${file}.bam" > "${file}.stats"

    # 2)
    samtools sort "${file}.bam" > \
    "${BLUEPRINT_PROCESSED_FULL_FILE_PATH}/${file}.sorted.bam"
    cd "${BLUEPRINT_PROCESSED_FULL_FILE_PATH}" || delete_logs 1

    # Need to use the -h option here to keep the headers 
    # so that the next samtools view can function properly
    samtools view -q "${minimum_tolerated_phred_score}" -h "${file}.sorted.bam" | \
    samtools sort /dev/stdin -o "${file}.sorted.filtered.bam"

    # The -h option here it to ensure the idxstats can be completed in step 4.  
    # The -F 1796 exludes reads with the following flags: 
    # a) unmapped reads b) non-primary alignment reads 
    # c) Reads that fail PCR/vendor checks d) Reads that are PCR/optical duplicates
    samtools view -F 1796 -h "${file}.sorted.filtered.bam" > \
    "${file}.sorted.filtered.noDuplicates.bam"

    # 3)
    rm "${file}.sorted.bam"
    rm "${file}.sorted.filtered.bam"

    # 4)
    samtools index "${file}.sorted.filtered.noDuplicates.bam"
    samtools idxstats "${file}.sorted.filtered.noDuplicates.bam" > \
    "${file}.sorted.filtered.noDuplicates.PerChromosomeStats.txt"
    samtools stats "${file}.sorted.filtered.noDuplicates.bam" > \
    "${file}.sorted.filtered.noDuplicates.stats"
done


## ======================= ##
##   LOG FILE MANAGEMENT   ##
## ======================= ##

echo "Job completed at:"
date -u
end_time=$(date +%s)
time_taken=$((end_time-start_time))
echo "Job took a total of: ${time_taken} seconds to complete."

delete_logs 0


