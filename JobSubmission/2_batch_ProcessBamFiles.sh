#!/bin/bash
# Export all environment variables to the batch job
#SBATCH --export=ALL
# Submit to the mrc queue for faster queue times
#SBATCH -p mrcq
# Tests have thus far shown a linear relationship between file size and time
# Current approximation is: [time (mins)] = 0.1 + 9*[file size (GB)]
#SBATCH --time=12:00:00 
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
#SBATCH --job-name=2_Processing

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
## AUTHOR: Sam Fletcher                                                             ||
## CONTACT: s.o.fletcher@exeter.ac.uk                                               ||
## CREATED: November 2023                                                           ||
## =================================================================================##
## PREREQUISITES: All .bam files for a specific epigenetic mark must be in 1 folder ||
## =================================================================================##
## DEPENDENCIES: Samtools                                                           ||
## =================================================================================##
## INPUTS:                                                                          ||
## $1 -> Full (or relative) file path for configuation file directory               ||
## $2 -> Epigenetic mark to process                                                 ||
## $3 -> Phred score threshold value (default: 20)                                  ||
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
    echo "==================================================================="
    echo "Purpose: Processes .bam files by removing duplicates,"
    echo "filtering out poor quality reads and sorting."
    echo "Author: Sam Fletcher"
    echo "Contact: s.o.fletcher@exeter.ac.uk"
    echo "Dependencies: Samtools"
    echo "Inputs:"
    echo "\$1 -> Full (or relative) file path for configuation file directory"
    echo "\$2 -> Name of epigenetic mark"
    echo "\$3 -> Phred score threshold value (default: 20)"
    echo "Optional:"
    echo "Specify --array in sbatch options, to set a custom array size."
    echo "==================================================================="
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
# [epigenetic mark name]~[job id]~[array id]~[date]-[time]

mv "${SLURM_SUBMIT_DIR}/temp${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log" \
"${LOG_FILE_PATH}/$2~${SLURM_ARRAY_JOB_ID}~${SLURM_ARRAY_TASK_ID}~${timestamp:=}.log"
mv "${SLURM_SUBMIT_DIR}/temp${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err" \
"${LOG_FILE_PATH}/$2~${SLURM_ARRAY_JOB_ID}~${SLURM_ARRAY_TASK_ID}~$timestamp.err"

## =============== ##
##    VARIABLES    ##
## =============== ##

mark_name=$2
minimum_tolerated_phred_score=$3
RAW_FULL_FILE_PATH="${RAW_DIR}/${mark_name}"
PROCESSED_FULL_FILE_PATH="${PROCESSED_DIR}/${mark_name}"

## ====== DEFAULTS ====================================================================
if ! [[ "${minimum_tolerated_phred_score}" =~ ^[0-9]+$ ]]; then
    minimum_tolerated_phred_score=20
    echo "Phred score threshold given is invalid, "\
    "using the default value of ${minimum_tolerated_phred_score}."
fi
## ====================================================================================

## ===================== ##
##    FILE MANAGEMENT    ##
## ===================== ##

cd "${RAW_FULL_FILE_PATH}" || \
{ >&2 echo "ERROR: [\${RAW_FULL_FILE_PATH} - ${RAW_FULL_FILE_PATH}] \
doesn't exist, make sure you typed the epigenetic mark correctly and that you \
have ran 1_MoveFilesToSingleDirectory.sh first."; finishing_statement 1; }

# Get base name of the files in 3 steps
# 1) Find all of the .bam files in the mark directory
# 2) Remove the "./" at the start of these file paths
# 3) Remove the .bam found at the end of these file names using substition
list_of_files=$(find . -type f -name "*.bam" | cut -d "/" -f 2 | sed 's/.bam//')
mkdir -p "${PROCESSED_FULL_FILE_PATH}"


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


if [[ "${SLURM_ARRAY_TASK_ID}" -eq "${SLURM_ARRAY_TASK_COUNT}" ]]; then
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

module purge
module load SAMtools

# The processing is split into 4 steps:
# 1) Create an index file, an index stats file and a stats file for the original files
# 2) Sort the .bam files and remove reads with a phred score that is below
#    $minimum_tolerated_phred_score [Note that blueprint files have already processed 
#    to remove reads with phred score below 15]   
# 3) Delete intermediate files
# 4) Create an index file, an index stats file and a stats file for the processed files

for file in ${files_to_process}; do
    echo "Processing ${file}.bam..."
    success=0
    cd "${RAW_FULL_FILE_PATH}" || finishing_statement 1
    # 1)
    samtools index "${file}.bam" 
    samtools idxstats "${file}.bam" > "${file}.PerChromosomeStats.txt"
    samtools stats "${file}.bam" > "${file}.stats"

    # 2)
    samtools sort "${file}.bam" > \
    "${PROCESSED_FULL_FILE_PATH}/${file}.sorted.bam"
    if [[ $? == 1 ]]; then 
        { >&2 echo "Sorting failed for ${file}."; }
        success=1
    fi
    
    cd "${PROCESSED_FULL_FILE_PATH}" || finishing_statement 1

    # Need to use the -h option here to keep the headers 
    # so that the next samtools view can function properly
    samtools view -q "${minimum_tolerated_phred_score}" -h "${file}.sorted.bam" | \
    samtools sort /dev/stdin -o "${file}.sorted.filtered.bam"
    if (( PIPESTATUS[0] != 0 || PIPESTATUS[1] != 0 )); then 
        { >&2 echo "Filtering failed for ${file}."; }
        success=1
    fi

    # The -h option here it to ensure the idxstats can be completed in step 4.  
    # The -F 1796 exludes reads with the following flags: 
    # a) unmapped reads b) non-primary alignment reads 
    # c) Reads that fail PCR/vendor checks d) Reads that are PCR/optical duplicates
    samtools view -F 1796 -h "${file}.sorted.filtered.bam" > \
    "${file}.sorted.filtered.noDuplicates.bam"
    if [[ $? == 1 ]]; then 
        { >&2 echo "Remove duplicates failed for ${file}."; }
        success=1
    fi

    # 3)
    rm "${file}.sorted.bam"
    rm "${file}.sorted.filtered.bam"

    # 4)
    samtools index "${file}.sorted.filtered.noDuplicates.bam"
    samtools idxstats "${file}.sorted.filtered.noDuplicates.bam" > \
    "${file}.sorted.filtered.noDuplicates.PerChromosomeStats.txt"
    samtools stats "${file}.sorted.filtered.noDuplicates.bam" > \
    "${file}.sorted.filtered.noDuplicates.stats"

    if [[ success -eq 0 ]]; then
        echo "${file}.bam processed successfully."
    else
        { echo -e "${file}.bam had an error during processing. Check error log:\n\
        ${LOG_FILE_PATH}/$1~${SLURM_ARRAY_JOB_ID\
        }~${SLURM_ARRAY_TASK_ID}~$timestamp.err"; }
    fi
done

finishing_statement 0


