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
## -c|--config= -> Full/relative file path for configuation file directory          ||
## -m|--mark=   -> Epigenetic mark name                                             ||
## -p|--phred=  -> Phred score threshold value (default: 20)                        ||
## =================================================================================##
## OUTPUTS:                                                                         ||
## Processed .bam files                                                             ||
## Index files for raw .bam files and processed .bam files                          ||
## Per chromosome stats and general stats for raw and processed .bam files          ||
## =================================================================================##

## ===================== ##
##   ARGUMENT PARSING    ##
## ===================== ##

usage() {
cat <<EOF
=======================================================================
2_batch_ProcessBamFiles.sh
=======================================================================
Purpose: Processes .bam files by removing duplicates,
filtering out poor quality reads and sorting.
Author: Sam Fletcher
Contact: s.o.fletcher@exeter.ac.uk
Dependencies: Samtools
Inputs:
-c|--config= -> Full/relative file path for configuation file directory
-m|--mark=   -> Epigenetic mark name
-p|--phred=  -> Phred score threshold value (default: 20)
=======================================================================
EOF
    exit 0
}

needs_argument() {
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
    c | config )  needs_argument; configuration_directory="$OPTARG" ;;
    m | mark )    needs_argument; mark_name="$OPTARG" ;;
    p | phred )   needs_argument; minimum_tolerated_phred_score="$OPTARG" ;;
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
# [epigenetic mark name]~[job id]~[array id]~[date]-[time]

mv "${SLURM_SUBMIT_DIR}/temp${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log" \
"${LOG_FILE_PATH}/$2~${SLURM_ARRAY_JOB_ID}~${SLURM_ARRAY_TASK_ID}~${timestamp:=}.log"
mv "${SLURM_SUBMIT_DIR}/temp${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err" \
"${LOG_FILE_PATH}/$2~${SLURM_ARRAY_JOB_ID}~${SLURM_ARRAY_TASK_ID}~$timestamp.err"

## =============== ##
##    VARIABLES    ##
## =============== ##

RAW_FULL_FILE_PATH="${RAW_DIR}/${mark_name}"
PROCESSED_FULL_FILE_PATH="${PROCESSED_DIR}/${mark_name}"

## ====== DEFAULTS =============================================================
if ! [[ "${minimum_tolerated_phred_score}" =~ ^[0-9]+$ ]]; then
    minimum_tolerated_phred_score=20
    echo "Phred score threshold given is invalid, "\
    "using the default value of ${minimum_tolerated_phred_score}."
fi
## =============================================================================

## ===================== ##
##    FILE MANAGEMENT    ##
## ===================== ##

cd "${RAW_FULL_FILE_PATH}" || \
{ >&2 echo "ERROR: [\${RAW_FULL_FILE_PATH} - ${RAW_FULL_FILE_PATH}] \
doesn't exist, make sure you typed the epigenetic mark correctly and that you \
have ran 1_MoveFilesToSingleDirectory.sh first."; finishing_statement 1; }

# This is in place in case a previous run of this script timed out or was
# cancelled
rm ./*samtools*tmp*.bam

list_of_files=$(find . -type f -name "*.bam")
mkdir -p "${PROCESSED_FULL_FILE_PATH}"

cd "${PROCESSED_FULL_FILE_PATH}" || finishing_statement 1

# We want to remove any files from previous runs in case this script previously
# timed out or was cancelled (helps with later scripts) 
rm ./*

## =============================== ##
##    PARALLEL PROCESSING LOGIC    ##
## =============================== ##

# Split the directory into chunks which are determined by the array id.
total_number_of_files=$(echo "${list_of_files}" | wc -w)

# In the event that the number of files is not a multiple of the array size some
# files won't be processed if each array element processes the same number of files.
# The remaining files are processed in the highest indexed array using logic below.

number_files_for_each_array=$((total_number_of_files / SLURM_ARRAY_TASK_MAX))
start_file_index=$((SLURM_ARRAY_TASK_ID * number_files_for_each_array))

remainder=$((total_number_of_files % SLURM_ARRAY_TASK_COUNT)) 
left_over_files=$((remainder + number_files_for_each_array))


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
    tail -$number_files_for_each_array )
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
    echo "Processing ${file}..."

    base_name=$(basename "${file}" .bam)

    # Processing includes multiple steps, we use the 'success' variable to
    # keep track of any failures in the pipeline (overall errors for a file).
    # to make error files easier to read.
    success=0
    cd "${RAW_FULL_FILE_PATH}" || finishing_statement 1
    # 1)
    samtools index "${base_name}.bam" 
    samtools idxstats "${base_name}.bam" > "${base_name}.PerChromosomeStats.txt"
    samtools stats "${base_name}.bam" > "${base_name}.stats"

    # 2)
    samtools sort "${base_name}.bam" > \
    "${PROCESSED_FULL_FILE_PATH}/${base_name}.sorted.bam"
    if [[ $? == 1 ]]; then 
        { >&2 echo "Sorting failed for ${file}."; }
        success=1
    fi
    
    cd "${PROCESSED_FULL_FILE_PATH}" || finishing_statement 1

    # Need to use the -h option here to keep the headers 
    # so that the next samtools view can function properly
    samtools view -q "${minimum_tolerated_phred_score}" -h "${base_name}.sorted.bam" | \
    samtools sort /dev/stdin -o "${base_name}.sorted.filtered.bam"
    if (( PIPESTATUS[0] != 0 || PIPESTATUS[1] != 0 )); then 
        { >&2 echo "Filtering failed for ${file}."; }
        success=1
    fi

    # The -h option here it to ensure the idxstats can be completed in step 4.  
    # The -F 1796 exludes reads with the following flags: 
    # a) unmapped reads b) non-primary alignment reads 
    # c) Reads that fail PCR/vendor checks d) Reads that are PCR/optical duplicates
    samtools view -F 1796 -h "${base_name}.sorted.filtered.bam" > \
    "${base_name}.sorted.filtered.noDuplicates.bam"
    if [[ $? == 1 ]]; then 
        { >&2 echo "Remove duplicates failed for ${file}."; }
        success=1
    fi

    # 3)
    rm "${base_name}.sorted.bam"
    rm "${base_name}.sorted.filtered.bam"

    # 4)
    samtools index "${base_name}.sorted.filtered.noDuplicates.bam"
    samtools idxstats "${base_name}.sorted.filtered.noDuplicates.bam" > \
    "${base_name}.sorted.filtered.noDuplicates.PerChromosomeStats.txt"
    samtools stats "${base_name}.sorted.filtered.noDuplicates.bam" > \
    "${base_name}.sorted.filtered.noDuplicates.stats"

    if [[ success -eq 0 ]]; then
        echo "${file} processed successfully."
    else
        { echo -e "${file} had an error during processing. Check error log:\n\
        ${LOG_FILE_PATH}/$2~${SLURM_ARRAY_JOB_ID\
        }~${SLURM_ARRAY_TASK_ID}~$timestamp.err"; }
    fi
done

finishing_statement 0


