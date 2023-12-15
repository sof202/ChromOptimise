#!/bin/bash
# export all environment variables to the batch job
#SBATCH --export=ALL
# submit to the mrc queue for faster queue times
#SBATCH -p mrcq
# Time for a single bam file of size: 1134MB was: ~10 minutes (scale up accordingly)
#SBATCH --time=01:00:00 
#SBATCH -A Research_Project-MRC190311 
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=16
# Make sure it is not higher than the number of files being processed
# as this ends up with all files being processed by max index array
#SBATCH --array=1-4
# specify bytes memory to reserve
#SBATCH --mem=10G 
# Send an email after the job is done
#SBATCH --mail-type=END
# Temporary log file, later to be removed
#SBATCH --output=temp%a.log
# Temporary error file, later to be removed
#SBATCH --error=temp%a.err
#SBATCH --job-name=Processing

## -------------------------------------------------------------------------------------------- ##
##                                                                                              ##
##                                            PREAMBLE                                          ##
##                                                                                              ##
## -------------------------------------------------------------------------------------------- ##
##                                            PURPOSE                                           ##
##             1_RawBamFiles contains .bam files for each epigentic mark, however these         ##
##     have yet to be processed, they still contain duplicates, multimapped fragments and low   ##
##                        quality reads. This script processes the bam files.                   ##
## -------------------------------------------------------------------------------------------- ##
##                        AUTHOR: Sam Fletcher s.o.fletcher@exeter.ac.uk                        ##
##                                     CREATED: November 2023                                   ##
## -------------------------------------------------------------------------------------------- ##
##                                         PREREQUISITES                                        ##
##                All .bam files for a specific epigenetic mark must be in 1 folder             ##
## -------------------------------------------------------------------------------------------- ##
##                                          DEPENDENCIES                                        ##
##                                            SAMtools                                          ##
## -------------------------------------------------------------------------------------------- ##
##                                             INPUTS                                           ##
##                               $1 -> Epigenetic mark to process                               ##
##                               $2 -> Phred score threshold value                              ##
## -------------------------------------------------------------------------------------------- ##
##                                            OUTPUTS                                           ##
##                                      Processed .bam files                                    ##
##                      Index files for raw .bam files and processed .bam files                 ##
##             Per chromosome stats and general stats for raw and processed .bam files          ##
## -------------------------------------------------------------------------------------------- ##

## ------------------------ ##
##    HELP FUNCTIONALITY    ##
## ------------------------ ##

if [[ "$1" == "--help" || "$1" == "-h" ]]; then
    echo "====================================================="
    echo "Purpose: Processes .bam files by removing duplicates,"
    echo "filtering out poor quality reads and sorting."
    echo "Author: Sam Fletcher"
    echo "Contact: s.o.fletcher@exeter.ac.uk"
    echo "Dependencies: Samtools"
    echo "Inputs:"
    echo "\$1 -> Name of epigenetic mark"
    echo "\$2 -> Phred score threshold value"
    echo "====================================================="
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
# [epigenetic mark name]~[job id]~[array id]~[date]-[time]
# This requires a hard link as you cannot rename log files
# whilst running the script without a wrapper function
LOG_FILE_PATH="${LOG_DIR}/$USER/$SLURM_JOB_NAME/"
timestamp=$(date -u +%Y.%m.%d-%H:%M)
ln "${SLURM_SUBMIT_DIR}/temp${SLURM_ARRAY_TASK_ID}.log" \
"${LOG_FILE_PATH}/$1~${SLURM_ARRAY_JOB_ID}~${SLURM_ARRAY_TASK_ID}~$timestamp.log"
ln "${SLURM_SUBMIT_DIR}/temp${SLURM_ARRAY_TASK_ID}.err" \
"${LOG_FILE_PATH}/$1~${SLURM_ARRAY_JOB_ID}~${SLURM_ARRAY_TASK_ID}~$timestamp.err"

# Open the main blueprint directory
BLUEPRINT_MARK_NAME=$1
BLUEPRINT_FULL_FILE_PATH="${RAW_DIR}/${BLUEPRINT_MARK_NAME}"
BLUEPRINT_PROCESSED_FULL_FILE_PATH="${PROCESSED_DIR}/${BLUEPRINT_MARK_NAME}"
MINIMUM_TOLERATED_PHRED_SCORE=$2

# Check if the directory with the epigenetic mark actually exists
if [ -d "${BLUEPRINT_FULL_FILE_PATH}" ]; then
    echo "Changing directory to: ${BLUEPRINT_FULL_FILE_PATH}" 
    cd "${BLUEPRINT_FULL_FILE_PATH}" || exit 1
else
    echo "Directory does not exist yet. Make sure you typed the epigenetic mark correctly and that you have ran 1_MoveFilesToSingleDirectory.sh first"
    echo "Aborting..."

    # Remove temporary log files
    rm "${SLURM_SUBMIT_DIR}/temp${SLURM_ARRAY_TASK_ID}.log"
    rm "${SLURM_SUBMIT_DIR}/temp${SLURM_ARRAY_TASK_ID}.err"
    exit 1
fi

# Set a default value for minimum phred score in case one is not given
if [ -z "${MINIMUM_TOLERATED_PHRED_SCORE}" ]; then
    MINIMUM_TOLERATED_PHRED_SCORE=20
    echo "No Phred score threshold was given, using default value of 20."
fi
echo "Processing .bam files using Phred score threshold of: ${MINIMUM_TOLERATED_PHRED_SCORE} for epigenetic mark: ${BLUEPRINT_MARK_NAME}"


## ----------------------------- ##
##    PROCESSING OF BAM FILES    ##
## ----------------------------- ##

# We purge all of the modules to ensure there are no conflicts and load the SAMtools module for the processing steps
module purge
module load SAMtools

# Get base name of the files in 3 steps
# First we find all of the .bam files in the mark directory
# We then remove the "./" at the start of these file paths
# Then we remove the .bam found at the end of these file names using a substition (search for .bam and replace with the empty string)
List_Of_Files=$(find . -type f -name "*.bam" | cut -d "/" -f 2 | sed 's/.bam//')
mkdir -p "${BLUEPRINT_PROCESSED_FULL_FILE_PATH}"


## ------------------------------- ##
##    PARALLEL PROCESSING LOGIC    ##
## ------------------------------- ##


# Effectively what we are doing here is splitting the directory into chunks which are determined by the array id.
Total_Number_Of_Files=$(echo "${List_Of_Files}" | wc -w)

Number_Of_Files_For_Each_Array=$((Total_Number_Of_Files / SLURM_ARRAY_TASK_COUNT))
Remainder=$((Total_Number_Of_Files % SLURM_ARRAY_TASK_COUNT)) # This is for those final files that won't be picked up otherwise if the number of files is not a multiple of the array size
Start_File_Index=$((SLURM_ARRAY_TASK_ID * Number_Of_Files_For_Each_Array))
Left_Over_Files=$((Remainder + Number_Of_Files_For_Each_Array))

# In the event that the number of files is not a multiple of the array size, we need to pick up the last few files in the highest Id array using the remainder calculated above
# Head and tail logic is used to extract files from the list. The same could be done with 'cut' if you wish. I think this is more transparent.
# The 'Total_Number_Of_Files' variable is rather redundant. But again, I think this way of displaying the code is a little more transparent. 
if [ "${SLURM_ARRAY_TASK_ID}" -eq "${SLURM_ARRAY_TASK_COUNT}" ]; then
    Files_To_Process=$(find . -type f -name "*.bam" | cut -d "/" -f 2 | sed 's/.bam//' | tail -$Left_Over_Files)
else
    Files_To_Process=$(find . -type f -name "*.bam" | cut -d "/" -f 2 | sed 's/.bam//' | head -$Start_File_Index | tail -$Number_Of_Files_For_Each_Array )
fi

echo "Processing the following files:"
echo "${Files_To_Process}"


## ---------------------- ##
##    PROCESSING STAGE    ##
## ---------------------- ##

# The processing is characterised into 4 stages
# 1) Create an index file, an index stats file and a stats file for the original files (for comparisons)
# 2) sort the .bam files, remove reads with a phred score that is below the threshold given by user (default is 20) 
# [Note that blueprint files have already processed to remove reads with phred score below 15]
# 3) Delete intermediate files
# 4) Create an index file, an index stats file and a stats file for the processed files (for comparisons)
for file in ${Files_To_Process}; do
    cd "${BLUEPRINT_FULL_FILE_PATH}" || exit 1
    # 1)
    samtools index "${file}.bam" 
    samtools idxstats "${file}.bam" > "${file}.PerChromosomeStats.txt"
    samtools stats "${file}.bam" > "${file}.stats"

    # 2)
    samtools sort "${file}.bam" > "${BLUEPRINT_PROCESSED_FULL_FILE_PATH}/${file}.sorted.bam"
    cd "${BLUEPRINT_PROCESSED_FULL_FILE_PATH}" || exit 1 
    samtools view -q "${MINIMUM_TOLERATED_PHRED_SCORE}" -h "${file}.sorted.bam" | samtools sort /dev/stdin -o "${file}.sorted.filtered.bam"
    # Need to use the -h option here to keep the headers so that the next samtools view can function properly
    samtools view -F 1796 -h "${file}.sorted.filtered.bam" > "${file}.sorted.filtered.noDuplicates.bam"
    # The -h option here it to ensure the idxstats can be completed in step 4.  
    # The -F 1796 exludes the following flags: 
    # 1) unmapped reads 2) non-primary alignment reads 3) Reads that fail PCR/vendor checks 4) Reads that are PCR/optical duplicates

    # 3)
    rm "${file}.sorted.bam"
    rm "${file}.sorted.filtered.bam"

    # 4)
    samtools index "${file}.sorted.filtered.noDuplicates.bam"
    samtools idxstats "${file}.sorted.filtered.noDuplicates.bam" > "${file}.sorted.filtered.noDuplicates.PerChromosomeStats.txt"
    samtools stats "${file}.sorted.filtered.noDuplicates.bam" > "${file}.sorted.filtered.noDuplicates.stats"
done


## ----------------------- ##
##   LOG FILE MANAGEMENT   ##
## ----------------------- ##

# Finishing message
echo "Job completed at:"
date -u
end_time=$(date +%s)
time_taken=$((end_time-start_time))
echo "Job took a total of: ${time_taken} seconds to complete"

# Removing temporary log files
rm "${SLURM_SUBMIT_DIR}/temp${SLURM_ARRAY_TASK_ID}.log"
rm "${SLURM_SUBMIT_DIR}/temp${SLURM_ARRAY_TASK_ID}.err"


