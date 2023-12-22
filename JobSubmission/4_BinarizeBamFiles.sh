#!/bin/bash
# Export all environment variables to the batch job
#SBATCH --export=ALL
# Submit to the mrc queue for faster queue times
#SBATCH -p mrcq 
# Consult information/Processing_Times.md for expected time
#SBATCH --time=01:00:00 
#SBATCH -A Research_Project-MRC190311 
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=16
# Memory consumption is generally fairly low
#SBATCH --mem=4G 
# Send an email after the job is done
#SBATCH --mail-type=END 
# Temporary log file, later to be removed
#SBATCH --output=temp%j.log
# Temporary error file, later to be removed
#SBATCH --error=temp%j.err
#SBATCH --job-name=Binarization

## =================================================================================##
##                                                                                  ||
##                                     PREAMBLE                                     ||
##                                                                                  ||
## =================================================================================##
## PURPOSE:                                                                         ||
## The Subsampled files for each mark will now be binarized through the use of      ||
## ChromHMM's BinarizeBam command. This script is to be ran after all of the        ||
## epigenetic marks that one wants to inspect have been subsampled. The .bam files  ||
## need to be binarized so that they can be used by ChromHMM's LearnModel command.  ||
## =================================================================================##
## AUTHOR: Sam Fletcher                                                             ||
## CONTACT: s.o.fletcher@exeter.ac.uk                                               ||
## CREATED: November 2023                                                           ||
## =================================================================================##
## PREREQUISITES: Run: 3_SubsampleBamFiles.sh                                       ||
## =================================================================================##
## DEPENDENCIES: Java, ChromHMM                                                     ||
## =================================================================================##
## INPUTS:                                                                          ||
## $1 -> Bin size to be used by BinarizeBam command (default: 200)                  ||
## $2 -> Sample size used in 3_SubsampleBamFiles.sh                                 ||
## =================================================================================##
## OUTPUTS:                                                                         ||
## Binary signal files for every chromosome in the dataset except for mitochondrial ||
## DNA                                                                              ||
## =================================================================================##

## ======================== ##
##    HELP FUNCTIONALITY    ##
## ======================== ##

if [[ "$1" == "--help" || "$1" == "-h" ]]; then
    echo "================================================================"
    echo "Purpose: Creates a 'cell mark file table' and uses"
    echo "ChromHMM's BinarizeBam command to binarize bam files."
    echo "Author: Sam Fletcher"
    echo "Contact: s.o.fletcher@exeter.ac.uk"
    echo "Dependencies: Java, ChromHMM"
    echo "Inputs:"
    echo "\$1 -> Bin size to be used by BinarizeBam command (default: 200)"
    echo "\$2 -> Sample size used in 3_SubsampleBamFiles.sh"
    echo "================================================================"
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
/scripts/integrative/ChromHMM_OptimumStates/config/config.txt"

# Rename the output and error files to have format:
# BinSize-[bin size]~[job id]~[date]-[time]
# This requires a hard link as you cannot rename log files
# whilst running the script without a wrapper function
LOG_FILE_PATH="${LOG_DIR}/$SLURM_JOB_NAME/$USER"
mkdir -p "${LOG_FILE_PATH}"
timestamp=$(date -u +%Y.%m.%d-%H:%M)

ln "${SLURM_SUBMIT_DIR}/temp${SLURM_JOB_ID}.log" \
"${LOG_FILE_PATH}/BinSize-$1~${SLURM_JOB_ID}~$timestamp.log"
ln "${SLURM_SUBMIT_DIR}/temp${SLURM_JOB_ID}.err" \
"${LOG_FILE_PATH}/BinSize-$1~${SLURM_JOB_ID}~$timestamp.err"

## ============================= ##
##    VARIABLES AND FUNCTIONS    ##
## ============================= ##

bin_size=$1
sample_size=$2

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
if ! [[ "${bin_size}" =~ ^[0-9]+$ ]]; then
    bin_size=200
    echo "Invalid bin size given, using the default value of ${bin_size} instead."
fi

# 'Intelligently' find the sample size using first file name in subsampled directory
if [[ -z "${sample_size}" ]]; then
    cd "${SUBSAMPLED_DIR}" || \
    { >&2 echo "ERROR: \${SUBSAMPLED_DIR} - ${SUBSAMPLED_DIR} doesn't exist, make \
    sure config.txt is pointing to the correct directory"; finishing_statement 1; }

    sample_size=$(find . -type f -name "Subsampled*" | head -1 | cut -d "." -f 3)
    echo -e "WARNING: No sample size was given.\n\
    Assuming that ${sample_size} is the desired sample size..."
fi
if [[ -z "${sample_size}" ]]; then
    { >&2 echo -e "ERROR: No sample size even after fail safe. Please run\n\
    3_SubsampleBamFiles.sh before running this script" ;}
fi
## ====================================================================================

## ================== ##
##   FILE EXISTENCE   ##
## ================== ##

cd "${SUBSAMPLED_DIR}" || \
{ >&2 echo "ERROR: \${SUBSAMPLED_DIR} - ${SUBSAMPLED_DIR} doesn't exist, make \
sure config.txt is pointing to the correct directory"; finishing_statement 1; }

if [[ -z "$(ls -A)" ]]; then
    { >&2 echo "ERROR: \${SUBSAMPLED_DIR} - ${SUBSAMPLED_DIR} is empty.
    Ensure that 3_SubsampleBamFiles.sh has been ran before this script."; }

    finishing_statement 1
fi

## =================================== ##
##    CREATE A CELL MARK FILE TABLE    ##
## =================================== ##

# As we have already merged files in 3_SubsampleBamFiles.sh, 
# the cell mark file table will just be 1 file for each mark that has been processed
# The only added level of complexity is obtaining the mark name from the file names.

rm -f "cellmarkfiletable.txt" 
for file in *.bam; do
    echo -ne "Mature_Neutorphil_SampleSize_${sample_size}_BinSize_${bin_size}\t" >>\
    "cellmarkfiletable.txt"
    # The subsampled files are named: subsampled.[SampleSize].[mark_name].bam. 
    # Below extracts the mark name
    mark_name=$(echo "$file" | cut -d "." -f 3) 
    echo -ne "${mark_name}\t" >> "cellmarkfiletable.txt"
    echo "$file" >> "cellmarkfiletable.txt"
done

## ================================= ##
##    BINARIZATION USING CHROMHMM    ##
## ================================= ##

echo "Binarizing subsampled bam files found in \${SUBSAMPLED_DIR} - \
${SUBSAMPLED_DIR} with sample size: ${sample_size} using a bin size \
of: ${bin_size}."

cd "${BINARY_DIR}" || \
{ echo "ERROR: \${BINARY_DIR} - ${BINARY_DIR} doesn't exist, \
make sure config.txt is pointing to the correct directory"; finishing_statement 1; }
rm ./*.txt*

module purge
module load Java

# Binarize the files in the subsampled directory.
# The blueprint data uses GChr37 assembly (which is equivalent to UCSC's hg19).
java -mx4G -jar "${CHROMHMM_MAIN_DIR}/ChromHMM.jar" BinarizeBam -b "${bin_size}" \
-gzip "${CHROMHMM_CHROM_SIZES}/hg19.txt" "${SUBSAMPLED_DIR}" \
"${SUBSAMPLED_DIR}/cellmarkfiletable.txt" "${BINARY_DIR}"

# Optional: One may not want to keep the mitochondrial DNA in the analysis
# rm "${BINARY_DIR}\
# /Mature_Neutorphil_SampleSize_${sample_size}_BinSize_${bin_size}_chrM_binary.txt.gz"

finishing_statement 0