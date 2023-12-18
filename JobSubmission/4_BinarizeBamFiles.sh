#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -p mrcq # submit to the mrc queue for faster queue times
#SBATCH --time=01:00:00 # In testing, binarizing 6 files (total size being 8GB), a total of 22 minutes ellapsed for this binarization process to complete
#SBATCH -A Research_Project-MRC190311 
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=16 
#SBATCH --mem=100G # specify bytes memory to reserve, lots is required for this script as the binarization process stores a lot of data in memory
#SBATCH --mail-type=END # Send an email after the job is done
# Temporary log file, later to be removed
#SBATCH --output=temp%j.log
# Temporary error file, later to be removed
#SBATCH --error=temp%j.err
#SBATCH --job-name=Binarization

## -------------------------------------------------------------------------------------------- ##
##                                                                                              ##
##                                            PREAMBLE                                          ##
##                                                                                              ##
## -------------------------------------------------------------------------------------------- ##
##                                            PURPOSE                                           ##
##    The Subsampled files for each mark will now be binarized through the use of ChromHMM's    ##
##   BinarizeBam command. This script is to be ran after all of the epigenetic marks that one   ##
## wants to inspect have been subsampled. The .bam files need to be binarized so that they can  ##
##    be used by ChromHMM's LearnModel command. This is a crucial step in the pipeline as the   ##
##  bin size is decided here, which will likely have a major effect on the quality of the data. ##
## -------------------------------------------------------------------------------------------- ##
##                        AUTHOR: Sam Fletcher s.o.fletcher@exeter.ac.uk                        ##
##                                     CREATED: November 2023                                   ##
## -------------------------------------------------------------------------------------------- ##
##                                         PREREQUISITES                                        ##
##                                  Run: 3_SubsampleBamFiles.sh                                 ##
## -------------------------------------------------------------------------------------------- ##
##                                          DEPENDENCIES                                        ##
##                                              Java                                            ##
##                                            ChromHMM                                          ##
## -------------------------------------------------------------------------------------------- ##
##                                             INPUTS                                           ##
##                        $1 -> Bin size to be used by BinarizeBam command                      ##
##                        $2 -> Sample size used in 3_SubsampleBamFiles.sh                      ##
## -------------------------------------------------------------------------------------------- ##
##                                            OUTPUTS                                           ##
##       Binary signal files for every chromosome in the dataset except for mitochondrial DNA   ##
## -------------------------------------------------------------------------------------------- ##

## ------------------------ ##
##    HELP FUNCTIONALITY    ##
## ------------------------ ##

if [[ "$1" == "--help" || "$1" == "-h" ]]; then
    echo "======================================================================================================="
    echo "Purpose: Creates a 'cell mark file table' and uses ChromHMM's BinarizeBam command to binarize bam files"
    echo "Author: Sam Fletcher"
    echo "Contact: s.o.fletcher@exeter.ac.uk"
    echo "Dependencies: Java, ChromHMM"
    echo "Inputs:"
    echo "\$1 -> Bin size to be used by BinarizeBam command"
    echo "\$2 -> Sample size used in 3_SubsampleBamFiles.sh"
    echo "======================================================================================================="
    exit 3
fi

## ------------ ##
##    SET UP    ##
## ------------ ##

echo "Job '$SLURM_JOB_NAME' started at:"
date -u

start_time=$(date +%s)

# Activate config.txt to access all file paths
# CHANGE THIS TO YOUR OWN CONFIG FILE
echo "Loading config file..."
source "/lustre/projects/Research_Project-MRC190311\
/scripts/integrative/blueprint/config/config.txt"

# Rename the output and error files to have format:
# BinSize-[bin size]~[job id]~[date]-[time]
# This requires a hard link as you cannot rename log files
# whilst running the script without a wrapper function
LOG_FILE_PATH="${LOG_DIR}/$USER/$SLURM_JOB_NAME/"
mkdir -p "${LOG_FILE_PATH}"
timestamp=$(date -u +%Y.%m.%d-%H:%M)

ln "${SLURM_SUBMIT_DIR}/temp${SLURM_JOB_ID}.log" \
"${LOG_FILE_PATH}/BinSize-$1~${SLURM_JOB_ID}~$timestamp.log"
ln "${SLURM_SUBMIT_DIR}/temp${SLURM_JOB_ID}.err" \
"${LOG_FILE_PATH}/BinSize-$1~${SLURM_JOB_ID}~$timestamp.err"



# Initialise variables
BIN_SIZE=$1
SAMPLE_SIZE=$2

# Check if the 3_SubsampledBamFiles directory is empty. If it is, then abort to avoid errors and display helpful information to the user
cd "${SUBSAMPLED_DIR}" || { echo "Subsample directory doesn't exist, make sure config.txt is pointing to the correct directory"; exit 1; }
if [ -z "$(ls -A)" ]; then
    echo "3_SubsampledBamFiles is empty, ensure that 3_SubsampleBamFiles.sh has been ran before this script."
    echo "Aborting..."

    # Remove temporary log files
    rm "${SLURM_SUBMIT_DIR}/temp${SLURM_JOB_ID}.log"
    rm "${SLURM_SUBMIT_DIR}/temp${SLURM_JOB_ID}.err"

    exit 1
fi
## SET DEFAULTS FOR VARIABLES ##
## ============================================================================================================================================ ##
if [ -z "$BIN_SIZE" ]; then
    echo "No bin size was given, using the default value of 200 instead."
    BIN_SIZE=200
fi

# Obtains a sample size default intelligently by searching through the subsample directory
if [ -z "${SAMPLE_SIZE}" ]; then
    echo "No sample size was given. Assuming that the first file in the subsampled directory uses the correct sample size..."
    cd "${SUBSAMPLED_DIR}" || { echo "Subsample directory doesn't exist, make sure config.txt is pointing to the correct directory"; exit 1; }
    SAMPLE_SIZE=$(find . -type f -name "Sub*" | head -1 | cut -d "." -f 3)
fi
## ============================================================================================================================================ ##


echo "Binarizing subsampled bam files found in 3_SubsampledBamFiles with sample size: ${SAMPLE_SIZE} using a bin size of: ${BIN_SIZE}."

## ------------------------------------------------ ##
##    CREATE A CELL MARK FILE TABLE FOR CHROMHMM    ##
## ------------------------------------------------ ##

# As we have already merged files in 3_SubsampleBamFiles.sh, the cell mark file table will just be 1 file for each mark that has been processed
# The only added level of complexity is obtaining the mark name from the file names. The cell type for all blueprint data (currently) is the same, so no need to do anything there

rm "cellmarkfiletable.txt" # Reset the cell mark file table in case it was already made in a previous run through of the pipeline
for file in *.bam; do
    echo -ne "Mature_Neutorphil_SampleSize_${SAMPLE_SIZE}_BinSize_${BIN_SIZE}\t" >> "cellmarkfiletable.txt" #echo -ne is used here to ensure that no newline is inserted and instead a tab is inserted after the cell name (as this file is supposed to be tsv).
    Mark_Name=$(echo "$file" | cut -d "." -f 3) # In 3_SubsampleBamFiles.sh, the output files were named: subsampled.(Mark_Name).bam. Thus we need to extract the 2nd entry in this dot separated file name
    echo -ne "${Mark_Name}\t" >> "cellmarkfiletable.txt"
    echo "$file" >> "cellmarkfiletable.txt"
done




## --------------------------------- ##
##    BINARIZATION USING CHROMHMM    ##
## --------------------------------- ##

# Delete any binary files that have been created in previous run throughs
cd "${BINARY_DIR}" || { echo "Binary directory doesn't exist, make sure config.txt is pointing to the correct directory"; exit 1; }
rm ./*.txt*

module purge
module load Java

# Binarize the files in the subsampled directory using the user specified (or default) bin size. The blueprint data is using GChr37 (which is equivalent to UCSC's hg19).
java -mx30G -jar "${CHROMHMM_MAIN_DIR}/ChromHMM.jar" BinarizeBam -b "${BIN_SIZE}" -gzip "${CHROMHMM_CHROM_SIZES}/hg19.txt" "${SUBSAMPLED_DIR}" "${SUBSAMPLED_DIR}/cellmarkfiletable.txt" "${BINARY_DIR}"

# Delete Mitochondrial DNA binary file
rm "${BINARY_DIR}/Mature_Neutorphil_SampleSize_${SAMPLE_SIZE}_BinSize_${BIN_SIZE}_chrM_binary.txt.gz"



## ----------------------- ##
##   LOG FILE MANAGEMENT   ##
## ----------------------- ##


# Finishing message
echo "Job completed at:"
date -u
end_time=$(date +%s)
time_taken=$(("$end_time"-"$start_time"))
echo "Job took a total of: ${time_taken} seconds to complete"

# Removing temporary log files
rm "${SLURM_SUBMIT_DIR}/temp${SLURM_JOB_ID}.log"
rm "${SLURM_SUBMIT_DIR}/temp${SLURM_JOB_ID}.err"