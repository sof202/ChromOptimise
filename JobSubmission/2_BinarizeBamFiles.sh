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
#SBATCH --job-name=2_Binarization

## ===========================================================================##
##                                                                            ||
##                                  PREAMBLE                                  ||
##                                                                            ||
## ===========================================================================##
## PURPOSE:                                                                   ||
## The subsampled files for each mark will now be binarized through the use   ||
## of ChromHMM's BinarizeBam command. This script is to be ran after all of   ||
## the epigenetic marks that one wants to inspect have been subsampled. The   ||
## .bam files need to be binarized before they can be used by ChromHMM's      ||
## LearnModel command.                                                        ||
## ===========================================================================##
## AUTHOR: Sam Fletcher                                                       ||
## CONTACT: s.o.fletcher@exeter.ac.uk                                         ||
## CREATED: November 2023                                                     ||
## ===========================================================================##
## PREREQUISITES: Run: 1_SubsampleBamFiles.sh                                 ||
## ===========================================================================##
## DEPENDENCIES: Java, ChromHMM                                               ||
## ===========================================================================##
## INPUTS:                                                                    ||
## -c|--config=     -> Full/relative file path for configuation file directory||
## -b|--binsize=    -> Bin size to be used by BinarizeBam command             ||
##                     (default: 200)                                         ||
## -s|--samplesize= -> Sample size used in 3_SubsampleBamFiles.sh             ||
## -a|--assembly=   -> The assembly to use (default: hg19)                    ||
## ===========================================================================##
## OUTPUTS:                                                                   ||
## Binary signal files for every chromosome in the dataset                    ||
## ===========================================================================##

## ===================== ##
##   ARGUMENT PARSING    ##
## ===================== ##

usage() {
cat <<EOF
================================================================================
2_BinarizeBamFiles
================================================================================
Purpose: Creates a 'cell mark file table' and uses ChromHMM's BinarizeBam 
command to binarize bam files.
Author: Sam Fletcher
Contact: s.o.fletcher@exeter.ac.uk
Dependencies: Java, ChromHMM
Inputs:
-c|--config=     -> Full/relative file path for configuation file directory
-b|--binsize=    -> Bin size to be used by BinarizeBam command (default: 200)
-s|--samplesize= -> Sample size used in 3_SubsampleBamFiles.sh
-a|--assembly=   -> The assembly to use (default: hg19)
================================================================================
EOF
    exit 0
}

needs_argument() {
    # Required check in case user uses -a -b or -b -a (no argument given).
    if [[ -z "$OPTARG" || "${OPTARG:0:1}" == - ]]; then usage; fi
}

if [[ ! $1 =~ -.* ]]; then usage; fi

while getopts c:b:s:a:-: OPT; do
    # Adds support for long options by reformulating OPT and OPTARG
    # This assumes that long options are in the form: "--long=option"
    if [ "$OPT" = "-" ]; then
        OPT="${OPTARG%%=*}"
        OPTARG="${OPTARG#"$OPT"}"
        OPTARG="${OPTARG#=}"
    fi
    case "$OPT" in
        c | config )       needs_argument; configuration_directory="$OPTARG" ;;
        b | binsize )      needs_argument; bin_size="$OPTARG" ;;
        s | samplesize )   needs_argument; sample_size="$OPTARG" ;;
        a | assembly )     needs_argument; assembly="$OPTARG" ;;
        \? )               usage ;;  # Illegal short options are caught by getopts
        * )                usage ;;  # Illegal long option
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


# Temporary log files are moved like this as SLURM cannot create directories.
# The alternative would be forcing the user to create the file structure
# themselves and using full file paths in the SLURM directives (bad)
mv "${SLURM_SUBMIT_DIR}/temp${SLURM_JOB_ID}.log" \
"${LOG_FILE_PATH}/BinSize-${bin_size:=200}~${SLURM_JOB_ID}~${timestamp:=}.log"
mv "${SLURM_SUBMIT_DIR}/temp${SLURM_JOB_ID}.err" \
"${LOG_FILE_PATH}/BinSize-${bin_size:=200}~${SLURM_JOB_ID}~$timestamp.err"

## =============== ##
##    VARIABLES    ##
## =============== ##

## ====== DEFAULTS =============================================================
if ! [[ "${bin_size}" =~ ^[0-9]+$ ]]; then
    bin_size=200
    echo "Invalid bin size given, using the default value of ${bin_size}" \
    "instead."
fi

# 'Intelligently' find the sample size using first file name in the 
# subsampled directory
if [[ -z "${sample_size}" ]]; then
    cd "${SUBSAMPLED_DIR}" || \
    { >&2 echo "ERROR: \
    [\${SUBSAMPLED_DIR} - ${SUBSAMPLED_DIR}] doesn't exist, make "\
    "sure FilePaths.txt is pointing to the correct directory"
    finishing_statement 1; }

    sample_size=$(find . -type f -name "Subsampled*" | \
    head -1 | \
    cut -d "." -f 3)

    echo -e "WARNING: No sample size was given.\n"\
    "Assuming that ${sample_size} is the desired sample size..."
fi

# If no sample size was found, then the sampling script likely hasn't ran yet
if [[ -z "${sample_size}" ]]; then
    { >&2 echo -e "ERROR: No sample size even after fail safe. Please run "\
    "3_SubsampleBamFiles.sh before running this script."
    finishing_statement 1; }
fi

if [[ -z "${assembly}" ]]; then
    assembly=hg19
    echo "No assembly was given, using the default value of" \
    "${assembly} instead."
fi
## =============================================================================

## ================== ##
##   FILE EXISTENCE   ##
## ================== ##

cd "${SUBSAMPLED_DIR}" || \
{ >&2 echo "ERROR: \
[\${SUBSAMPLED_DIR} - ${SUBSAMPLED_DIR}] doesn't exist, make \
sure FilePaths.txt is pointing to the correct directory"
finishing_statement 1; }

if [[ -z "$(find . -type f -name "Subsampled.${sample_size}*")" ]]; then
    { >&2 echo -e "ERROR: [\${SUBSAMPLED_DIR} - ${SUBSAMPLED_DIR}] is empty.\n"\
    "Ensure that 3_SubsampleBamFiles.sh has been ran before this script."
    finishing_statement 1; }
fi

## ================================= ##
##    CREATE CELL MARK FILE TABLE    ##
## ================================= ##

# As we have already merged files in 3_SubsampleBamFiles.sh, the cell mark
# file table will just be 1 file for each mark that has been processed. The
# only added level of complexity is obtaining the mark name from the file names.

rm -f "bam_cellmarkfiletable.txt" 
if [ $(find . -name "*${sample_size}*.bam" | head -1) ]; then
    for file in *"${sample_size}"*.bam; do
        # We're assuming here that there is only one cell type inspected
        echo -ne "ChromOptimise\t" >> \
        "bam_cellmarkfiletable.txt"
        # The subsampled files are named: subsampled.[SampleSize].[mark_name].bam. 
        # Below extracts the mark name
        mark_name=$(echo "$file" | cut -d "." -f 3) 
        echo -ne "${mark_name}\t" >> "bam_cellmarkfiletable.txt"
        echo "$file" >> "bam_cellmarkfiletable.txt"
    done
fi

rm -f "bed_cellmarkfiletable.txt" 
if [ $(find . -name "*${sample_size}*.bed" | head -1) ]; then
    for file in *"${sample_size}"*.bed; do
        echo -ne "ChromOptimise\t" >> \
        "bed_cellmarkfiletable.txt"
        mark_name=$(echo "$file" | cut -d "." -f 3) 
        echo -ne "${mark_name}\t" >> "bed_cellmarkfiletable.txt"
        echo "$file" >> "bed_cellmarkfiletable.txt"
    done
fi

## ================================= ##
##    BINARIZATION USING CHROMHMM    ##
## ================================= ##
#
cd "${BINARY_DIR}" || \
{ echo "ERROR: [\${BINARY_DIR} - ${BINARY_DIR}] doesn't exist, \
make sure FilePaths.txt is pointing to the correct directory"
finishing_statement 1; }

rm -rf "BinSize_${bin_size}_SampleSize_${sample_size}"
mkdir -p "BinSize_${bin_size}_SampleSize_${sample_size}"
cd "BinSize_${bin_size}_SampleSize_${sample_size}" || finishing_statement 1

module purge
module load Java

if [[ -s "${SUBSAMPLED_DIR}/bam_cellmarkfiletable.txt" ]]; then
    echo "Binarizing subsampled bam files found in [\${SUBSAMPLED_DIR} - "\
    "${SUBSAMPLED_DIR}] with sample size: ${sample_size} using a bin size "\
    "of: ${bin_size}."

    java -mx4G \
    -jar "${CHROMHMM_MAIN_DIR}/ChromHMM.jar" \
    BinarizeBam \
    -b "${bin_size}" \
    -gzip \
    "${CHROMHMM_CHROM_SIZES}/${assembly}.txt" \
    "${SUBSAMPLED_DIR}" \
    "${SUBSAMPLED_DIR}/bam_cellmarkfiletable.txt" \
    "${BINARY_DIR}/BinSize_${bin_size}_SampleSize_${sample_size}/bam"
fi

if [[ -s "${SUBSAMPLED_DIR}/bed_cellmarkfiletable.txt" ]]; then
    echo "Binarizing subsampled bed files found in [\${SUBSAMPLED_DIR} - "\
    "${SUBSAMPLED_DIR}] with sample size: ${sample_size} using a bin size "\
    "of: ${bin_size}."

    java -mx4G \
    -jar "${CHROMHMM_MAIN_DIR}/ChromHMM.jar" \
    BinarizeBed \
    -b "${bin_size}" \
    -gzip \
    "${CHROMHMM_CHROM_SIZES}/${assembly}.txt" \
    "${SUBSAMPLED_DIR}" \
    "${SUBSAMPLED_DIR}/bed_cellmarkfiletable.txt" \
    "${BINARY_DIR}/BinSize_${bin_size}_SampleSize_${sample_size}/bed"
fi

java -mx4G \
-jar "${CHROMHMM_MAIN_DIR}/ChromHMM.jar" \
MergeBinary \
-gzip \
"${BINARY_DIR}/BinSize_${bin_size}_SampleSize_${sample_size}" \
"${BINARY_DIR}/BinSize_${bin_size}_SampleSize_${sample_size}"

rm -r "${BINARY_DIR}/BinSize_${bin_size}_SampleSize_${sample_size}/bam"
rm -r "${BINARY_DIR}/BinSize_${bin_size}_SampleSize_${sample_size}/bed"

# Non-autosomal chromosomes are not factored into ldsc step of the pipeline
# so we minimise their impact by deleting their associated binary files
find "${BINARY_DIR}/BinSize_${bin_size}_SampleSize_${sample_size}" \
-type f \
! -name "ChromOptimise_chr[0-9]_binary.txt.gz" \
-a ! -name "ChromOptimise_chr[0-9][0-9]_binary.txt.gz" \
-exec rm {} \;

finishing_statement 0
