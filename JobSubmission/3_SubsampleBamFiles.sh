#!/bin/bash
# Export all environment variables to the batch job
#SBATCH --export=ALL
# Submit to the mrc queue for faster queue times
#SBATCH -p mrcq 
# Consult information/Processing_Times.md for expected time
#SBATCH --time=24:00:00 
#SBATCH -A Research_Project-MRC190311 
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=16 
# Predicted that memory consumption will rise massively when merging lots of files
#SBATCH --mem=10G 
# Send an email after the job is done
#SBATCH --mail-type=END 
# Temporary log file, later to be removed
#SBATCH --output=temp%j.log
# Temporary error file, later to be removed
#SBATCH --error=temp%j.err
#SBATCH --job-name=3_Merging_and_Subsampling

## =================================================================================##
##                                                                                  ||
##                                     PREAMBLE                                     ||
##                                                                                  ||
## =================================================================================##
## PURPOSE:                                                                         ||
## Obtain a sample of the bam files. The bam files have varying sizes due to the    ||
## number of reads. Sampling the files produced in 2_ProcessBamFiles.sh will lead   ||
## to samples with the same number of files that contain a different number of      ||
## reads. This removes the reproducability of the proceedure. To get around         ||
## this, this script merges all of the processed .bam files and subsequently        ||
## samples this larger file randomly. To save space, the merged file is deleted.    ||
## =================================================================================##
## AUTHOR: Sam Fletcher                                                             ||
## CONTACT: s.o.fletcher@exeter.ac.uk                                               ||
## CREATED: November 2023                                                           ||
## =================================================================================##
## PREREQUISITES: Run 2_batch_ProcessBamFiles.sh                                    ||
## =================================================================================##
## DEPENDENCIES: Samtools                                                           ||
## =================================================================================##
## INPUTS:                                                                          ||
## -c|--config= -> Full (or relative) file path for configuation file directory     ||
## -m|--mark= -> Epigenetic mark to process                                         ||
## -s|--samplesize= -> Sample size as a percentage (default : 50)                   ||
## =================================================================================##
## OUTPUTS:                                                                         ||
## Subsampled .bam files                                                            ||
## =================================================================================##

## ===================== ##
##   ARGUMENT PARSING    ##
## ===================== ##

usage() {
cat <<EOF
================================================================================
3_SubsampleBamFiles.sh
================================================================================
Purpose: Merges and subsamples processed .bam files present in specified folder.
Author: Sam Fletcher
Contact: s.o.fletcher@exeter.ac.uk
Dependencies: Samtools
Inputs:
-c|--config=     -> Full/relative file path for configuation file directory
-m|--mark=       -> Epigenetic mark to process
-s|--samplesize= -> Sample size as a percentage (default : 50) 
================================================================================
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
        c | config )       needs_argument; configuration_directory="$OPTARG" ;;
        m | mark )         needs_argument; mark_name="$OPTARG" ;;
        s | samplesize )   needs_argument; sample_size="$OPTARG" ;;
        \? )               usage ;;  # Illegal short options are caught by getopts
        * )                usage ;;  # bad long option
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
# [epigenetic mark name]~[Sample size]~[job id]~[date]-[time]

mv "${SLURM_SUBMIT_DIR}/temp${SLURM_JOB_ID}.log" \
"${LOG_FILE_PATH}/$2~$3~${SLURM_JOB_ID}~${timestamp:=}.log"
mv "${SLURM_SUBMIT_DIR}/temp${SLURM_JOB_ID}.err" \
"${LOG_FILE_PATH}/$2~$3~${SLURM_JOB_ID}~$timestamp.err"

## =============== ##
##    VARIABLES    ##
## =============== ##

PROCESSED_FULL_FILE_PATH="${PROCESSED_DIR}/${mark_name}"

if [[ -z "${mark_name}" ]]; then
    { >&2 echo -e "ERROR: No epigenetic mark name given.\n\
    Ensure that the first argument is the name of a processed epigenetic mark." ; \
    finishing_statement 1; }
fi

## ====== DEFAULTS ====================================================================
if ! [[ "${sample_size}" =~ ^[0-9]+$ ]]; then
    sample_size=50
    echo "Invalid sample size was given, using the default value of: ${sample_size}%."
fi

if [[ "${sample_size}" -gt 100 || "${sample_size}" -le 0 ]]; then
    echo "Sample size must be greater than 0 and less than or equal to 100."
    sample_size=50
    echo "Using the default value of: ${sample_size}% instead."
fi
## ====================================================================================

## ========================= ##
##   MERGING OF .BAM FILES   ##
## ========================= ##

cd "${PROCESSED_FULL_FILE_PATH}" || \
{ >&2 echo "ERROR: [\${PROCESSED_DIR}/\${mark_name} - ${PROCESSED_DIR}/${mark_name}] \
doesn't exist, make sure that you typed the epigenetic mark correctly and that \
FilePaths.txt is pointing to the correct directory."; finishing_statement 1; }

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
{ >&2 echo "ERROR: [\${SUBSAMPLED_DIR} - ${SUBSAMPLED_DIR} doesn't exist,] \
make sure FilePaths.txt is pointing to the correct directory"; finishing_statement 1; }

sample_size_decimal=$(echo "scale=2; $sample_size /100" | bc)
echo "Subsampling merged .bam file with sample size ${sample_size}%..."

# Ensure headers are kept in subsampled file to avoid errors later in pipeline
samtools view -s "${sample_size_decimal}" -h "${output_file_path}" \
> "Subsampled.${sample_size}.${mark_name}.bam"

rm "FullMerged.${mark_name}.bam"
cd "${PROCESSED_FULL_FILE_PATH}" || finishing_statement 1
rm List_Of_Bam_Files_To_Merge.txt

finishing_statement 0