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
#SBATCH --job-name=1_Binarization

usage() {
cat <<EOF
================================================================================
1_BinarizeBamFiles
================================================================================
Purpose: Creates a 'cell mark file table' and uses ChromHMM's BinarizeBam 
command to binarize bam files.
Author: Sam Fletcher
Contact: s.o.fletcher@exeter.ac.uk
Dependencies: Java, ChromHMM
Inputs:
-c|--config=    -> Full/relative file path for configuation file directory
-i|--input=     -> Input directory holding bed/bam files
-b|--binsize=   -> Bin size to be used by BinarizeBam command (default: 200)
-a|--assembly=  -> The assembly to use (default: hg19)
================================================================================
EOF
    exit 0
}

needs_argument() {
    # Required check in case user uses -a -b or -b -a (no argument given).
    if [[ -z "$OPTARG" || "${OPTARG:0:1}" == - ]]; then usage; fi
}

if [[ ! $1 =~ -.* ]]; then usage; fi

while getopts c:b:i:a:-: OPT; do
    # Adds support for long options by reformulating OPT and OPTARG
    # This assumes that long options are in the form: "--long=option"
    if [ "$OPT" = "-" ]; then
        OPT="${OPTARG%%=*}"
        OPTARG="${OPTARG#"$OPT"}"
        OPTARG="${OPTARG#=}"
    fi
    case "$OPT" in
        c | config )       needs_argument; configuration_directory="$OPTARG" ;;
        i | input )        needs_argument; input_directory="$OPTARG" ;;
        b | binsize )      needs_argument; bin_size="$OPTARG" ;;
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

source "${WRAPPER_SCRIPT}" || \
{ echo "The log file management script does not exist in the specified \
location: ${configuration_directory}"; exit 1; }

## =============== ##
##    VARIABLES    ##
## =============== ##

## ====== DEFAULTS =============================================================
if ! [[ "${bin_size}" =~ ^[0-9]+$ ]]; then
    bin_size=200
    echo "Invalid bin size given, using the default value of ${bin_size}" \
    "instead."
fi


if [[ -z "${assembly}" ]]; then
    assembly=hg19
    echo "No assembly was given, using the default value of" \
    "${assembly} instead."
fi
## =============================================================================

## ================================= ##
##    CREATE CELL MARK FILE TABLE    ##
## ================================= ##

# As we have already merged files in 3_SubsampleBamFiles.sh, the cell mark
# file table will just be 1 file for each mark that has been processed. The
# only added level of complexity is obtaining the mark name from the file names.

bam_CellMarkFileTable="${BINARY_DIR}/bed_cellmarkfiletable.txt"
rm -f "${bam_CellMarkFileTable}" 
for file in "${input_directory}"/*.bam; do
    # We're assuming here that there is only one cell type inspected
    echo -ne "${CELL_TYPE}\t" >> \
    "${bam_CellMarkFileTable}"
    # The subsampled files are named: subsampled.[SampleSize].[mark_name].bam. 
    # Below extracts the mark name
    mark_name=$(echo "$file" | cut -d "." -f 3) 
    echo -ne "${mark_name}\t" >> "${bam_CellMarkFileTable}"
    echo "$file" >> "${bam_CellMarkFileTable}"
done

bed_CellMarkFileTable="${BINARY_DIR}/bed_cellmarkfiletable.txt"
rm -f "${bed_CellMarkFileTable}" 
for file in "${input_directory}"/*.bed; do
    echo -ne "${CELL_TYPE}\t" >> \
    "${bed_CellMarkFileTable}"
    mark_name=$(echo "$file" | cut -d "." -f 3) 
    echo -ne "${mark_name}\t" >> "${bed_CellMarkFileTable}"
    echo "$file" >> "${bed_CellMarkFileTable}"
done

## ================================= ##
##    BINARIZATION USING CHROMHMM    ##
## ================================= ##

rm -rf "${BINARY_DIR}/BinSize_${bin_size}"
mkdir -p "${BINARY_DIR}/BinSize_${bin_size}"

module purge
module load Java

if [[ -s "${bam_CellMarkFileTable}" ]]; then
    echo "Binarizing bam files found in: "\
    "${input_directory} using a bin size of: ${bin_size}."

    java -mx4G \
    -jar "${CHROMHMM_MAIN_DIR}/ChromHMM.jar" \
    BinarizeBam \
    -b "${bin_size}" \
    -gzip \
    "${CHROMHMM_CHROM_SIZES}/${assembly}.txt" \
    "${input_directory}" \
    "${BINARY_DIR}/bam_cellmarkfiletable.txt" \
    "${BINARY_DIR}/BinSize_${bin_size}/bam"
fi

if [[ -s "${SUBSAMPLED_DIR}/bed_cellmarkfiletable.txt" ]]; then
    echo "Binarizing bed files found in: "\
    "${input_directory} using a bin size of: ${bin_size}."

    java -mx4G \
    -jar "${CHROMHMM_MAIN_DIR}/ChromHMM.jar" \
    BinarizeBed \
    -b "${bin_size}" \
    -gzip \
    "${CHROMHMM_CHROM_SIZES}/${assembly}.txt" \
    "${input_directory}" \
    "${BINARY_DIR}/bam_cellmarkfiletable.txt" \
    "${BINARY_DIR}/BinSize_${bin_size}/bed"
fi

java -mx4G \
-jar "${CHROMHMM_MAIN_DIR}/ChromHMM.jar" \
MergeBinary \
-gzip \
"${BINARY_DIR}/BinSize_${bin_size}" \
"${BINARY_DIR}/BinSize_${bin_size}"

## ============ ##
##   CLEAN UP   ##
## ============ ##

rm -r "${BINARY_DIR}/BinSize_${bin_size}/bam" \
    "${BINARY_DIR}/BinSize_${bin_size}/bed" \
    "${bam_CellMarkFileTable}" "${bed_CellMarkFileTable}"

# Non-autosomal chromosomes are not factored into ldsc step of the pipeline
# so we minimise their impact by deleting their associated binary files
find "${BINARY_DIR}/BinSize_${bin_size}" \
-type f \
! -name "ChromOptimise_chr[0-9]_binary.txt.gz" \
-a ! -name "ChromOptimise_chr[0-9][0-9]_binary.txt.gz" \
-exec rm {} \;

finishing_statement 0
