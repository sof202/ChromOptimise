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
#SBATCH --mail-type=END 
#SBATCH --output=temp%j.log
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
Dependencies: Java, ChromHMM, gzip
Inputs:
\$1 -> Full/relative file path for configuation file directory
================================================================================
EOF
    exit 0
}

if [[ $# -eq 0 ]]; then usage; fi

## ============ ##
##    SET UP    ##
## ============ ##

configuration_directory=$1

source "${configuration_directory}/Config.txt" || \
{ echo "The configuration file does not exist in the specified location: \
${configuration_directory}"; exit 1; }

source "${WRAPPER_SCRIPT}" || exit 1

## =============== ##
##    VARIABLES    ##
## =============== ##

## ====== DEFAULTS =============================================================
if ! [[ "${BIN_SIZE}" =~ ^[0-9]+$ ]]; then
    BIN_SIZE=200
    echo "Invalid bin size given, using the default value of ${BIN_SIZE}" \
    "instead."
fi


if [[ -z "${ASSEMBLY}" ]]; then
    ASSEMBLY=hg19
    echo "No assembly was given, using the default value of" \
    "${ASSEMBLY} instead."
fi
## =============================================================================

## ================================= ##
##    CREATE CELL MARK FILE TABLE    ##
## ================================= ##

bam_CellMarkFileTable="${BINARY_DIR}/bed_cellmarkfiletable.txt"
rm -f "${bam_CellMarkFileTable}" 
for file in "${INPUT_DIRECTORY}"/*.bam; do
    echo -ne "${CELL_TYPE}\t" >> \
    "${bam_CellMarkFileTable}"

    # Assumed that mark names are the base file names
    file_name=$(basename "$file") 
    mark_name=${file_name%.*}

    echo -ne "${mark_name}\t" >> "${bam_CellMarkFileTable}"
    echo "$file" >> "${bam_CellMarkFileTable}"
done

bed_CellMarkFileTable="${BINARY_DIR}/bed_cellmarkfiletable.txt"
rm -f "${bed_CellMarkFileTable}" 
for file in "${INPUT_DIRECTORY}"/*.bed; do
    echo -ne "${CELL_TYPE}\t" >> \
    "${bed_CellMarkFileTable}"

    file_name=$(basename "$file") 
    mark_name=${file_name%.*}

    echo -ne "${mark_name}\t" >> "${bed_CellMarkFileTable}"
    echo "$file" >> "${bed_CellMarkFileTable}"
done

## ================================= ##
##    BINARIZATION USING CHROMHMM    ##
## ================================= ##

rm -rf "${BINARY_DIR}/BinSize_${BIN_SIZE}"
mkdir -p "${BINARY_DIR}/BinSize_${BIN_SIZE}"

module purge
module load Java

if [[ -s "${bam_CellMarkFileTable}" ]]; then
    echo "Binarizing bam files found in: "\
    "${INPUT_DIRECTORY} using a bin size of: ${BIN_SIZE}."

    java \
        -mx4G \
        -jar "${CHROMHMM_MAIN_DIR}/ChromHMM.jar" \
        BinarizeBam \
        -b "${BIN_SIZE}" \
        -gzip \
        "${CHROMHMM_CHROM_SIZES}/${ASSEMBLY}.txt" \
        "${INPUT_DIRECTORY}" \
        "${BINARY_DIR}/bam_cellmarkfiletable.txt" \
        "${BINARY_DIR}/BinSize_${BIN_SIZE}/bam"
fi

if [[ -s "${bed_CellMarkFileTable}" ]]; then
    echo "Binarizing bed files found in: "\
    "${INPUT_DIRECTORY} using a bin size of: ${BIN_SIZE}."

    java \
        -mx4G \
        -jar "${CHROMHMM_MAIN_DIR}/ChromHMM.jar" \
        BinarizeBed \
        -b "${BIN_SIZE}" \
        -gzip \
        "${CHROMHMM_CHROM_SIZES}/${ASSEMBLY}.txt" \
        "${INPUT_DIRECTORY}" \
        "${BINARY_DIR}/bed_cellmarkfiletable.txt" \
        "${BINARY_DIR}/BinSize_${BIN_SIZE}/bed"
fi

java \
    -mx4G \
    -jar "${CHROMHMM_MAIN_DIR}/ChromHMM.jar" \
    MergeBinary \
    -gzip \
    "${BINARY_DIR}/BinSize_${BIN_SIZE}" \
    "${BINARY_DIR}/BinSize_${BIN_SIZE}"

## ============ ##
##   CLEAN UP   ##
## ============ ##

rm -r "${BINARY_DIR}/BinSize_${BIN_SIZE}/bam" \
    "${BINARY_DIR}/BinSize_${BIN_SIZE}/bed" \
    "${bam_CellMarkFileTable}" "${bed_CellMarkFileTable}"

# Non-autosomal chromosomes are not factored into ldsc step of the pipeline
# so we minimise their impact by deleting their associated binary files
find "${BINARY_DIR}/BinSize_${BIN_SIZE}" \
    -type f \
    ! -name "*_chr[0-9]_binary.txt.gz" \
    -a ! -name "*_chr[0-2][0-9]_binary.txt.gz" \
    -exec rm {} \;

finishing_statement 0
