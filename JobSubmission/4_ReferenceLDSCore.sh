#!/bin/bash
# Export all environment variables to the batch job
#SBATCH --export=ALL 
# Submit to the mrc queue for faster queue times
#SBATCH -p mrcq 
# SNP assignment takes ~20-30 minutes for chr 1, ldsc takes <15 minutes
#SBATCH --time=2:00:00
#SBATCH -A Research_Project-MRC190311 
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=16
#SBATCH --array=1-22
# memory consumption of SNPassignment is somewhere around 1GB per array element
# for ldsc, it can be much higher when using lots of categories
#SBATCH --mem=100G 
# Send an email after the job is done
#SBATCH --mail-type=END 
# Temporary log file, later to be removed
#SBATCH --output=temp%A_%a.log
# Temporary error file, later to be removed
#SBATCH --error=temp%A_%a.err
#SBATCH --job-name=4_ReferenceLDSCore

usage() {
cat <<EOF
===========================================================================
4_ReferenceLDSCore
===========================================================================
Purpose: Generates annotation files based on baseline annotations (LDSC)
and ChromHMM state annotations. Obtains LDSCores from said annotation.
Author: Sam Fletcher
Contact: s.o.fletcher@exeter.ac.uk
Dependencies: R, Bedtools, LDSC, gwas traits, 1000 genomes files
Inputs:
\$1 -> Full/relative file path for configuation file directory
===========================================================================
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

input_directory="${OPTIMUM_STATES_DIR}/BinSize_${BIN_SIZE}_models_${NUMBER_OF_MODELS}"

if [[ -z "${OPTIMUM_NUMBER_OF_STATES}" ]]; then
    if [[ -z "$(ls -A "${input_directory}")" ]]; then
        { >&2 echo -e "ERROR: No files found in: ${input_directory}.\n"\
        "Please run 3_OptimumNumberOfStates.sh before this script."
        finishing_statement 1; }
    fi
    optimum_state_file="${input_directory}/OptimumNumberOfStates.txt"
    OPTIMUM_NUMBER_OF_STATES=$(tail -1 "${optimum_state_file}" | \
    cut -d: -f2 | \
    tr -d ' ')
    echo "Optimum model size found was: ${OPTIMUM_NUMBER_OF_STATES}"
fi

## =================== ##
##   FILE MANAGEMENT   ##
## =================== ##

output_directory="${LD_ASSESSMENT_DIR}\
/BinSize_${BIN_SIZE}_models_${NUMBER_OF_MODELS}"

if [[ "${SLURM_ARRAY_TASK_ID}" -eq 1 ]]; then
    rm -rf "${output_directory:?}"
    mkdir -p \
        "${output_directory}/annotation" \
        "${output_directory}/heritability" \
        "${output_directory}/plots/ChromOptimise_Categories" \
        "${output_directory}/plots/All_Categories" 
fi

# We sleep here to ensure files are not removed prematurely
sleep 5

full_model_directory="${MODEL_DIR}/BinSize_${BIN_SIZE}_models_${NUMBER_OF_MODELS}"

full_binary_directory="${BINARY_DIR}/BinSize_${BIN_SIZE}"

# We ignore non-autosomal chromosomes as 1000 genomes doesn't provide this data
# Hence our chromosomes are just 1-22 (the array indices)
chromosome=${SLURM_ARRAY_TASK_ID}

# The production of the annotation file requires some intermediary files
# so we create a temporary directory to hold these
temporary_directory="${output_directory}/temp_${chromosome}"
mkdir -p "${temporary_directory}"

dense_bed_file=$(\
find "${full_model_directory}" \
-name "*_${OPTIMUM_NUMBER_OF_STATES}_dense.bed")

binary_file=$(\
    find "${full_binary_directory}" \
-name "*_chr${chromosome}_binary*")

bim_file=$(\
    find "${LD_PLINK_DIR}" \
-name "*.${chromosome}.bim")

## ============================ ##
##   ANNOTATION FILE CREATION   ##
## ============================ ##

cd "${RSCRIPTS_DIR}" || \
{ >&2 echo "ERROR: [\${RSCRIPTS_DIR} - ${RSCRIPTS_DIR}] doesn't exist, \
make sure FilePaths.txt is pointing to the correct directory"
finishing_statement 1; }

## ---------------------- ##
##  CONVERT TO BED FORMAT ## 
## ---------------------- ##

module purge
module load R/4.2.1-foss-2022a

Rscript BinarytoBed.R \
    <(zcat "${binary_file}") \
    "${BIN_SIZE}" \
    "chr${chromosome}" \
    "${temporary_directory}/binary-${chromosome}.bed"

Rscript BimtoBed.R \
    <(cat "${bim_file}") \
    "${temporary_directory}/SNP_positions-${chromosome}.bed"

## ----------------------------------- ##
##   FIND STATE AND MARK ASSIGNMENTS   ##
## ----------------------------------- ##

module purge
module load BEDTools/2.29.2-GCC-9.3.0

bedtools intersect -wb \
    -a "${temporary_directory}/SNP_positions-${chromosome}.bed" \
    -b "${dense_bed_file}" | \
    awk '{print $7}' > \
    "${temporary_directory}/state_assignments-${chromosome}.txt"

# We get the mark names at the top of the file for the Rscript that appends
# these columns to the annotation file later for convenience
zcat "${binary_file}" | \
    awk 'NR==2 {for(i=1; i<=NF; i++) \
    printf "CELL_TYPE_%s%s", $i, (i==NF ? "\n" : "\t")}' > \
    "${temporary_directory}/mark_assignments-${chromosome}.txt"

sed -i "s/CELL_TYPE/${CELL_TYPE}/g" \
    "${temporary_directory}/mark_assignments-${chromosome}.txt"

bedtools intersect -wb \
    -a "${temporary_directory}/SNP_positions-${chromosome}.bed" \
    -b "${temporary_directory}/binary-${chromosome}.bed" | \
    awk '{ for (i=7; i<=NF; i++) printf "%s%s", $i, (i<NF ? "\t" : "\n") }' >> \
    "${temporary_directory}/mark_assignments-${chromosome}.txt"

## ----------------------- ##
##   GENERATE ANNOTATION   ##
## ----------------------- ##

module purge
module load R/4.2.1-foss-2022a

baseline_annot="${LD_BASELINE_DIR}/baselineLD.${chromosome}.annot.gz"

Rscript CreateAnnotationFile.R \
    <(zcat "${baseline_annot}") \
    <(cat "${temporary_directory}/state_assignments-${chromosome}.txt") \
    <(cat "${temporary_directory}/mark_assignments-${chromosome}.txt") \
    "${OPTIMUM_NUMBER_OF_STATES}" \
    "${output_directory}/annotation/ChromOptimise.${chromosome}.annot"
    "${CELL_TYPE}"

rm -rf "${temporary_directory}"

module purge

## ======================= ##
##   REFERENCE LD SCORES   ##
## ======================= ##

source "${CONDA_SHELL}" || \
{ echo "conda.sh does not exist in specified location: \
[\${CONDA_SHELL} - ${CONDA_SHELL}]"; exit 1; }
conda activate "${LDSC_ENVIRONMENT}"

plink_prefix=$(\
find "${LD_PLINK_DIR}" -type f -name "*22.bim" -print0 | \
xargs -0 basename | \
sed "s/22\..*//" \
)

python \
    "${LD_SOFTWARE_DIR}/ldsc.py" \
    --l2 \
    --bfile      "${LD_PLINK_DIR}/${plink_prefix}${chromosome}" \
    --ld-wind-cm 1 \
    --annot      "${output_directory}/annotation/ChromOptimise.${chromosome}.annot" \
    --out        "${output_directory}/annotation/ChromOptimise.${chromosome}"

conda deactivate

## ============================ ##
##   PARTITIONED HERITABILITY   ##
## ============================ ##

# This job is being ran as an array, which means that the memory of the job
# is split among each array task. This memory allocation is not dynamic and so
# at this point in the program the task only has (max memory)/22 GB of memory.
# This is not enough to handle the partitioned heritability (unless you have
# ~220GB of memory available). Hence at this point we run a new script with
# sbatch

if [[ ${SLURM_ARRAY_TASK_ID} -eq 1 ]]; then
    cd "${JOBSUBMISSION_DIR}" || \
    { echo "Could not find the JobSubmission directory in \
    ${JOBSUBMISSION_DIR}/Jobsubmission. Please check your configuration file." \
    finishing_statement 0; }

    sbatch \
        --dependency=afterok:"${SLURM_ARRAY_JOB_ID}" \
        5_PartitionedHeritability.sh \
        "${configuration_directory}"

    finishing_statement 0
else
    finishing_statement 0
fi
