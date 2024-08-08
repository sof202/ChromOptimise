#!/bin/bash
# Export all environment variables to the batch job
#SBATCH --export=ALL 
# Submit to the mrc queue for faster queue times
#SBATCH -p mrcq 
# The time of this script increases linearly with the number of gwas traits
# considered
#SBATCH --time=24:00:00
#SBATCH -A Research_Project-MRC190311 
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=16
# According to ldsc, roughly 8 GB are required for 50 categories, we increase
# this here as we have 50 baseline categories plus each of our states
#SBATCH --mem=50G 
# Send an email after the job is done
#SBATCH --mail-type=END 
# Temporary log file, later to be removed
#SBATCH --output=temp%j.log
# Temporary error file, later to be removed
#SBATCH --error=temp%j.err
#SBATCH --job-name=5_Heritability

usage() {
cat <<EOF
===========================================================================
5_PartitionedHeritability
===========================================================================
Purpose: Determines and plots partitioned heritability using LDSC
Author: Sam Fletcher
Contact: s.o.fletcher@exeter.ac.uk
Dependencies: R, LDSC, gwas traits, 1000 genomes files
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

ld_directory="${LD_ASSESSMENT_DIR}/BinSize_${BIN_SIZE}_models_${NUMBER_OF_MODELS}"

## ============================ ##
##   PARTITIONED HERITABILITY   ##
## ============================ ##

rm -rf \
    "${ld_directory}/heritability" \
    "${ld_directory}/plots/Chromoptimise_Categories" \
    "${ld_directory}/plots/All_Categories" 

mkdir -p \
    "${ld_directory}/heritability" \
    "${ld_directory}/plots/Chromoptimise_Categories" \
    "${ld_directory}/plots/All_Categories" 

source "${CONDA_SHELL}" || \
{ echo "conda.sh does not exist in specified location: \
[\${CONDA_SHELL} - ${CONDA_SHELL}]"; exit 1; }
conda activate "${LDSC_ENVIRONMENT}"

weights_prefix=$(\
find "${LD_WEIGHTS_DIR}" -type f -name "*22.l2*" -print0 | \
xargs -0 basename | \
sed "s/22\..*//" \
)

frq_prefix=$(\
find "${LD_FRQ_DIR}" -type f -name "*22.frq*" -print0 | \
xargs -0 basename | \
sed "s/22\..*//" \
)

gwas_traits=$(\
find "${LD_GWAS_TRAITS_DIR}" -name "*.sumstats*"\
)

for file_name in ${gwas_traits}; do
    output_file=$(basename "${file_name}" .sumstats.gz)

    python \
        "${LD_SOFTWARE_DIR}/ldsc.py" \
        --h2          "${file_name}" \
        --ref-ld-chr  "${ld_directory}/annotation/ChromOptimise." \
        --w-ld-chr    "${LD_WEIGHTS_DIR}/${weights_prefix}" \
        --frqfile-chr "${LD_FRQ_DIR}/${frq_prefix}" \
        --overlap-annot \
        --out         "${ld_directory}/heritability/${output_file}"
done

## ================= ##
##   VISUALISATION   ##
## ================= ##

conda deactivate
module purge
module load R/4.2.1-foss-2022a

Rscript "${RSCRIPTS_DIR}/HeritabilityPlots.R" \
    <(find "${ld_directory}/heritability" -name "*.results") \
    "${CELL_TYPE}" \
    "${ld_directory}/plots"
