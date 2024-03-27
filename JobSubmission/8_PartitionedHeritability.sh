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
#SBATCH --job-name=8_Heritability

## ===========================================================================##
##                                                                            ||
##                                  PREAMBLE                                  ||
##                                                                            ||
## ===========================================================================##
## PURPOSE:                                                                   ||
## This script will use ldsc to calculate the partitioned heritability for    ##
## each gwas trait. It is automatically called by 7_ReferenceLDSCore.sh.      ##
## After generating the results files, heatmaps of enrichment and bar plots   ##
## for the p-values of these enrichments                                      ##
## ===========================================================================##
## AUTHOR: Sam Fletcher                                                       ||
## CONTACT: s.o.fletcher@exeter.ac.uk                                         ||
## CREATED: March 2023                                                        ||
## ===========================================================================##
## PREREQUISITES: Run 7_ReferenceLDSCore.sh                                   ||
## ===========================================================================##
## DEPENDENCIES: R, LDSC, gwas traits (BED files), conda                      ||
##               1000 genomes files (plink files, weights)                    ||
## ===========================================================================##
## INPUTS:                                                                    ||
## -c|--config=     -> Full/relative file path for configuation file directory||
## -b|--binsize=    -> The bin size used in 4_BinarizeBamFiles                ||
## -s|--samplesize= -> The sample size used in 3_SubsampleBamFiles            ||
## -n|--nummodels=  -> Number of models to learn (default: 4)                 ||
## ===========================================================================##
## OUTPUTS:                                                                   ||
## Results files for paritioned heritability for each GWAS trait              ||
## Heatmap of enrichments from partitioned heritability for each GWAS trait   ||
## Bar plots for enrichment p values for each GWAS trait                      ||
## ===========================================================================##


## ===================== ##
##   ARGUMENT PARSING    ##
## ===================== ##

usage() {
cat <<EOF
===========================================================================
8_PartitionedHeritability
===========================================================================
Purpose: Determines and plots partitioned heritability using LDSC
Author: Sam Fletcher
Contact: s.o.fletcher@exeter.ac.uk
Dependencies: R, LDSC, gwas traits, 1000 genomes files
Inputs:
-c|--config=     -> Full/relative file path for configuation file directory
-b|--binsize=    -> The bin size used in 4_BinarizeBamFiles
-s|--samplesize= -> The sample size used in 3_SubsampleBamFiles
-n|--nummodels=  -> Number of models to learn (default: 4)
===========================================================================
EOF
    exit 0
}

needs_argument() {
    # Required check in case user uses -a -b or -b -a (no argument given).
    if [[ -z "$OPTARG" || "${OPTARG:0:1}" == - ]]; then usage; fi
}

if [[ ! $1 =~ -.* ]]; then usage; fi

while getopts c:b:s:n:-: OPT; do
    # Adds support for long options by reformulating OPT and OPTARG
    # This assumes that long options are in the form: "--long=option"
    if [ "$OPT" = "-" ]; then
        OPT="${OPTARG%%=*}"
        OPTARG="${OPTARG#"$OPT"}"
        OPTARG="${OPTARG#=}"
    fi
    case "$OPT" in
        c | config )      needs_argument; configuration_directory="$OPTARG" ;;
        b | binsize )     needs_argument; bin_size="$OPTARG" ;;
        s | samplesize )  needs_argument; sample_size="$OPTARG" ;;
        n | nummodels )   needs_argument; number_of_models="$OPTARG" ;;
        \? )              usage ;;  # Illegal short options are caught by getopts
        * )               usage ;;  # Illegal long option
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
echo "Configuration file used with Rscripts: \
${configuration_directory}/config.R"
echo ""
cat "${configuration_directory}/config.R"
echo ""

source "${configuration_directory}/LogFileManagement.sh" || \
{ echo "The log file management script does not exist in the specified \
location: ${configuration_directory}"; exit 1; }


# Temporary log files are moved like this as SLURM cannot create directories.
# The alternative would be forcing the user to create the file structure
# themselves and using full file paths in the SLURM directives (bad)
mv "${SLURM_SUBMIT_DIR}/temp${SLURM_JOB_ID}.log" \
"${LOG_FILE_PATH}/${SLURM_JOB_ID}~${timestamp:=}.log"
mv "${SLURM_SUBMIT_DIR}/temp${SLURM_JOB_ID}.err" \
"${LOG_FILE_PATH}/${SLURM_JOB_ID}~$timestamp.err"


## =============== ##
##    VARIABLES    ##
## =============== ##

<<<<<<< HEAD
=======
if [[ -z "${gwas_pattern}" || "${gwas_pattern}" == "'*'" ]]; then
    gwas_pattern=""
    echo "No glob pattern was given for selecting gwas traits." \
    "All gwas traits will be considered."
fi

>>>>>>> 4748ca36ac08b1e0cae0911a11658a5da26cfdaa
if [[ -z "${bin_size}" || -z "${sample_size}" || -z "${number_of_models}" ]]; then
    # If the user doesn't put in all of these options, our best hope is to look
    # for the first approximate match
    ld_directory=$( \
    find "${LD_ASSESSMENT_DIR}" -type d \
    -name "BinSize_*${bin_size}*_SampleSize_*${sample_size}*_*${number_of_models}*" | \
    head -1)

    bin_size=$(basename "${ld_directory}" | cut -d_ -f2)
    sample_size=$(basename "${ld_directory}" | cut -d_ -f4)
    number_of_models=$(basename "${ld_directory}" | cut -d_ -f5)
else
ld_directory="${LD_ASSESSMENT_DIR}\
/BinSize_${bin_size}_SampleSize_${sample_size}_${number_of_models}"
fi

## ============================ ##
##   PARTITIONED HERITABILITY   ##
## ============================ ##

module purge
module load Anaconda3/2020.02

source "${CONDA_SHELL}/profile.d/conda.sh" || \
{ echo "profile.d/conda.sh does not exist in specified location: \
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

    # Despite the fact that there are no overlapping annotations, we still
    # need to parse frq files as otherwise ldsc does not output .results files
    python \
    "${LD_SOFTWARE_DIR}/ldsc.py" \
    --h2          "${file_name}" \
    --ref-ld-chr  "${ld_directory}/annotation/ChromOptimise." \
    --w-ld-chr    "${LD_WEIGHTS_DIR}/${weights_prefix}" \
    --frqfile-chr "${LD_FRQ_DIR}/${frq_prefix}" \
    --overlap-annot \
    --out         "${ld_directory}/heritability/${output_file}"
done

## ====================== ##
##   DATA VISUALISATION   ##
## ====================== ##

conda deactivate
module purge
module load R/4.2.1-foss-2022a

cd "${RSCRIPTS_DIR}" || \
{ >&2 echo "ERROR: [\${RSCRIPTS_DIR} - ${RSCRIPTS_DIR}] doesn't exist, \
make sure FilePaths.txt is pointing to the correct directory"
finishing_statement 1; }

Rscript HeritabilityPlots.R \
<(find "${ld_directory}/heritability" -name "*.results") \
"${ld_directory}/plots"
