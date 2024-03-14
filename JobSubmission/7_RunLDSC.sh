#!/bin/bash
# Export all environment variables to the batch job
#SBATCH --export=ALL 
# Submit to the mrc queue for faster queue times
#SBATCH -p mrcq 
# SNP assignment takes quite a while
#SBATCH --time=48:00:00
#SBATCH -A Research_Project-MRC190311 
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=16
# Memory consumption is very low in testing
# Consult information/Memory_Profiling.md for expected memory usage
#SBATCH --mem=1G 
# Send an email after the job is done
#SBATCH --mail-type=END 
# Temporary log file, later to be removed
#SBATCH --output=temp%j.log
# Temporary error file, later to be removed
#SBATCH --error=temp%j.err
#SBATCH --job-name=7_LDSC

## ===========================================================================##
##                                                                            ||
##                                  PREAMBLE                                  ||
##                                                                            ||
## ===========================================================================##
## PURPOSE:                                                                   ||
## Creates partitioned heritability results using ldsc (linkage               ||
## disequilibrium scores) for a selection of gwas traits. The partitions      ||
## obviously being the state assignments from ChromHMM. These results files   ||
## are then converted into heatmaps. The outputs are useful in sense checking ||
## that a state assignment is biologically relevant (and not just             ||
## statistically relevant)                                                    ||
##                                                                            ||
## In order for this to work, you will need to download 1000 genomes project  ||
## files from: https://zenodo.org/records/10515792                            ||
## ===========================================================================##
## AUTHOR: Sam Fletcher                                                       ||
## CONTACT: s.o.fletcher@exeter.ac.uk                                         ||
## CREATED: March 2023                                                        ||
## ===========================================================================##
## PREREQUISITES: Run: 5_batch_CreateIncrementalModels.sh                     ||
##                     6_OptimalNumberOfStates.sh                             ||
##                Create ldsc conda environment (follow instructions on       ||
##                GitHub page: https://github.com/bulik/ldsc)                 ||
## ===========================================================================##
## DEPENDENCIES: R, LDSC, gwas traits (BED files), conda                      ||
##               1000 genomes files (plink files, weights)                    ||
## ===========================================================================##
## INPUTS:                                                                    ||
## -c|--config=     -> Full/relative file path for configuation file directory||
## -g|--gwas=       -> The glob pattern used for selecting gwas traits to use ||
##                     for heritability analysis, your input will be wrapped  ||
##                     in "*"s                                                ||
## -t|--state=      -> (Override) Optimal number of states. Use this if you   ||
##                     want to look at the results for a different model.     ||
## -b|--binsize=    -> The bin size used in 4_BinarizeBamFiles                ||
## -s|--samplesize= -> The sample size used in 3_SubsampleBamFiles            ||
## -n|--nummodels=  -> Number of models to learn (default: 4)                 ||
## ===========================================================================##
## OUTPUTS:                                                                   ||
## Heatmap of enrichments from partitioned heritability for each GWAS trait   ||
## Bar plots for enrichment p values for each GWAS trait                      ||
## ===========================================================================##


## ===================== ##
##   ARGUMENT PARSING    ##
## ===================== ##

usage() {
cat <<EOF
===========================================================================
7_RunLDSC
===========================================================================
Purpose: Determines and plots partitioned heritability using LDSC
Author: Sam Fletcher
Contact: s.o.fletcher@exeter.ac.uk
Dependencies: R, LDSC, gwas traits, 1000 genomes files
Inputs:
-c|--config=     -> Full/relative file path for configuation file directory
-g|--gwas=       -> The glob pattern used for selecting gwas traits to use 
                    for heritability analysis, your input will be wrapped
                    in "*"s
-t|--state=      -> (Override) Optimal number of states. Use this if you 
                    want to look at the results for a different model.
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

while getopts c:t:g:b:s:n:-: OPT; do
    # Adds support for long options by reformulating OPT and OPTARG
    # This assumes that long options are in the form: "--long=option"
    if [ "$OPT" = "-" ]; then
        OPT="${OPTARG%%=*}"
        OPTARG="${OPTARG#"$OPT"}"
        OPTARG="${OPTARG#=}"
    fi
    case "$OPT" in
        c | config )     needs_argument; configuration_directory="$OPTARG" ;;
        t | state )      needs_argument; model_size="$OPTARG" ;;
        g | gwas )       needs_argument; gwas_pattern="$OPTARG" ;;
        b | binsize )    needs_argument; bin_size="$OPTARG" ;;
        s | samplesize ) needs_argument; sample_size="$OPTARG" ;;
        n | nummodels)   needs_argument; number_of_models="$OPTARG" ;;
        \? )             usage ;;  # Illegal short options are caught by getopts
        * )              usage ;;  # Illegal long option
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

if [[ -z "${gwas_pattern}" ]]; then
    echo "No glob pattern was given for selecting gwas traits." \
    "All gwas traits will be considered."
fi

if [[ -z "${bin_size}" || -z "${sample_size}" || -z "${number_of_models}" ]]; then
    # If the user doesn't put in all of these options, our best hope is to look
    # for the first approximate match
    input_directory=$( \
    find "${OPTIMUM_STATES_DIR}" -type d \
    -name "BinSize_*${bin_size}*_SampleSize_*${sample_size}*_*${number_of_models}*" | \
    head -1)

    bin_size=$(basename "${input_directory}" | cut -d_ -f2)
    sample_size=$(basename "${input_directory}" | cut -d_ -f4)
    number_of_models=$(basename "${input_directory}" | cut -d_ -f5)
else
input_directory="${OPTIMUM_STATES_DIR}\
/BinSize_${bin_size}_SampleSize_${sample_size}_${number_of_models}"
fi

if [[ -z "${model_size}" ]]; then
    if [[ -z "$(ls -A "${input_directory}")" ]]; then
        { >&2 echo -e "ERROR: No files found in: ${input_directory}.\n"\
        "Please run 6_OptimumNumberOfStates.sh before this script."
        finishing_statement 1; }
    fi
    optimum_state_file="${input_directory}/OptimumNumberOfStates.txt"
    model_size=$(tail -1 "${optimum_state_file}" | \
    cut -d: -f2 | \
    tr -d ' ')
    echo "Optimum model size found was: ${model_size}"
fi

output_directory="${LD_ASSESSMENT_DIR}\
/BinSize_${bin_size}_SampleSize_${sample_size}_${model_size}"

rm -rf "${output_directory:?}"
mkdir -p "${output_directory:?}/annotation"
mkdir -p "${output_directory:?}/heritability"
mkdir -p "${output_directory:?}/plots"

## ============================ ##
##   ANNOTATION FILE CREATION   ##
## ============================ ##

# We need to classify all of the SNPs from the 1000 genomes data into our
# state assignments from ChromHMM to create the partitions (annotation files)

module purge
module load R/4.2.1-foss-2022a

full_model_directory="${MODEL_DIR}\
/BinSize_${bin_size}_SampleSize_${sample_size}_${number_of_models}"

dense_bed_file=$(\
find "${full_model_directory}" \
-name "*_${model_size}_dense.bed")

cd "${RSCRIPTS_DIR}" || \
{ >&2 echo "ERROR: [\${RSCRIPTS_DIR} - ${RSCRIPTS_DIR}] doesn't exist, \
make sure FilePaths.txt is pointing to the correct directory"
finishing_statement 1; }

# We ignore non-autosomal chromosomes as 1000 genomes doesn't provide this data
for chromosome in {1..22}; do
    assignment_start=$(date +%s)
    bim_file="${LD_PLINK_DIR}/${PLINK_PREFIX}.${chromosome}.bim"

    Rscript SNPAssignment.R \
    "${model_size}" \
    <(awk -v chromosome="chr${chromosome}" \
    '$1 == chromosome {print $2, $3, $4}' "${dense_bed_file}") \
    <(awk '{print $1, $2, $3, $4}' "${bim_file}") \
    "${output_directory}/annotation/ChromHMM.${chromosome}.annot"

    assignment_end=$(date +%s)
    time_taken=$((assignment_end - assignment_start))
    echo "Chromosome ${chromosome} took: ${time_taken} seconds"
done

## ======================= ##
##   REFERENCE LD SCORES   ##
## ======================= ##

module purge
module load Anaconda3/2020.02

source "${CONDA_SHELL}/profile.d/conda.sh" || \
{ echo "profile.d/conda.sh does not exist in specified location: \
[\${CONDA_SHELL} - ${CONDA_SHELL}]"; exit 1; }
conda activate "${LDSC_ENVIRONMENT}"


for chromosome in {1..22}; do
    python \
    "${LD_SOFTWARE_DIR}/ldsc.py" \
    --l2 \
    --bfile      "${LD_PLINK_DIR}/${PLINK_PREFIX}.${chromosome}" \
    --ld-wind-cm 1 \
    --annot      "${output_directory}/annotation/ChromHMM.${chromosome}.annot" \
    --out        "${output_directory}/annotation/ChromHMM.${chromosome}"
done

## ============================ ##
##   PARTITIONED HERITABILITY   ##
## ============================ ##

gwas_traits=$(find "${LD_GWAS_TRAITS_DIR}" -name "*${gwas_pattern}*.sumstats*")

for file_name in "${gwas_traits[@]}"; do
    output_file=$(basename "${file_name}" .sumstats.gz)

    # No need for overlap frq files here as state assignments are necessarily
    # distinct categories
	python \
	"${LD_SOFTWARE_DIR}/ldsc.py" \
	--h2          "${file_name}" \
	--ref-ld-chr  "${output_directory}/annotation/ChromHMM." \
	--w-ld-chr    "${WEIGHTS_DIR}/${WEIGHTS_PREFIX}." \
	--out         "${output_directory}/heritability/${output_file}"
done

## ====================== ##
##   DATA VISUALISATION   ##
## ====================== ##

conda deactivate
module purge
module module load R/4.2.1-foss-2022a

Rscript HeritabilityHeatmap.R \
<(find "${output_directory}/heritability" -name "*${gwas_pattern}*.results") \
"${output_directory}/plots"
