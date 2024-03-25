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
## Annotation files for each chromosome (annotations being all baseline       ||
## categories and state assignments from ChromHMM).                           ||
## LD scores for each category for each chromosome.                           ||
## ===========================================================================##


## ===================== ##
##   ARGUMENT PARSING    ##
## ===================== ##

usage() {
cat <<EOF
===========================================================================
7_ReferenceLDSCore
===========================================================================
Purpose: Generates reference LD scores after generating annotation files
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
        c | config )      needs_argument; configuration_directory="$OPTARG" ;;
        t | state )       needs_argument; model_size="$OPTARG" ;;
        g | gwas )        needs_argument; gwas_pattern="$OPTARG" ;;
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
mv "${SLURM_SUBMIT_DIR}/temp${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log" \
"${LOG_FILE_PATH}/${SLURM_ARRAY_JOB_ID}~${SLURM_ARRAY_TASK_ID}~${timestamp:=}.log"
mv "${SLURM_SUBMIT_DIR}/temp${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err" \
"${LOG_FILE_PATH}/${SLURM_ARRAY_JOB_ID}~${SLURM_ARRAY_TASK_ID}~$timestamp.err"


## =============== ##
##    VARIABLES    ##
## =============== ##

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
/BinSize_${bin_size}_SampleSize_${sample_size}_${number_of_models}"

# If jobs get queued, some files will be deleted prematurely
if [[ "${SLURM_ARRAY_TASK_ID}" -eq 1 ]]; then
    rm -rf "${output_directory:?}"
    mkdir -p "${output_directory}/annotation" \
    "${output_directory}/heritability" \
    "${output_directory}/plots/State_Categories" \
    "${output_directory}/plots/All_Categories" 
fi

# We sleep here to ensure files are not removed prematurely
sleep 5

## =================== ##
##   FILE MANAGEMENT   ##
## =================== ##

full_model_directory="${MODEL_DIR}\
/BinSize_${bin_size}_SampleSize_${sample_size}_${number_of_models}"

full_binary_directory="${BINARY_DIR}\
/BinSize_${bin_size}_SampleSize_${sample_size}"

# We ignore non-autosomal chromosomes as 1000 genomes doesn't provide this data
# Hence our chromosomes are just 1-22 (the array indices)
chromosome=${SLURM_ARRAY_TASK_ID}

# The production of the annotation file requires some intermediary files
# so we create a temporary directory to hold these
temporary_directory="${output_directory}/temp_${chromosome}"
mkdir -p "${temporary_directory}"

dense_bed_file=$(\
find "${full_model_directory}" \
-name "*_${model_size}_dense.bed")

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

module purge
module load R/4.2.1-foss-2022a

Rscript BinarytoBed.R \
<(zcat "${binary_file}") \
"${bin_size}" \
"chr${chromosome}" \
"${temporary_directory}/binary-${chromosome}.bed"

Rscript BimtoBed.R \
<(cat "${bim_file}") \
"${temporary_directory}/SNP_positions-${chromosome}.bed"

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
awk 'NR==2' > \
"${temporary_directory}/mark_assignments-${chromosome}.txt"

bedtools intersect -wb \
-a "${temporary_directory}/SNP_positions-${chromosome}.bed" \
-b "${temporary_directory}/binary-${chromosome}.bed" | \
awk '{ for (i=7; i<=NF; i++) printf "%s%s", $i, (i<NF ? "\t" : "\n") }' >> \
"${temporary_directory}/mark_assignments-${chromosome}_temp.txt"

module purge
module load R/4.2.1-foss-2022a

baseline_annot="${LD_BASELINE_DIR}/baselineLD.${chromosome}.annot.gz"

Rscript CreateAnnotationFile.R \
<(zcat "${baseline_annot}") \
<(cat "${temporary_directory}/state_assignments-${chromosome}.txt") \
<(cat "${temporary_directory}/mark_assignments-${chromosome}.txt") \
"${model_size}" \
"${output_directory}/annotation/ChromHMM.${chromosome}.annot"

rm -rf "${temporary_directory}"

## ======================= ##
##   REFERENCE LD SCORES   ##
## ======================= ##

module purge
module load Anaconda3/2020.02

source "${CONDA_SHELL}/profile.d/conda.sh" || \
{ echo "profile.d/conda.sh does not exist in specified location: \
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
--annot      "${output_directory}/annotation/ChromHMM.${chromosome}.annot" \
--out        "${output_directory}/annotation/ChromHMM.${chromosome}"


## ============================ ##
##   PARTITIONED HERITABILITY   ##
## ============================ ##

# We only need to calculate partitioned heritability (and produce the plots)
# once, so if the number of chromosomes completed isn't yet 22, we exit early.
chromosomes_completed=$(\
find "${output_directory}/annotation/" \
-name "ChromHMM*.log" -size 0+ | \
wc -l)

if [[ "${chromosomes_completed}" -ne 22 ]]; then
    echo "Not all ld scores have been calculated in this array."
    echo "Exiting early..."
    finishing_statement 0
fi

# This job is being ran as an array, which means that the memory of the job
# is split among each array task. This memory allocation is not dynamic and so
# at this point in the program the task only has (max memory)/22 GB of memory.
# This is not enough to handle the partitioned heritability (unless you have
# ~220GB of memory available). Hence at this point we run a new script with
# sbatch

cd "${SCRIPTS_DIR}/JobSubmission" || { echo "Could not find the JobSubmission \
directory in ${SCRIPTS_DIR}/Jobsubmission. Please check your configuration \
file."; finishing_statement 0; }

sbatch 8_PartitionedHeritability.sh \
--config="${configuration_directory}" \
--gwas="${gwas_pattern}" \
--binsize="${bin_size}" \
--samplesize="${sample_size}" \
--nummodels="${number_of_models}"

finishing_statement 0
