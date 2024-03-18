#!/bin/bash
# Export all environment variables to the batch job
#SBATCH --export=ALL
# Submit to the mrc queue for faster queue times
#SBATCH -p mrcq
# Consult information/Processing_Times.md for expected time
# Forward backwards algorithm has time complexity of O(N^2T), where T is the number of
# genomic bins and N is the number of states
#SBATCH --time=01:00:00 
#SBATCH -A Research_Project-MRC190311 
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=16
# Make sure array is not higher than the number of models being learned
# as this ends up with all models being processed by max index array
#SBATCH --array=1-4
# Memory consumption doesn't appear to be dependent on number of states
# Consult information/Memory_Profiling.md for expected memory usage
#SBATCH --mem=5G
# Send an email after the job is done
#SBATCH --mail-type=END 
# Temporary log file, later to be removed
#SBATCH --output=temp%A_%a.log
# Temporary error file, later to be removed
#SBATCH --error=temp%A_%a.err
#SBATCH --job-name=5_Model_Learning

## ===========================================================================##
##                                                                            ||
##                                  PREAMBLE                                  ||
##                                                                            ||
## ===========================================================================##
## PURPOSE:                                                                   ||
## Uses ChromHMM's LearnModel command to generate models with increasing      ||
## numbers of states. The idea here is to increment the number of states in   ||
## the HMM so that they can later be compared when testing for redundant      ||
## states, allowing for an 'intelligent' choice for the number of states in   ||
## the model.                                                                 ||
##                                                                            ||
## IMPORTANT NOTE: The number of states in any one model cannot exceed 2^k,   ||
## where k is the number of marks in the binary files. This is because the    ||
## 'information' method is being used by ChromHMM's LearnModel command for    ||
## reproducability. This 2^k limit is a hard cap, but depending on the data,  ||
## a smaller soft cap may exist.                                              ||
## ChromHMM outputs the maximum number of states allowed if the cap is        ||
## exceeded, check the error logs for this message.                           ||
## ===========================================================================##
## AUTHOR: Sam Fletcher                                                       ||
## CONTACT: s.o.fletcher@exeter.ac.uk                                         ||
## CREATED: November 2023                                                     ||
## ===========================================================================##
## PREREQUISITES: Run 4_BinarizeBamFiles.sh                                   ||
## ===========================================================================##
## DEPENDENCIES: Java, ChromHMM                                               ||
## ===========================================================================##
## INPUTS:                                                                    ||
## -c|--config=     -> Full/relative file path for configuation file directory||
## -n|--nummodels=  -> Number of models to learn (default: 4)                 ||
## -b|--binsize=    -> The bin size used in 4_BinarizeBamFiles                ||
## -s|--samplesize= -> The sample size used in 3_SubsampleBamFiles            ||
## -a|--assembly=   -> The assembly to use (default: hg19)                    ||
## ===========================================================================##
## OUTPUTS:                                                                   ||
## Emission parameter matrix for models (.png, .txt and .svg)                 ||
## Transition parameter matrix for models (.png, .txt and .svg)               ||
## Full model files                                                           ||
## The overlap and fold enrichment with existing genomic annotations for      ||
##   models                                                                   ||
## The maximum achieved estimated log likelihood achieved by each model       ||
## ===========================================================================##

## ===================== ##
##   ARGUMENT PARSING    ##
## ===================== ##

usage() {
cat <<EOF
===========================================================================
5_batch_CreateIncrementalModels
===========================================================================
Purpose: Uses ChromHMM's LearnModel command to generate several models
with increasing numbers of states.
Author: Sam Fletcher
Contact: s.o.fletcher@exeter.ac.uk
Dependencies: Java, ChromHMM
Inputs:
-c|--config=     -> Full/relative file path for configuation file directory
-n|--nummodels=  -> Number of models to learn (default: 4)
-b|--binsize=    -> The bin size used in 4_BinarizeBamFiles
-s|--samplesize= -> The sample size used in 3_SubsampleBamFiles
-a|--assembly=   -> The assembly to use (default: hg19)
===========================================================================
EOF
    exit 0
}

needs_argument() {
    # Required check in case user uses -a -b or -b -a (no argument given).
    if [[ -z "$OPTARG" || "${OPTARG:0:1}" == - ]]; then usage; fi
}

if [[ ! $1 =~ -.* ]]; then usage; fi

while getopts c:n:i:b:s:a:-: OPT; do
    # Adds support for long options by reformulating OPT and OPTARG
    # This assumes that long options are in the form: "--long=option"
    if [ "$OPT" = "-" ]; then
        OPT="${OPTARG%%=*}"
        OPTARG="${OPTARG#"$OPT"}"
        OPTARG="${OPTARG#=}"
    fi
    case "$OPT" in
        c | config )       needs_argument; configuration_directory="$OPTARG" ;;
        n | nummodels)     needs_argument; number_of_models_to_generate="$OPTARG" ;;
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
mv "${SLURM_SUBMIT_DIR}/temp${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log" \
"${LOG_FILE_PATH}/${SLURM_ARRAY_JOB_ID}~${SLURM_ARRAY_TASK_ID}~${timestamp:=}.log"
mv "${SLURM_SUBMIT_DIR}/temp${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err" \
"${LOG_FILE_PATH}/${SLURM_ARRAY_JOB_ID}~${SLURM_ARRAY_TASK_ID}~$timestamp.err"

## =============== ##
##    VARIABLES    ##
## =============== ##

## ====== DEFAULTS =============================================================
if ! [[ "${number_of_models_to_generate}" =~ ^[0-9]+$ ]]; then
    number_of_models_to_generate=4
    echo "Value for 'number of models to generate' is invalid."
    echo "Using the default value of: ${number_of_models_to_generate} instead."
fi

if ! [[ "${bin_size}" =~ ^[0-9]+$  || "${sample_size}" =~ ^[0-9]+$ ]]; then
    cd "${BINARY_DIR}" || \
    { >&2 echo "ERROR: [\${BINARY_DIR} - ${BINARY_DIR}] doesn't exist, \
    make sure FilePaths.txt is pointing to  the correct directory."
    finishing_statement 1; }

    # The bin size and sample size are given in the file names of the binary
    # files. The following two commands find the first binary file in the
    # binary directory and extract the bin and sample size. Which, given
    # the directory structure, will be the 9th and 7th fields respectively
    # of the file path when delimited by underscores.
    bin_size=$(find . -type f -name "*.txt*.gz" | head -1 | cut -d "_" -f 9)
    sample_size=$(find . -type f -name "*.txt*.gz" | head -1 | cut -d "_" -f 7)
    
    echo "Bin size or sample size given is invalid."
    echo "Using the following values instead."
    echo "Bin size: ${bin_size}. Sample Size: ${sample_size}."
fi


if [[ -z "${assembly}" ]]; then
    assembly=hg19
    echo "No assembly was given, using the default value of" \
    "${assembly} instead."
fi
# ==============================================================================

## =============================== ##
##   CLEAN UP AND ERROR CATCHING   ##
## =============================== ##

full_binary_path="${BINARY_DIR}/BinSize_${bin_size}_SampleSize_${sample_size}"

cd "${full_binary_path}" || \
{ >&2 echo -e "ERROR: Binary directory for bin/sample size is empty.\n" \
"Ensure that 4_BinarizeBamFiles.sh has been ran before this script."
finishing_statement 1; }

# Clean up from previous runs of script
output_directory="${MODEL_DIR:?}\
/BinSize_${bin_size}_SampleSize_${sample_size}_${number_of_models_to_generate}"

sleep 5
# Avoid stale file handles
if [[ "${SLURM_ARRAY_TASK_ID}" -eq 1 ]]; then
    rm -rf "${output_directory}"
    mkdir -p "${output_directory}/Likelihood_Values"
    mkdir -p "${output_directory}/STATEBYLINE"
fi

cd "${output_directory}/Likelihood_Values" || finishing_statement 1


## ========================== ##
##   PARALLELISATION SET UP   ##
## ========================== ##

number_of_models_per_array=$((number_of_models_to_generate / SLURM_ARRAY_TASK_COUNT))
remainder=$((number_of_models_to_generate % SLURM_ARRAY_TASK_COUNT))

# By default, we assume that the user wants to learn every model from 2 upwards
states_increment=1

# If the number of models to be generated isn't a multiple of the size of the
# array, the array with the smallest id will learn the left over smaller models
# (spreading the larger models evenly across the other array elements).
if [[ "${SLURM_ARRAY_TASK_ID}" -eq 1 ]]; then
    starting_number_of_states=2
    ending_number_of_states=$(( \
    2 + states_increment*(remainder+number_of_models_per_array-1) \
    ))
else
    starting_number_of_states=$(( \
    (((SLURM_ARRAY_TASK_ID-1)*number_of_models_per_array) + remainder)*\
    states_increment + 2 \
    )) 
    ending_number_of_states=$(( \
    (((SLURM_ARRAY_TASK_ID)*number_of_models_per_array) + remainder -1 )*\
    states_increment + 2 \
    ))
fi

## ===================== ##
##    MODEL GENERATION   ##
## ===================== ##

module purge
module load Java

sequence=$(\
seq "$starting_number_of_states" "$states_increment" "$ending_number_of_states"\
)

echo "Learning models using a bin size of ${bin_size}..." 

# Main loop
for numstates in ${sequence}; do
    echo "Learning model with: ${numstates} states..."

    # -noautoopen used so html files are not opened after model learning finshes.
    # -printstatebyline used to get the state assignment for isolation metrics
    java -mx4G \
    -jar "${CHROMHMM_MAIN_DIR}/ChromHMM.jar" LearnModel \
    -noautoopen \
    -printstatebyline \
    -b "${bin_size}" \
    "${full_binary_path}" "${output_directory}" "${numstates}" "${assembly}" > \
    "ChromHMM.output.numstates.${numstates}.txt"

    echo "Writing estimated log likelihood to: likelihoods.txt"
    echo "Estimated Log Likelihood for ${numstates} states: " >> \
    "likelihoods.txt"

    # grep selects terminal logs that are not associated with writing to files.
    grep "  " "ChromHMM.output.numstates.${numstates}.txt" | \
    tail -1 | \
    awk '{print $2}' >> \
    "likelihoods.txt"

    # Instead of storing chromHMM's log file in a separate location, it is
    # easier to just store the log in the existing log file
    grep "  " "ChromHMM.output.numstates.${numstates}.txt" >> \
    "${LOG_FILE_PATH}/${SLURM_ARRAY_JOB_ID}~${SLURM_ARRAY_TASK_ID}~${timestamp}.log"

    rm "ChromHMM.output.numstates.${numstates}.txt"
done

## ========================= ##
##   RENAMING OUTPUT FILES   ##
## ========================= ##

# ChromHMM's output file names are not particularly descriptive for this
# pipeline. A key feature of this pipeline is the bin and sample size used
# The below code is in place to rename our output files so that this information
# can be inferred from the file names.

cd "${output_directory}" || finishing_statement 1

# html files are not required for subsequent analysis
rm ./*.html

finishing_statement 0
