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

## =================================================================================##
##                                                                                  ||
##                                     PREAMBLE                                     ||
##                                                                                  ||
## =================================================================================##
## PURPOSE:                                                                         ||
## Uses ChromHMM's LearnModel command to generate models with increasing numbers of ||
## states. The idea here is to increment the number of states in the hidden Markov  ||
## model so that they can later be compared when testing for redundant states,      ||
## allowing for an 'intelligent' choice for the number of states in the model.      ||
##                                                                                  ||
## IMPORTANT NOTE: The number of states in any one model cannot exceed 2^k, where   ||
## k is the number of marks in the binary files. This is because the 'information'  ||
## method is being used by ChromHMM's LearnModel command for reproducability. This  ||
## 2^k limit is a hard cap, but depending on the data, a smaller soft cap may       ||
## exist. ChromHMM outputs the maximum number of states allowed if the cap is       ||
## exceeded, check the error logs for this message.                                 ||
## =================================================================================##
## AUTHOR: Sam Fletcher                                                             ||
## CONTACT: s.o.fletcher@exeter.ac.uk                                               ||
## CREATED: November 2023                                                           ||
## =================================================================================##
## PREREQUISITES: Run 4_BinarizeBamFiles.sh                                         ||
## =================================================================================##
## DEPENDENCIES: Java, ChromHMM                                                     ||
## =================================================================================##
## INPUTS:                                                                          ||
## $1 -> Number of models to learn (default: 4)                                     ||
## $2 -> The increment to use between model sizes (default: 1)                      ||
## $3 -> The bin size to use                                                        ||
## $4 -> The sample size to use                                                     ||
## $5 -> The assembly to use (default: hg19)                                        ||
## =================================================================================##
## OUTPUTS:                                                                         ||
## Emission parameter matrix for models (.png, .txt and .svg)                       ||
## Transition parameter matrix for models (.png, .txt and .svg)                     ||
## Full model files                                                                 ||
## The overlap and fold enrichment with existing genomic annotations for models     ||
## The maximum achieved estimated log likelihood achieved by each model             ||
## =================================================================================##

## ======================== ##
##    HELP FUNCTIONALITY    ##
## ======================== ##

if [[ "$1" == "--help" || "$1" == "-h" ]]; then
    echo "======================================================================"
    echo "Purpose: Uses ChromHMM's LearnModel command to generate several models"
    echo "with increasing numbers of states."
    echo "Author: Sam Fletcher"
    echo "Contact: s.o.fletcher@exeter.ac.uk"
    echo "Dependencies: Java, ChromHMM"
    echo "Inputs:"
    echo "\$1 -> Number of models to learn (default: 4)"
    echo "\$2 -> The increment to use between model sizes (default: 1)"
    echo "\$3 -> The bin size to use"
    echo "\$4 -> The sample size to use"
    echo "\$5 -> The assembly to use (default: hg19)"
    echo "Optional:"
    echo "Specify --array in sbatch options, to set a custom array size."
    echo "======================================================================"
    exit 0
fi

## ============ ##
##    SET UP    ##
## ============ ##

# CHANGE THESE TO YOUR OWN CONFIG FILES
source "/lustre/projects/Research_Project-MRC190311/scripts/integrative\
/ChromOptimise/configuration/FilePaths.txt"
source "/lustre/projects/Research_Project-MRC190311/scripts/integrative\
/ChromOptimise/configuration/LogFileManagement.sh"

# Output and error files renamed to:
# [job id]~[array id]~[date]-[time]

mv "${SLURM_SUBMIT_DIR}/temp${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log" \
"${LOG_FILE_PATH}/${SLURM_ARRAY_JOB_ID}~${SLURM_ARRAY_TASK_ID}~${timestamp:=}.log"
mv "${SLURM_SUBMIT_DIR}/temp${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err" \
"${LOG_FILE_PATH}/${SLURM_ARRAY_JOB_ID}~${SLURM_ARRAY_TASK_ID}~$timestamp.err"

## =============== ##
##    VARIABLES    ##
## =============== ##

number_of_models_to_generate=$1
states_increment=$2
bin_size=$3
sample_size=$4
assembly=$5

## ====== DEFAULTS ====================================================================
if ! [[ "${number_of_models_to_generate}" =~ ^[0-9]+$ ]]; then
    number_of_models_to_generate=4
    echo "Value for 'number of models to generate' is invalid."
    echo "Using the default value of: ${number_of_models_to_generate} instead."
fi

if ! [[ "${states_increment}" =~ ^[0-9]+$ ]]; then
    states_increment=1
    echo "State increment given is invalid."
    echo "Using the default value of ${states_increment} instead."
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
    echo "No assembly was given, using the default value of ${assembly} instead."
fi
# =====================================================================================

## =============================== ##
##   CLEAN UP AND ERROR CATCHING   ##
## =============================== ##

full_binary_path="${BINARY_DIR}/BinSize_${bin_size}_SampleSize_${sample_size}"

cd "${full_binary_path}" || \
{ >&2 echo -e "ERROR: Binary directory for bin/sample size is empty.\n" \
"Ensure that 4_BinarizeBamFiles.sh has been ran before this script."
finishing_statement 1; }

# Clean up from previous runs of script
cd "${MODEL_DIR}" || \
{ >&2 echo "ERROR: [\${MODEL_DIR} - ${MODEL_DIR}] doesn't exist, \
make sure FilePaths.txt is pointing to the correct directory."
finishing_statement 1; }
rm -f ./*

cd STATEBYLINE || { finishing_statement 1; }
rm -f ./*


## ========================== ##
##   PARALLELISATION SET UP   ##
## ========================== ##

number_of_models_per_array=$((number_of_models_to_generate / SLURM_ARRAY_TASK_COUNT))
remainder=$((number_of_models_to_generate % SLURM_ARRAY_TASK_COUNT))

# If the number of models to be generated isn't a multiple of the size of the array,
# the array with the smallest id will learn the left over smaller models
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

mkdir -p "${OPTIMUM_STATES_DIR}/Likelihood_Values"
cd "${OPTIMUM_STATES_DIR}/Likelihood_Values" || finishing_statement 1

# Job is to be submitted as an array so we only want to 
# remake likelihood files for one of the array tasks, not all of them.
if [[ "${SLURM_ARRAY_TASK_ID}" -eq 1 ]]; then
    rm -f "likelihood.BinSize.${bin_size}.SampleSize.${sample_size}.txt"
    touch "likelihood.BinSize.${bin_size}.SampleSize.${sample_size}.txt"
fi

sequence=$(\
seq "$starting_number_of_states" "$states_increment" "$ending_number_of_states"\
)

echo "Learning models using a bin size of ${bin_size}..." 

# Main loop
for numstates in ${sequence}; do
    echo "Learning model with: ${numstates} states..."
    # -noautoopen used so html files are not opened after model learning finshes.
    # -nobed used as genome browser files and segmentation files are not required.
    # -printstatebyline used to get the state assignment for isolation metrics
    java -mx4G \
    -jar "${CHROMHMM_MAIN_DIR}/ChromHMM.jar" LearnModel \
    -noautoopen \
    -nobed \
    -printstatebyline \
    -b "${bin_size}" \
    "${full_binary_path}" "${MODEL_DIR}" "${numstates}" "${assembly}" > \
    "ChromHMM.Output.BinSize.${bin_size}.numstates.${numstates}.txt"

    echo -n "Writing estimated log likelihood to: "
    echo "likelihood.BinSize.${bin_size}.SampleSize.${sample_size}.txt"
    echo -n "Estimated Log Likelihood for ${numstates} states: " >> \
    "likelihood.BinSize.${bin_size}.SampleSize.${sample_size}.txt"

    # grep removes the terminal logs associated with writing to files. 
    # The tail and awk locate the final estimated log likelihood.
    grep "  " "ChromHMM.Output.BinSize.${bin_size}.numstates.${numstates}.txt" | \
    tail -1 | \
    awk '{print $2}' >> \
    "likelihood.BinSize.${bin_size}.SampleSize.${sample_size}.txt" 

    grep "  " "ChromHMM.Output.BinSize.${bin_size}.numstates.${numstates}.txt" >> \
    "${LOG_FILE_PATH}/${SLURM_ARRAY_JOB_ID}~${SLURM_ARRAY_TASK_ID}~${timestamp}.log"

    rm "ChromHMM.Output.BinSize.${bin_size}.numstates.${numstates}.txt"
done

## ========================= ##
##   RENAMING OUTPUT FILES   ##
## ========================= ##

cd "${MODEL_DIR}" || \
{ >&2 echo "ERROR: [\${MODEL_DIR} - ${MODEL_DIR}] doesn't exist, \
make sure FilePaths.txt is pointing to the correct directory."
finishing_statement 1; }

# html files are not required for subsequent analysis
rm ./*.html

# Find all files in model directory, but not the state assignment files
# We specifically want the model files that have yet to be renamed
# renamed files' first word will be in captials, so this find statement will
# not find the renamed files.
files_to_rename=$(find . -maxdepth 1 -type f -name "*ions_*")

# files will be named [emissions/transitions] (file start) followed by
# information about the file (file middle) followed by the number of states
# and file extension (file end)
file_middle="_BinSize_${bin_size}_SampleSize_${sample_size}_States_"

for file in $files_to_rename; do
    # files are originally in the form 
    # [file type]_[state_number].[file extension]
    # We want to split this into the type of the file (file start) and the 
    # number of states and file extension (file end)
    file_start=$(echo "$file" | cut -d "_" -f 1)
    file_end=$(echo "$file" | cut -d "_" -f 2)

    # We use the ^ expansion of file start to ensure that the files are not 
    # renamed multiple times
    mv "$file" "${file_start^}${file_middle}${file_end}"
done

finishing_statement 0