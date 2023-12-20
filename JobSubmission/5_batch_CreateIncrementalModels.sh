#!/bin/bash
# Export all environment variables to the batch job
#SBATCH --export=ALL
# Submit to the mrc queue for faster queue times
#SBATCH -p mrcq
# Tests: 1MB of binary data, 3 states: 55s, 4 states: 105s, 5 states: 575s,
# 6 states: 763s, 7 states: 903s, 8 states: 1500s
# Forward backwards algorithm has time complexity of N^2T, where T is the number of 
# Genomic bins and N is the number of states
#SBATCH --time=01:00:00 
#SBATCH -A Research_Project-MRC190311 
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=16
# Make sure array is not higher than the number of models being learned
# as this ends up with all models being processed by max index array
#SBATCH --array=1-4
# Tests: 1MB of binary data. 2 states: 181 MB 3 states: 145 MB, 4 states: 165 MB,
# 5 states: 128 MB, 6 states 165 MB, 7 states 168 MB 
# Memory consumption doesn't seem to be dependent on number of states
#SBATCH --mem=5G
# Send an email after the job is done
#SBATCH --mail-type=END 
# Temporary log file, later to be removed
#SBATCH --output=temp%A_%a.log
# Temporary error file, later to be removed
#SBATCH --error=temp%A_%a.err
#SBATCH --job-name=Model_Learning

## =================================================================================##
##                                                                                  ||
##                                     PREAMBLE                                     ||
##                                                                                  ||
## =================================================================================##
## PURPOSE:                                                                         ||
## Uses ChromHMM's LearnModel command to generate models with increasing numbers of ||
## states. The idea here is to increment the number of states in the hidden Markov  ||
## model so that they  can later be compared to test for redundant states, allowing ||
## for an 'intelligent' choice for the number of states in the model.               ||
##                                                                                  ||
## IMPORTANT NOTE: The number of states in any one model cannot exceed 2^k, where   ||
## k is the number of marks in the binary files. This is because the 'init' method  ||
## is being used by ChromHMM's LearnModel command for reproducability. This 2^k     ||
## limit is a hard cap, but depending on the data, a smaller soft cap may exist.    ||
## ChromHMM outputs the maximum number of states allowed if the cap is exceeded,    ||
## check the error logs for this message.                                           ||
## =================================================================================##
## AUTHOR: Sam Fletcher s.o.fletcher@exeter.ac.uk                                   ||
## CREATED: November 2023                                                           ||
## =================================================================================##
## PREREQUISITES: All .bam files for a specific epigenetic mark must be in 1 folder ||
## =================================================================================##
## DEPENDENCIES: Java, ChromHMM                                                     ||
## =================================================================================##
## INPUTS:                                                                          ||
## $1 -> Number of models to generate                                               ||
## $2  -> The increment to use between model sizes                                  ||
## $3 -> Bin size, WARNING: Use the same bin size as was used in                    ||
##       4_BinarizeBamFiles.sh                                                      ||
## $4 -> Sample Size, WARNING: Use the same sample size as was used in              ||
##       3_SubsampleBamFiles.sh                                                     ||
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
    echo "with increasing numbers of states"
    echo "Author: Sam Fletcher"
    echo "Contact: s.o.fletcher@exeter.ac.uk"
    echo "Dependencies: Java, ChromHMM"
    echo "Inputs:"
    echo "\$1 -> Number of models to generate"
    echo "\$2 -> The increment to use between model sizes"
    echo "\$3 -> Bin size, WARNING: Use the same bin size as was used in" 
    echo "4_BinarizeBamFiles.sh"
    echo "\$4 -> Sample Size, WARNING: Use the same sample size as was used in"
    echo "3_SubsampleBamFiles.sh"
    echo "Optional:"
    echo "Specify --array in sbatch options, to set a custom array size."
    echo "======================================================================"
    exit 0
fi

## ============ ##
##    SET UP    ##
## ============ ##

echo "Job '$SLURM_JOB_NAME' started at:"
date -u

start_time=$(date +%s)

# Activate config.txt to access all file paths
# CHANGE THIS TO YOUR OWN CONFIG FILE
echo "Loading config file..."
source "/lustre/projects/Research_Project-MRC190311\
/scripts/integrative/blueprint/config/config.txt"

# Rename the output and error files to have format:
# [job id]~[array id]~[date]-[time]
# This requires a hard link as you cannot rename log files
# whilst running the script without a wrapper function
LOG_FILE_PATH="${LOG_DIR}/$SLURM_JOB_NAME/$USER"
mkdir -p "${LOG_FILE_PATH}"
timestamp=$(date -u +%Y.%m.%d-%H:%M)

ln "${SLURM_SUBMIT_DIR}/temp${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log" \
"${LOG_FILE_PATH}/${SLURM_ARRAY_JOB_ID}~${SLURM_ARRAY_TASK_ID}~$timestamp.log"
ln "${SLURM_SUBMIT_DIR}/temp${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err" \
"${LOG_FILE_PATH}/${SLURM_ARRAY_JOB_ID}~${SLURM_ARRAY_TASK_ID}~$timestamp.err"


cd "${BINARY_DIR}" || { echo "Binary directory doesn't exist, \
make sure config.txt is pointing to the correct directory"; exit 1; }
if [ -z "$(ls -A)" ]; then
    echo "4_BinarizedFiles is empty."
    echo "Ensure that 4_BinarizeBamFiles.sh has been ran before this script."
    echo "Aborting..."

    rm "${SLURM_SUBMIT_DIR}/temp${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log"
    rm "${SLURM_SUBMIT_DIR}/temp${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err"

    exit 1
fi

## ========================= ##
##    VARIABLE ASSIGNMENT    ##
## ========================= ##

NUMBER_OF_MODELS_TO_GENERATE=$1
STATE_INCREMENT=$2
BIN_SIZE=$3
SAMPLE_SIZE=$4

if [ -z "${NUMBER_OF_MODELS_TO_GENERATE}" ]; then
    echo "Number of models to generate was not given."
    echo "Using the default value of 4 instead."
    NUMBER_OF_MODELS_TO_GENERATE=4
fi

if [ -z "${STATE_INCREMENT}" ]; then
    echo "The value for the state increment was not given."
    echo "Using the default value of 1 instead."
    STATE_INCREMENT=1
fi

# 'intelligently' set bin size default by searching through the subsample directory
if [ -z "${BIN_SIZE}" ]; then
    echo "The value for the bin size was not given." 
    echo -n "Assuming that the first file in the subsampled directory "
    echo "uses the correct bin size..."

    cd "${BINARY_DIR}" || { echo "Binary directory doesn't exist, \
    make sure config.txt is pointing to the correct directory"; exit 1; }

    BIN_SIZE=$(find . -type f -name "*.txt*.gz" | head -1 | cut -d "_" -f 6)
fi
# 'intelligently' set sample size default by searching through the subsample directory
if [ -z "${SAMPLE_SIZE}" ]; then
    echo "No sample size was given."
    echo -n "Assuming that the first file in the subsampled directory "
    echo "uses the correct sample size..."

    cd "${BINARY_DIR}" || { echo "Binary directory doesn't exist, \
    make sure config.txt is pointing to the correct directory"; exit 1; }

    SAMPLE_SIZE=$(find . -type f -name "*.txt*.gz" | head -1 | cut -d "_" -f 4)
fi

echo -n "Generating ${NUMBER_OF_MODELS_TO_GENERATE} models for the binary files "
echo "found in the 4_BinarizedFiles directory."
echo -n "Models will have states starting at 2 and going up in increments of: "
echo "${STATE_INCREMENT}."
echo "ChromHMM's LearnModel command will use the option '-b ${BIN_SIZE}'."

## ===================== ##
##   ADDITIONAL SET UP   ##
## ===================== ##

cd "${MODEL_DIR}" || { echo "Model directory doesn't exist, \
make sure config.txt is pointing to the correct directory"; exit 1; }
rm -f ./*

cd "${OPTIMUM_STATES_DIR}" || { echo "Optimum states directory doesn't exist, \
make sure config.txt is pointing to the correct directory"; exit 1; }
mkdir -p "Likelihood_Values"
cd "Likelihood_Values" || exit 1

## ========================== ##
##   PARALLELISATION SET UP   ##
## ========================== ##

number_of_models_per_array=$((NUMBER_OF_MODELS_TO_GENERATE / SLURM_ARRAY_TASK_COUNT))
remainder=$((NUMBER_OF_MODELS_TO_GENERATE % SLURM_ARRAY_TASK_COUNT))

# If the number of models to be generated isn't a multiple of the size of the array,
# the array with the smallest id will learn the left over smaller models
# (spreading the larger models evenly across the other array elements).
if [ "${SLURM_ARRAY_TASK_ID}" -eq 1 ]; then
    starting_number_of_states=2
    ending_number_of_states=$(( \
    2 + STATE_INCREMENT*(remainder+number_of_models_per_array-1) \
    ))
else
    starting_number_of_states=$(( \
    (((SLURM_ARRAY_TASK_ID-1)*number_of_models_per_array) + remainder)*\
    STATE_INCREMENT + 2 \
    )) 
    ending_number_of_states=$(( \
    (((SLURM_ARRAY_TASK_ID)*number_of_models_per_array) + remainder -1 )*\
    STATE_INCREMENT + 2 \
    ))
fi

## ===================== ##
##    MODEL GENERATION   ##
## ===================== ##

module purge
module load Java

# Job is to be submitted as an array so we only want to 
# remake likelihood files for one of the array tasks, not all of them.
if [[ "${SLURM_ARRAY_TASK_ID}" -eq 1 ]]; then
    rm -f "likelihood.BinSize.${BIN_SIZE}.SampleSize.${SAMPLE_SIZE}.txt"
    touch "likelihood.BinSize.${BIN_SIZE}.SampleSize.${SAMPLE_SIZE}.txt"
fi

sequence=$(\
seq "$starting_number_of_states" "$STATE_INCREMENT" "$ending_number_of_states"\
)

for numstates in ${sequence}; do
    echo "Learning model with: ${numstates} states..."
    # The -nobed option is here as we have no need for the genome browser files
    # or segmentation files.
    java -mx4G \
    -jar "${CHROMHMM_MAIN_DIR}/ChromHMM.jar" LearnModel \
    -noautoopen \
    -nobed \
    -b "${BIN_SIZE}" \
    "${BINARY_DIR}" "${MODEL_DIR}" "${numstates}" hg19 > \
    "ChromHMM.Output.BinSize.${BIN_SIZE}.numstates.${numstates}.txt"

    echo -n "Writing estimated log likelihood to: "
    echo "likelihood.BinSize.${BIN_SIZE}.SampleSize.${SAMPLE_SIZE}.txt"
    echo -n "Estimated Log Likelihood for ${numstates} states: " >> \
    "likelihood.BinSize.${BIN_SIZE}.SampleSize.${SAMPLE_SIZE}.txt"

    # grep removes the terminal logs associated with writing to files. 
    # The tail and awk locate the final estimated log likelihood.
    grep "       " "ChromHMM.Output.BinSize.${BIN_SIZE}.numstates.${numstates}.txt" | \
    tail -1 | \
    awk '{print $2}' >> \
    "likelihood.BinSize.${BIN_SIZE}.SampleSize.${SAMPLE_SIZE}.txt" 

    rm "ChromHMM.Output.BinSize.${BIN_SIZE}.numstates.${numstates}.txt"
done

## ========================= ##
##   RENAMING OUTPUT FILES   ##
## ========================= ##

cd "${MODEL_DIR}" || { echo "Model directory doesn't exist, \
make sure config.txt is pointing to the correct directory"; exit 1; }
emission_files_to_rename=$(find . -type f -name "emissions*")
for file in $emission_files_to_rename; do
    file_ending=$(echo "$file" | cut -d "_" -f 2) 
    mv "$file" \
    "Emissions_BinSize_${BIN_SIZE}_SampleSize_${SAMPLE_SIZE}_NumberOfStates_${file_ending}" 
done

transistion_files_to_rename=$(find . -type f -name "transitions*")
for file in $transistion_files_to_rename; do
    file_ending=$(echo "$file" | cut -d "_" -f 2)
    mv "$file" \
    "Transitions_BinSize_${BIN_SIZE}_SampleSize_${SAMPLE_SIZE}_NumberOfStates_${file_ending}" 
done

model_files_to_rename=$(find . -type f -name "model*")
for file in $model_files_to_rename; do
    file_ending=$(echo "$file" | cut -d "_" -f 2)
    mv "$file" \
    "Model_BinSize_${BIN_SIZE}_SampleSize_${SAMPLE_SIZE}_NumberOfStates_${file_ending}" 
done


## ======================= ##
##   LOG FILE MANAGEMENT   ##
## ======================= ##

echo "Job completed at:"
date -u
end_time=$(date +%s)
time_taken=$((end_time-start_time))
echo "Job took a total of: ${time_taken} seconds to complete"

rm "${SLURM_SUBMIT_DIR}/temp${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log"
rm "${SLURM_SUBMIT_DIR}/temp${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err"