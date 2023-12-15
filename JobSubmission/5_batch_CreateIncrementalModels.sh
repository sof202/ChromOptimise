#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -p mrcq # submit to the mrc queue for faster queue times
#SBATCH --time=01:00:00 # This can be a lengthy process. In a mock run that used information from 6 bam files, the end took 15 minutes
#SBATCH -A Research_Project-MRC190311 
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=16 
#SBATCH --array=1-4 # Change the number of arrays as you see fit. Default is 4. Ensure that the number of models being learned is greater than or equal to the number of arrays
#SBATCH --mem=10G # specify bytes memory to reserve
#SBATCH --mail-type=END # Send an email after the job is done
#SBATCH --output=/lustre/projects/Research_Project-MRC190311/scripts/integrative/blueprint/LogFiles/Model_Learning/temp%a.o #Put output file in log files with with temporary name
#SBATCH --error=/lustre/projects/Research_Project-MRC190311/scripts/integrative/blueprint/LogFiles/Model_Learning/temp%a.e #Put error file in log files with with temporary name
#SBATCH --job-name=Model_Learning

## -------------------------------------------------------------------------------------------- ##
##                                                                                              ##
##                                            PREAMBLE                                          ##
##                                                                                              ##
## -------------------------------------------------------------------------------------------- ##
##                                            PURPOSE                                           ##
##         Uses ChromHMM's LearnModel command to generate several models with increasing        ##
##           numbers of states. The idea here is to increment the number of states in           ##
##   the hidden Markov model so that they  can later be compared to test for redundant states,  ##
##        allowing  for an 'intelligent' choice for the number of states in the model.          ##
##                                                                                              ##
##                IMPORTANT NOTE: The number of states in any one model cannot                  ##
##       exceed 2^k, where k is the number of marks in the binary files. This is because        ##
##     the 'init' method is being used by ChromHMM's LearnModel command for reproducability.    ##
##    This 2^k limit is a hard cap, but depending on the data, a smaller soft cap may exist.    ##
##     ChromHMM outputs the maximum number of states allowed if the cap is exceeded, check      ##
##                              the error logs for this message.                                ##
## -------------------------------------------------------------------------------------------- ##
##                        AUTHOR: Sam Fletcher s.o.fletcher@exeter.ac.uk                        ##
##                                     CREATED: November 2023                                   ##
## -------------------------------------------------------------------------------------------- ##
##                                         PREREQUISITES                                        ##
##                                   Run: 4_BinarizeBamFiles.sh                                 ##
## -------------------------------------------------------------------------------------------- ##
##                                          DEPENDENCIES                                        ##
##                                              Java                                            ##
##                                            ChromHMM                                          ##
## -------------------------------------------------------------------------------------------- ##
##                                             INPUTS                                           ##
##                               $1 -> Number of models to generate                             ##
##                        $2  -> The increment to use between model sizes                       ##
##      $3 -> Bin size, WARNING: Use the same bin size as was used in 4_BinarizeBamFiles.sh     ##
##  $4 -> Sample Size, WARNING: Use the same sample size as was used in 3_SubsampleBamFiles.sh  ##
## -------------------------------------------------------------------------------------------- ##
##                                            OUTPUTS                                           ##
##        The emission parameter matrix of the models produced in .png, .txt and .svg form      ##
##       The transition parameter matrix of the models produced in .png, .txt and .svg form     ##
##               The full model files, including the two parameter matrices above               ##
##                     The overlap and fold enrichment files for each model                     ##
##            The maximum achieved estimated log likelihood achieved by each model              ##
## -------------------------------------------------------------------------------------------- ##

## ------------------------ ##
##    HELP FUNCTIONALITY    ##
## ------------------------ ##

if [[ "$1" == "--help" || "$1" == "-h" ]]; then
    echo "========================================================================================================"
    echo "Purpose: Uses ChromHMM's LearnModel command to generate several models with increasing numbers of states"
    echo "Author: Sam Fletcher"
    echo "Contact: s.o.fletcher@exeter.ac.uk"
    echo "Dependencies: Java, ChromHMM"
    echo "Inputs:"
    echo "\$1 -> Number of models to generate"
    echo "\$2 -> The increment to use between model sizes"
    echo "\$3 -> Bin size, WARNING: Use the same bin size as was used in 4_BinarizeBamFiles.sh"
    echo "\$4 -> Sample Size, WARNING: Use the same sample size as was used in 3_SubsampleBamFiles.sh"
    echo "========================================================================================================"
    exit 0
fi

## ------------ ##
##    SET UP    ##
## ------------ ##

# Rename the output and error files to have format: [Bin Size]~[job id]~[array id]~[date]-[time]
# This requires a hard link as you cannot rename log files whilst running the script without a wrapper function
LOG_FILE_PATH=/lustre/projects/Research_Project-MRC190311/scripts/integrative/blueprint/LogFiles/$SLURM_JOB_NAME/
mkdir -p "${LOG_FILE_PATH}"
cd "${LOG_FILE_PATH}" || exit 1
timestamp=$(date -u +%Y.%m.%d-%H:%M)
ln "temp${SLURM_ARRAY_TASK_ID}.e" "$3~${SLURM_ARRAY_JOB_ID}~${SLURM_ARRAY_TASK_ID}~$timestamp.e"
ln "temp${SLURM_ARRAY_TASK_ID}.o" "$3~${SLURM_ARRAY_JOB_ID}~${SLURM_ARRAY_TASK_ID}~$timestamp.o"

# Print start date/time
echo "Job '$SLURM_JOB_NAME' started at:"
date -u

# Get the start time for the program
start_time=$(date +%s)

# Activate config.txt to access all file paths
echo "Loading config file..."
source /lustre/projects/Research_Project-MRC190311/scripts/integrative/blueprint/config/config.txt

# Set up variables
NUMBER_OF_MODELS_TO_GENERATE=$1
STATE_INCREMENT=$2
BIN_SIZE=$3
SAMPLE_SIZE=$4

# Check for existence of files in the binarized directory
cd "${BINARY_DIR}" || { echo "Binary directory doesn't exist, make sure config.txt is pointing to the correct directory"; exit 1; }
if [ -z "$(ls -A)" ]; then
    echo "4_BinarizedFiles is empty, ensure that 4_BinarizeBamFiles.sh has been ran before this script."
    echo "Aborting..."

    # Remove temporary log files
    cd "${LOG_FILE_PATH}" || exit 1
    rm "temp${SLURM_ARRAY_TASK_ID}.e"
    rm "temp${SLURM_ARRAY_TASK_ID}.o"

    exit 1
fi

## SET DEFAULTS FOR VARIABLES ##
## ============================================================================================================================================ ##
if [ -z "${NUMBER_OF_MODELS_TO_GENERATE}" ]; then
    echo "Number of models to generate was not given. Using the default value of 4 instead"
    NUMBER_OF_MODELS_TO_GENERATE=4
fi

if [ -z "${STATE_INCREMENT}" ]; then
    echo "The value for the state increment was not given. Using the default value of 1 instead"
    STATE_INCREMENT=1
fi

# Obtains a bin size default 'intelligently' by searching through the subsample directory
if [ -z "${BIN_SIZE}" ]; then
    echo "The value for the bin size was not given.  Assuming that the first file in the subsampled directory uses the correct bin size..."
    cd "${BINARY_DIR}" || { echo "Binary directory doesn't exist, make sure config.txt is pointing to the correct directory"; exit 1; }
    BIN_SIZE=$(find . -type f -name "*.txt*.gz" | head -1 | cut -d "_" -f 6)
fi

# Obtains a sample size default 'intelligently' by searching through the subsample directory
if [ -z "${SAMPLE_SIZE}" ]; then
    echo "No sample size was given. Assuming that the first file in the subsampled directory uses the correct sample size..."
    cd "${BINARY_DIR}" || { echo "Binary directory doesn't exist, make sure config.txt is pointing to the correct directory"; exit 1; }
    SAMPLE_SIZE=$(find . -type f -name "*.txt*.gz" | head -1 | cut -d "_" -f 4)
fi
## ============================================================================================================================================ ##


echo -n "Generating ${NUMBER_OF_MODELS_TO_GENERATE} models for the binary files found in the 4_BinarizedFiles directory."
echo "Models will have states starting at 2 and going up in increments of: ${STATE_INCREMENT}. ChromHMM's LearnModel command will use the option '-b ${BIN_SIZE}'."

# Delete any model files that are currently in the 5_ModelFiles directory
cd "${MODEL_DIR}" || { echo "Model directory doesn't exist, make sure config.txt is pointing to the correct directory"; exit 1; }
rm -f ./*

cd "${OPTIMUM_STATES_DIR}" || { echo "Optimum states directory doesn't exist, make sure config.txt is pointing to the correct directory"; exit 1; }
mkdir -p "Likelihood_Values"
cd "Likelihood_Values" || exit 1

## -------------------------- ##
##   PARALLELISATION SET UP   ##
## -------------------------- ##

Number_Of_Models_Per_Array=$((NUMBER_OF_MODELS_TO_GENERATE / SLURM_ARRAY_TASK_COUNT))
Remainder=$((NUMBER_OF_MODELS_TO_GENERATE % SLURM_ARRAY_TASK_COUNT))

# To optimise things, if the number of models to be generated isn't a multiple of the number of arrays, the array with the smallest id will learn more models (the remaining ones), but these will all be the smaller models.
if [ "${SLURM_ARRAY_TASK_ID}" -eq 1 ]; then
    Starting_Number_Of_States=2
    Ending_Number_Of_States=$(( 2 + STATE_INCREMENT*(Remainder+Number_Of_Models_Per_Array-1) ))
else
    Starting_Number_Of_States=$(( (((SLURM_ARRAY_TASK_ID-1)*Number_Of_Models_Per_Array) + Remainder)*STATE_INCREMENT + 2 )) 
    Ending_Number_Of_States=$(( (((SLURM_ARRAY_TASK_ID)*Number_Of_Models_Per_Array) + Remainder -1 )*STATE_INCREMENT + 2 ))
fi

## --------------------- ##
##    MODEL GENERATION   ##
## --------------------- ##

module purge
module load Java

# Recreate any likelihood files if they already exist (to avoid artefacts in the plotting stage later in the pipeline).
# Job is submitted as an array so we only want to do this for one of the array tasks, not all of them.
if [[ "${SLURM_ARRAY_TASK_ID}" -eq 1 ]]; then
    rm -f "likelihood.BinSize.${BIN_SIZE}.SampleSize.${SAMPLE_SIZE}.txt"
    touch "likelihood.BinSize.${BIN_SIZE}.SampleSize.${SAMPLE_SIZE}.txt"
fi

for numstates in $(seq "$Starting_Number_Of_States" "$STATE_INCREMENT" "$Ending_Number_Of_States"); do
    echo "Learning model with: ${numstates} states..."
    # Please note that this script is as a SLURM array, and so the memory usage at this step could get out of hand if there are too many arrays and the memory usage of ChromHMM is too high for each.
    # The -nobrowser option is here as we have no need for the browser files generated by ChromHMM
    java -mx10G -jar "${CHROMHMM_MAIN_DIR}/ChromHMM.jar" LearnModel -nobrowser -b "${BIN_SIZE}" "${BINARY_DIR}" "${MODEL_DIR}" "${numstates}" hg19 > "ChromHMM.Output.BinSize.${BIN_SIZE}.numstates.${numstates}.txt"

    # Get the estimated log likelihood value from ChromHMM's output
    echo "Writing estimated log likelihood to likelihood.BinSize.${BIN_SIZE}.SampleSize.${SAMPLE_SIZE}.txt"
    echo -n "Estimated Log Likelihood for ${numstates} states: " >> "likelihood.BinSize.${BIN_SIZE}.SampleSize.${SAMPLE_SIZE}.txt"
    grep "       " "ChromHMM.Output.BinSize.${BIN_SIZE}.numstates.${numstates}.txt" | tail -1 | awk '{print $2}' >> "likelihood.BinSize.${BIN_SIZE}.SampleSize.${SAMPLE_SIZE}.txt" # grep removes the lines associated with writing to files. The tail and awk locate the final estimated log likelihood

    rm "ChromHMM.Output.BinSize.${BIN_SIZE}.numstates.${numstates}.txt"
done

## ------------------ ##
##   RENAMING FILES   ##
## ------------------ ##

cd "${MODEL_DIR}" || { echo "Model directory doesn't exist, make sure config.txt is pointing to the correct directory"; exit 1; }
Emission_Files_To_Rename=$(find . -type f -name "emissions*")
for file in $Emission_Files_To_Rename; do
    file_ending=$(echo "$file" | cut -d "_" -f 2) # obtains number of states and file extension name
    mv "$file" "Emissions_BinSize_${BIN_SIZE}_SampleSize_${SAMPLE_SIZE}_NumberOfStates_${file_ending}" 
done
Transistion_Files_To_Rename=$(find . -type f -name "transitions*")
for file in $Transistion_Files_To_Rename; do
    file_ending=$(echo "$file" | cut -d "_" -f 2)
    mv "$file" "Transitions_BinSize_${BIN_SIZE}_SampleSize_${SAMPLE_SIZE}_NumberOfStates_${file_ending}" 
done
Model_Files_To_Rename=$(find . -type f -name "model*")
for file in $Model_Files_To_Rename; do
    file_ending=$(echo "$file" | cut -d "_" -f 2)
    mv "$file" "Model_BinSize_${BIN_SIZE}_SampleSize_${SAMPLE_SIZE}_NumberOfStates_${file_ending}" 
done


## ----------------------- ##
##   LOG FILE MANAGEMENT   ##
## ----------------------- ##

# Finishing message
echo "Job completed at:"
date -u
end_time=$(date +%s)
time_taken=$((end_time-start_time))
echo "Job took a total of: ${time_taken} seconds to complete"

# Removing temporary log files
cd "${LOG_FILE_PATH}" || exit 1
rm "temp${SLURM_ARRAY_TASK_ID}.e"
rm "temp${SLURM_ARRAY_TASK_ID}.o"
