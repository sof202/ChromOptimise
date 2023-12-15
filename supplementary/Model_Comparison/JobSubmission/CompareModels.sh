#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -p mrcq # submit to the mrc queue for faster queue times
#SBATCH --time=01:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC190311 
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=16 
#SBATCH --mem=10G 
#SBATCH --mail-type=END # Send an email after the job is done
#SBATCH --output=/lustre/projects/Research_Project-MRC190311/scripts/integrative/blueprint/LogFiles/Model_Comparing/temp.o #Put output file in log files with with temporary name
#SBATCH --error=/lustre/projects/Research_Project-MRC190311/scripts/integrative/blueprint/LogFiles/Model_Comparing/temp.e #Put error file in log files with with temporary name
#SBATCH --job-name=Model_Comparing

## -------------------------------------------------------------------------------------------- ##
##                                                                                              ##
##                                            PREAMBLE                                          ##
##                                                                                              ##
## -------------------------------------------------------------------------------------------- ##
##                                            PURPOSE                                           ##
##   The models produced in 5_batch_CreateIncrementalModels.sh are compared using ChromHMM's    ##
##   CompareModels. The base model used for comparing is the most complex model. The emission   ##
##   file for the most complex model is deleted and the process is repeated for the next most   ##
##           complex model. This continues until all emission files have been deleted.          ##
## -------------------------------------------------------------------------------------------- ##
##                        AUTHOR: Sam Fletcher s.o.fletcher@exeter.ac.uk                        ##
##                                     CREATED: November 2023                                   ##
## -------------------------------------------------------------------------------------------- ##
##                                         PREREQUISITES                                        ##
##                            Run: 5_batch_CreateIncrementalModels.sh                           ##
## -------------------------------------------------------------------------------------------- ##
##                                          DEPENDENCIES                                        ##
##                                              Java                                            ##
##                                            ChromHMM                                          ##
## -------------------------------------------------------------------------------------------- ##
##                                            INPUTS                                            ##
##                                             NONE                                             ##
## -------------------------------------------------------------------------------------------- ##
##                                            OUTPUTS                                           ##
##                            Model comparison files in .txt format                             ##
## -------------------------------------------------------------------------------------------- ##

## ------------------------ ##
##    HELP FUNCTIONALITY    ##
## ------------------------ ##

if [[ "$1" == "--help" || "$1" == "-h" ]]; then
    echo "=========================================================================="
    echo "Purpose: Uses ChromHMM's CompareModels to compare model files sequentially"
    echo "Author: Sam Fletcher"
    echo "Contact: s.o.fletcher@exeter.ac.uk"
    echo "Dependencies: Java, ChromHMM"
    echo "Inputs:"
    echo "NONE"
    echo "=========================================================================="
    exit 0
fi

## ------------ ##
##    SET UP    ##
## ------------ ##

# Rename the output and error files to have format: [job id]~[date]-[time]
# This requires a hard link as you cannot rename log files whilst running the script without a wrapper function
LOG_FILE_PATH=/lustre/projects/Research_Project-MRC190311/scripts/integrative/blueprint/LogFiles/$SLURM_JOB_NAME/
cd "${LOG_FILE_PATH}" || exit 1
timestamp=$(date -u +%Y.%m.%d-%H:%M)
ln temp.e "${SLURM_JOB_ID}~$timestamp.e"
ln temp.o "${SLURM_JOB_ID}~$timestamp.o"

# Print start date/time
echo "Job '$SLURM_JOB_NAME' started at..."
date -u

# Get the start time for the program
start_time=$(date +%s)

# Activate config.txt to access all file paths
echo "Loading config file:"
source /lustre/projects/Research_Project-MRC190311/scripts/integrative/blueprint/config/config.txt

# Exit with warning if the model directory is empty
cd "${MODEL_DIR}" || { echo "Model directory doesn't exist, make sure config.txt is pointing to the correct directory"; exit 1; }
if [ -z "$(ls -A)" ]; then
    echo "5_ModelFiles is empty, ensure that 5_CreateIncrementalModels.sh has been ran before this script."
    echo "Aborting..."

    # Remove temporary log files
    cd "${LOG_FILE_PATH}" || exit 1
    rm temp.e
    rm temp.o

    exit 1
fi

# Create a temporary folder in the 6_ModelComparisonFiles directory
cd "${COMPARE_DIR}" || { echo "Comparison directory doesn't exist, make sure config.txt is pointing to the correct directory"; exit 1; }
mkdir -p temp
cd temp || exit 1
rm ./*


# Copy emission files to a temporary directory
cd "${MODEL_DIR}" || { echo "Model directory doesn't exist, make sure config.txt is pointing to the correct directory"; exit 1; }
Emission_Text_Files=$(find . -type f -name "Emission*.txt")

echo "Copying emission files to a temporary directory..."
for file in $Emission_Text_Files; do
    echo "Copying ${file}..."
    new_file_name=$(echo "$file" | tr '[:upper:]' '[:lower:]')
    # ChromHMM's CompareModels requires files to start with 'emissions', case sensitive.
    cp "$file" "${COMPARE_DIR}/temp/${new_file_name}"
done






## -------------------- ##
##   COMPARING MODELS   ##
## -------------------- ##

# Explanation of main loop below
# 1) Sort the files by the number of states
# 2) Create a comparison file using the most complex model as a base
# 3) Delete the most complex model

module purge
module load Java


for file in $Emission_Text_Files; do
    # 1)
    cd "${COMPARE_DIR}" || exit 1
    Most_Complex_Model_Number=$(find . -type f -name "emissions*.txt" | grep -oP "\d+(?=.txt)"| sort -g | tail -1) 
    # The grep command here is to extract the 'number of states' out of the emission files (all such files have the form 'emissions_x.txt', where x is the number of states).
    # The numbers are sorted geometrically and the largest value is extracted.
    Most_Complex_Model_File=$(find . -type f -name "emissions*.txt" | grep "${Most_Complex_Model_Number}.txt")
    echo "${Most_Complex_Model_File}"

    # 2)
    cd "${COMPARE_DIR}" || exit 1
    echo "Comparing ${Most_Complex_Model_File} to the present less complex emission files.."
    java -mx4G -jar "${CHROMHMM_MAIN_DIR}/ChromHMM.jar" CompareModels "${Most_Complex_Model_File}" temp/ "Comparison_To_${Most_Complex_Model_Number}_states" 

    # 3)
    rm "${Most_Complex_Model_File}"
done

rm -rf temp


## ----------------------- ##
##   LOG FILE MANAGEMENT   ##
## ----------------------- ##

#Finish message and time
echo "Job completed at:"
date -u
end_time=$(date +%s)
time_taken=$((end_time-start_time))
echo "Job took a total of: ${time_taken} seconds to complete"

# Remove temporary log files
cd "${LOG_FILE_PATH}" || exit 1
rm temp.e
rm temp.o



