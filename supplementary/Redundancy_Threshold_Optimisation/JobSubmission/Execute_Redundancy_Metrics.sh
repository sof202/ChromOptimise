#!/bin/bash
# Export all environment variables to the batch job
#SBATCH --export=ALL 
# Submit to the mrc queue for faster queue times
#SBATCH -p mrcq 
# Script usually takes less than three minutes
#SBATCH --time=00:10:00 
#SBATCH -A Research_Project-MRC190311
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
# Very low memory usage is expected from generating the plots
#SBATCH --mem=1G 
#SBATCH --mail-type=END # Send an email after the job is done
# Temporary log file, later to be removed
#SBATCH --output=temp%j.log
# Temporary error file, later to be removed
#SBATCH --error=temp%j.err
#SBATCH --job-name=Redundancy_Threshold_Plotting

## =================================================================================##
##                                                                                  ||
##                                     PREAMBLE                                     ||
##                                                                                  ||
## =================================================================================##
## PURPOSE:                                                                         ||
## Generate plots for the euclidean distances between the emission parameters (for  ||
## pairs of states) and the maximum transition probability towards each state in    ||
## the selected model.                                                              ||
## =================================================================================##
## AUTHOR: Sam Fletcher                                                             ||
## CONTACT: s.o.fletcher@exeter.ac.uk                                               ||
## CREATED: December 2023                                                           ||
## =================================================================================##
## PREREQUISITES: Run Generate_Big_Model.sh or ChromHMM's LearnModel command        ||
## =================================================================================##
## DEPENDENCIES: R                                                                  ||
## =================================================================================##
## INPUTS:                                                                          ||
## -c|--config= -> Full/relative file path for configuation file directory          
## -n|--size=   -> Size of model (default: 20)
## -s|--seed=   -> Random seed (default: 1)
## -o|--output  -> Path to directory containing the model files
##                 (default: \${BIG_MODELS_DIR} in FilePaths.txt)
## =================================================================================##
## OUTPUTS:                                                                         ||
## Histogram plot of Euclidean distances between emission parameters of pairs       ||
## of states. A threshold suggestion value is also given.                           ||
## A scatter plot of the maximum transition probabilitites towards each state.      ||
## =================================================================================##

## ===================== ##
##   ARGUMENT PARSING    ##
## ===================== ##

usage() {
cat <<EOF
=======================================================================
CompareModels
=======================================================================
Purpose: Generates plots to aid in thresholds used in config.R which
are used in 6_OptimumNumberOfStates.sh in determining redundant states.
Author: Sam Fletcher
Contact: s.o.fletcher@exeter.ac.uk
Dependencies: R
Inputs:
-c|--config= -> Full/relative file path for configuation file directory
-n|--size=   -> Size of model (default: 20)
-s|--seed=   -> Random seed (default: 1)
-o|--output  -> Path to directory containing the model files
                (default: \${BIG_MODELS_DIR} in FilePaths.txt)
=======================================================================
EOF
    exit 0
}

needs_argument() {
    # Required check in case user uses -a -b or -b -a (no argument given).
    if [[ -z "$OPTARG" || "${OPTARG:0:1}" == - ]]; then usage; fi
}

while getopts c:n:s:o:-: OPT; do
    # Adds support for long options by reformulating OPT and OPTARG
    # This assumes that long options are in the form: "--long=option"
    if [ "$OPT" = "-" ]; then
        OPT="${OPTARG%%=*}"
        OPTARG="${OPTARG#"$OPT"}"
        OPTARG="${OPTARG#=}"
    fi
    case "$OPT" in
        c | config )  needs_argument; configuration_directory="$OPTARG" ;;
        n | size )    needs_argument; model_size="$OPTARG" ;;
        s | seed )    needs_argument; seed="$OPTARG" ;;
        o | ouput )   needs_argument; model_file_dir="$OPTARG" ;;
        \? )          usage ;;  # Illegal short options are caught by getopts
        * )           usage ;;  # Illegal long option
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


# Output and error files renamed to:
# ModelSize-[model size]~[job id]~[date]-[time]

mv "${SLURM_SUBMIT_DIR}/temp${SLURM_JOB_ID}.log" \
"${LOG_FILE_PATH}/ModelSize-$2~${SLURM_JOB_ID}~${timestamp:=}.log"
mv "${SLURM_SUBMIT_DIR}/temp${SLURM_JOB_ID}.err" \
"${LOG_FILE_PATH}/ModelSize-$2~${SLURM_JOB_ID}~$timestamp.err"

## =============== ##
##    VARIABLES    ##
## =============== ##

## ====== DEFAULTS ====================================================================
if [[ -z "$model_file_dir" ]]; then
    model_file_dir="${BIG_MODELS_DIR}"
    echo "Model file directory was not given, using the default of: ${BIG_MODELS_DIR}"
fi

if ! [[ "${model_size}" =~ ^[0-9]+$ ]]; then
    model_size=20
    echo "Model size given is invalid, using default value of: ${model_size}."
fi

if ! [[ "$seed" =~ ^[0-9]+$ ]]; then
    seed=1
    echo "Random seed given is invalid, using defualt value of: ${seed}." 
fi
# =====================================================================================

## ================== ##
##   FILE EXISTANCE   ##
## ================== ##

cd "${model_file_dir}" ||  { >&2 echo "ERROR: ${model_file_dir} doesn't exist, \
ensure that the directory exists before running this script."; finishing_statement 1; }

if [[ -z $(find . -type f -name "emissions*") ]]; then
    { >&2 echo -e "ERROR: No model files were found in ${model_file_dir}.\n"\
    "Ensure that you have ran Generate_Big_Model.sh or ChromHMM's "\
    "LearnModel command before using this script."; finishing_statement 1; }  
fi

## ======== ##
##   MAIN   ##
## ======== ##

module purge
module load R/4.2.1-foss-2022a

cd "${SUPPLEMENTARY_DIR}/Redundancy_Threshold_Optimisation/Rscripts" ||  \
{ >&2 echo "ERROR: ${SUPPLEMENTARY_DIR}/Redundancy_Threshold_Optimisation/Rscripts \
doesn't exist, make sure [\${SUPPLEMENTARY_DIR} - ${SUPPLEMENTARY_DIR}] \
in FilePaths.txt is pointing to the correct directory"; finishing_statement 1; }

Rscript HistogramPlotForEuclideanDistances.R "${configuration_directory}/config.R" \
"${model_size}" "${seed}" "${model_file_dir}" 

Rscript ScatterPlotForTransitionMaxima.R "${configuration_directory}/config.R" \
"${model_size}" "${seed}" "${model_file_dir}"

cd "${RSCRIPTS_DIR}" || \
{ >&2 echo "ERROR: make sure [\${RSCRIPTS_DIR} - ${RSCRIPTS_DIR}] \
in FilePaths.txt is pointing to the correct directory"; finishing_statement 1; }

mkdir -p "${model_file_dir}/IsolationScores"
state_assignment_file=$(find "${model_file_dir}" -name "*${model_size}_${seed}_chr1_*")

Rscript IsolationScores.R "${configuration_directory}/config.R" \
"${state_assignment_file}" "${model_file_dir}/IsolationScores" 100


finishing_statement 0
