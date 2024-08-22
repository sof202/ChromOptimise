#!/bin/bash
# Export all environment variables to the batch job
#SBATCH --export=ALL 
# Submit to the mrc queue for faster queue times
#SBATCH -p mrcq 
# Script usually takes less than one minute
#SBATCH --time=00:10:00 
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
#SBATCH --job-name=3_Optimal_States

usage() {
cat <<EOF
===========================================================================
3_OptimalNumberOfStates
===========================================================================
Purpose: Determines the optimum number of states to use with your dataset
Author: Sam Fletcher
Contact: s.o.fletcher@exeter.ac.uk
Dependencies: R, gzip
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

if [[ -z "${CHROMOSOME_IDENTIFIER}" ]]; then
    CHROMOSOME_IDENTIFIER=1
    echo "No chromosome identifier was given, using the default value of:" \
    "${CHROMOSOME_IDENTIFIER} instead."
fi

## ================== ##
##   ERROR CATCHING   ##
## ================== ##

input_directory="${MODEL_DIR}/BinSize_${BIN_SIZE}_models_${NUMBER_OF_MODELS}"

if [[ -z "$(ls -A "${input_directory}")" ]]; then
    { >&2 echo -e "ERROR: No files found in: ${input_directory}.\n"\
    "Please run 2_CreateIncrementalModels.sh before this script."
    finishing_statement 1; }
fi

## =============== ##
##    MAIN LOOP    ##
## =============== ##

# Using the emission files here for the model numbers is arbitrary, model
# or transition files could have been used just the same.
# We reverse the order as we plan on working from the most complex model
# to the least until a model with no redundant states is found.
mapfile -t model_sizes < \
    <(find "${input_directory}" -type f -name "emissions*.txt" | \
    grep -oP "\d+(?=.txt)" | \
    sort -gr)

output_directory="${OPTIMUM_STATES_DIR}/BinSize_${BIN_SIZE}_models_${NUMBER_OF_MODELS}"

# Clear up the directory in case of repeat runs (with different thresholds)
rm -rf "${output_directory}"

mkdir -p \
    "${output_directory}/Euclidean_distances" \
    "${output_directory}/Flanking_states" \
    "${output_directory}/Isolation_scores"

module purge
module load "${R_MODULE}"

cd "${RSCRIPTS_DIR}" || \
{ >&2 echo "ERROR: [\${RSCRIPTS_DIR} - ${RSCRIPTS_DIR}] doesn't exist, \
make sure Config.txt is pointing to the correct directory"
finishing_statement 1; }

for model_number in "${model_sizes[@]}"; do
    # State assignments are named:
    # CellType_BinSize_ModelSize_Chromosome_statebyline.txt

    # We only look at one chromosome as the decision of how to handle the
    # following case is rather arbitrary (see wiki): 
    # "A state is not assigned on one chromosome but has dense assignment
    # on another"
    state_assignment_file=$(\
    find "${input_directory}" \
    -iname "${CELL_TYPE}*_${model_number}_chr${CHROMOSOME_IDENTIFIER}_*" \
    )
    emissions_file=$(\
    find "${input_directory}" -name "emissions_${model_number}.txt*" \
    )
    transitions_file=$(\
    find "${input_directory}" -name "transitions_${model_number}.txt*" \
    )

    # The existence of this file (provided the user doesn't delete them)
    # implies the existence of the emissions and transistions files.
    if [[ -z "$state_assignment_file" ]]; then
        { >&2 echo "ERROR: No state assignment file found for chromosome:" \
        "${CHROMOSOME_IDENTIFIER}, please check ${input_directory}/STATEBYLINE for" \
        "the existence of this state assignment file"; finishing_statement 1; }
    fi

    echo "Running SimilarEmissions.R for: ${model_number} states..."

    Rscript SimilarEmissions.R \
        "${emissions_file}" \
        "${output_directory}/Euclidean_distances" \
        FALSE

    echo "Running FlankingStates.R for: ${model_number} states..."

    Rscript FlankingStates.R \
        "${transitions_file}" \
        "${output_directory}/Flanking_states"

    echo "Running IsolationScores.R for: ${model_number} states..."

    Rscript IsolationScores.R \
        "${state_assignment_file}" \
        "${output_directory}/Isolation_scores" \
        "${model_number}" \
        100 \
        FALSE

    echo "Running RedundantStateChecker.R for: ${model_number} states..."

    Rscript RedundantStateChecker.R \
        "${configuration_directory}/Config.R" \
        "${model_number}" \
        "${output_directory}"

    redundant_states_found=$(tail -1 \
    "${output_directory}/Redundant_states_model-${model_number}.txt")

    if [[ "${redundant_states_found}" == "NONE" ]]; then
        echo "Model with ${model_number} states has no redundant states." >> \
        "${output_directory}/OptimumNumberOfStates.txt"
        break
    else
        echo -n "Model with ${model_number} states has redundant state(s):" >> \
        "${output_directory}/OptimumNumberOfStates.txt"

        echo "${redundant_states_found}" >> \
        "${output_directory}/OptimumNumberOfStates.txt"
    fi
done

optimum_states="${model_number}"

## ========================= ##
##   OPTIMUM STATES CHECK    ##
## ========================= ##

# If the largest model learned has no redundant states, this doesn't necessarily
# imply that it has the optimum number of states, perhaps a more complex model
# does. This section checks for this scenario.

if [[ $(wc -l < "${output_directory}/OptimumNumberOfStates.txt") -eq 1 ]]; then
cat >> "${output_directory}/OptimumNumberOfStates.txt" << EOF 
WARNING: Largest model learned has no redundant states.
${optimum_states} states may not be the optimum number of states.
Try increasing the size of the most complex model or increasing 
the thresholds in the Config.R file.
EOF
else
    echo "Optimum number of states for the data is: ${optimum_states}" >> \
    "${output_directory}/OptimumNumberOfStates.txt"
fi


next_smallest_model_number=$((optimum_states - 1))

# There is a small chance that the optimum states is 2, in which case no model
# with one state exists (the below would fail in this case)
if [[ "${next_smallest_model_number}" -gt 1 ]]; then
smaller_overlap_file=$(\
    find "${input_directory}" \
    -name "*${next_smallest_model_number}_overlap.txt" \
)
smaller_emissions_file=$(\
    find "${input_directory}" \
    -name "emissions_${next_smallest_model_number}.txt*" \
)
optimum_overlap_file=$(\
    find "${input_directory}" \
    -name "*${optimum_states}_overlap.txt" \
)
optimum_emissions_file=$(\
    find "${input_directory}" \
    -name "emissions_${optimum_states}.txt*" \
)

are_overlaps_a_subset=$(\
    "${RSCRIPTS_DIR}/CheckModelSubset.r" \
    "${configuration_directory}/config.r" \
    "${smaller_overlap_file}" \
    "${optimum_overlap_file}" \
    "overlap" \
)

are_emissions_a_subset=$(\
    "${RSCRIPTS_DIR}/CheckModelSubset.r" \
    "${configuration_directory}/config.r" \
    "${smaller_emissions_file}" \
    "${optimum_emissions_file}" \
    "emissions" \
)

if [[ "${are_overlaps_a_subset}" =~ "FALSE" ]] || [[ "${are_emissions_a_subset}" =~ "FALSE" ]]; then
cat >> "${output_directory}/OptimumNumberOfStates.txt" << EOF
WARNING: A model that is smaller than the determined 'optimum' model
($next_smallest_model_number states) has states that are not included in the
determined optimum model ($optimum_states states). Smaller models should be
subsets of larger models in general. If a smaller model has different states to
the model with the 'optimum' number of states, then the 'optimal' model misses
crucial information (despite its higher complexity). 

It is reccomended to look at both of these models and determine if you want to
keep the more complex model, or use the smaller model in subsequent analysis.
You can manually change the OPTIMUM_NUMBER_OF_STATES in Config.txt.
EOF
fi
fi


## ============== ##
##    PLOTTING    ##
## ============== ##

echo "Plotting the estimated log likelihoods of learned models against" \
"one another..."

Rscript PlotLikelihoods.R \
    "${input_directory}/Likelihood_Values/likelihoods.txt" \
    "${output_directory}"

# We plot the relative Bayesian information critereon for each model so the user
# can use a less complex model if they wish. A less complex model may have a
# similar BIC to a more complex model (meaning its roughly as accurate).
# This is up to the user's discretion as BIC is just a heuristic. 

# BIC requires the number of observations, which is the total number of lines  
# in each binary file minus 2 times the number of files 
# (due to the headers in each file).
full_binary_path="${BINARY_DIR}/BinSize_${BIN_SIZE}"

total_observations=0
for file in "${full_binary_path}"/*_binary.txt*; do
    observations=$(gzip -dc "$file" | wc -l)
    total_observations=$((total_observations + observations - 2))
done

optimum_number_of_states=$(\
tail -1 "${output_directory}/OptimumNumberOfStates.txt" | \
cut -d: -f2 | \
tr -d ' ')

echo "Processing the Bayesian information critereon of learned models..."
Rscript CalculateBIC.R \
    "${configuration_directory}/Config.R" \
    "${input_directory}/Likelihood_Values/likelihoods.txt" \
    "${optimum_number_of_states}" \
    "${total_observations}" \
    "${output_directory}" \

cp "${input_directory}/model_${optimum_number_of_states}.txt" \
    "${output_directory}/optimum_model.txt"

finishing_statement 0

