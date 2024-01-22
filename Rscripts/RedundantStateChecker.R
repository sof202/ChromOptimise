## ============================================================== ##
##                                                                ||
##                            PREAMBLE                            ||
##                                                                ||
## ============================================================== ##
## PURPOSE:                                                       ||
## Using the thresholds set in the config.txt file, this script   ||
## determinesthe presence of any redundant states in a hidden     ||
## Markov model produced by ChromHMM.                             ||
## The metrics for state similarity are:                          ||
##    1) Euclidean distance between emission parameters           ||
##    2) Maximum transition probability towards states            ||
## ============================================================== ##
## AUTHOR: Sam Fletcher                                           ||
## CONTACT: s.o.fletcher@exeter.ac.uk                             ||
## CREATED: December 2023                                         ||
## ============================================================== ##
## PREREQUISITES:                                                 ||
## Run 5_batch_CreateIncrementalModels.sh                         ||
## ============================================================== ##
## INPUTS:                                                        ||
## $1 -> Model size                                               ||
## $2 -> Bin size                                                 ||
## $3 -> Sample size                                              ||
## $4 -> Directory to place output files into                     ||
## ============================================================== ##
## OUTPUTS:                                                       ||
## Boolean response for the presence of redundant states          ||
## List of states with similar emission parameters                ||
## List of states with low transition probabilities towards them  ||
## List of states with high isolation scores (or no assignment)   ||
## List of states that are considered redundant                   ||
## ============================================================== ##


## ========== ##
##   SET UP   ##
## ========== ##

rm(list = ls())

library("stringr")

setwd("/lustre/projects/Research_Project-MRC190311/scripts/integrative")
source("ChromOptimise/configuration/config.R")
setwd(model_dir)

arguments <- commandArgs(trailingOnly = TRUE)
model_size <- as.numeric(arguments[1])
bin_size <- arguments[2]
sample_size <- arguments[3]
output_file_path <- arguments[4]

euclidean_distance <- function(vector_a, vector_b) {
  sqrt(sum((vector_a - vector_b) ^ 2))
}

## ================= ##
##   LOADING FILES   ##
## ================= ##

emissions_file <- paste0("Emissions_BinSize_", bin_size, "_SampleSize_",
                         sample_size, "_States_", model_size, ".txt")

emissions_data <- read.table(emissions_file, skip = 1)
emissions_data <- subset(emissions_data, select = -V1)


transitions_file <- paste0("Transitions_BinSize_", bin_size, "_SampleSize_",
                           sample_size, "_States_", model_size, ".txt")

transitions_data <- read.table(transitions_file, skip = 1)
transitions_data <- subset(transitions_data, select = -V1)

# Isolation scores are already in the output file path from IsolationScores.R
setwd(output_file_path)
isolation_data <- read.table("Isolation_Scores.txt", header = TRUE)

## =================================== ##
##   EUCLIDEAN DISTANCE CALCULATIONS   ##
## =================================== ##

euclidean_distances <- data.frame(state_pair = numeric(),
                                  euclidan_distance = double())

for (reference_state_index in 1:(model_size - 1)){
  for (comparison_state_index in (reference_state_index + 1):model_size){

    reference_state <- as.numeric(emissions_data[reference_state_index, ])
    comparison_state <- as.numeric(emissions_data[comparison_state_index, ])

    distance <- euclidean_distance(reference_state, comparison_state)

    state_pair <- paste0(reference_state_index, " ", comparison_state_index)

    euclidean_distances[nrow(euclidean_distances) + 1, ] <-
      c(state_pair, distance)
  }
}

## ================================================= ##
##   MAXIMUM TRANSITION PROBABILITIES CALCULATIONS   ##
## ================================================= ##

max_transition_towards_states <- data.frame(state = integer(),
                                            maximum_probability = double())

# Maximum value in each column gives the
# maximum transition probability towards each state
for (state in 1:model_size){
  max_transition_probability <- max(transitions_data[, state])

  max_transition_towards_states[nrow(max_transition_towards_states) + 1, ] <-
    c(state, max_transition_probability)
}

## ====================================== ##
##   IDENTIFICATION OF REDUNDANT STATES   ##
## ====================================== ##

## Similar state pairings ##
low_euclidean_distances <-
  euclidean_distances[euclidean_distances[, 2]
                      < emissions_threshold, ]

# Extract states from state pairs in low euclidean distances
similar_state_pairs <- low_euclidean_distances[, 1]
similar_states <- c()
for (state_pair in similar_state_pairs){
  similar_states <- append(similar_states, substring(state_pair, 1, 1))
  similar_states <- append(similar_states, substring(state_pair, 3, 3))
}
similar_states_list <- unique(similar_states)


## States with low max transition probability towards them ##
low_transition_probabilites <-
  max_transition_towards_states[max_transition_towards_states[, 2]
                                < transitions_threshold, ]


low_transition_states <- low_transition_probabilites[, 1]


## States with a high isolation score or no isolation score ##
isolated_states_data <- isolation_data[isolation_data$isolation_scores > isolation_threshold & !is.na(isolation_data$isolation_scores), ]

isolated_states <- isolated_states_data[, 1]

# Finds all states that have no isolation score (implying they are unassigned)
unassigned_states <- setdiff((1:model_size), isolation_data[[1]])
isolated_states <- append(isolated_states, unassigned_states)

# Finds all states that were only assigned once (resulting in NA)
single_assigned_states <- 
  isolation_data$states[isolation_data$isolation_scores %in% NA]
isolated_states <- append(isolated_states, single_assigned_states)


## Check for redundant states by using critereon ##
# (i) Similar state pair existence
# (ii) Low maximum transition probability
# (iii) High isolation
redundant_states <- c()
for (state in low_transition_states){
  if (state %in% similar_states_list & state %in% isolated_states) {
    redundant_states <- append(redundant_states, state)
  }
}

## =========== ##
##   OUTPUTS   ##
## =========== ##

output_file <- paste0("Redundant_States_Modelsize_", model_size, ".txt")
setwd(output_file_path)
separator <- "<------------------------------------------------------------>"

write("States with similar emission probabilities:\n",
      file = output_file)
write(separator, file = output_file, append = TRUE)

write("|Similar state pairs| |Euclidean distance between states|",
      file = output_file, append = TRUE)
write(separator, file = output_file, append = TRUE)

write.table(low_euclidean_distances, file = output_file,
            append = TRUE, row.names = FALSE, col.names = FALSE)
write(separator, file = output_file, append = TRUE)


write("\nStates with low transition probability towards them:\n"
      , file = output_file, append = TRUE)
write(separator, file = output_file, append = TRUE)

write("|State| |Maximum probability of transitioning towards state|",
      file = output_file, append = TRUE)
write(separator, file = output_file, append = TRUE)

write.table(low_transition_probabilites, file = output_file,
            append = TRUE, row.names = FALSE, col.names = FALSE)
write(separator, file = output_file, append = TRUE)


write("\nStates with high isolation score:\n"
      , file = output_file, append = TRUE)
write(separator, file = output_file, append = TRUE)

write.table(isolated_states_data, file = output_file,
            append = TRUE, row.names = FALSE, col.names = FALSE)
write(separator, file = output_file, append = TRUE)

write("\nStates with no/single assignment:\n"
      , file = output_file, append = TRUE)
write(separator, file = output_file, append = TRUE)

write.table(single_assigned_states, file = output_file,
            append = TRUE, row.names = FALSE, col.names = FALSE)
write.table(unassigned_states, file = output_file,
            append = TRUE, row.names = FALSE, col.names = FALSE)
write(separator, file = output_file, append = TRUE)


write("\nDetermined redundant states:", file = output_file, append = TRUE)
if (length(redundant_states) == 0) {
  write("NONE", file = output_file, append = TRUE)
} else {
  write(redundant_states, file = output_file, append = TRUE)
}
