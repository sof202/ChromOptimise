## ============================================================== ##
##                                                                ||
##                            PREAMBLE                            ||
##                                                                ||
## ============================================================== ##
## PURPOSE:                                                       ||
## Calculates the an estimated 'isolation score' for each state   ||
## in the state assignment created by ChromHMM. 'Isolation score' ||
## is given by 'the average number of bins separating two bins    ||
## with the same state assignment'                                ||
## ============================================================== ##
## AUTHOR: Sam Fletcher                                           ||
## CONTACT: s.o.fletcher@exeter.ac.uk                             ||
## CREATED: January 2022                                          ||
## ============================================================== ##
## PREREQUISITES:                                                 ||
## Run 5_batch_CreateIncrementalModels.sh                         ||
## ============================================================== ##
## INPUTS:                                                        ||
## $1 -> Location of configuation file                            ||
## $2 -> State assignments file                                   ||
## $3 -> Output file path                                         ||
## $4 -> Sample size for isolation score                          ||
## ============================================================== ##
## OUTPUTS:                                                       ||
## Text file containing isolation scores for each state           ||
## ============================================================== ##

## ========== ##
##   SET UP   ##
## ========== ##

rm(list = ls())

library(ggplot2)

arguments <- commandArgs(trailingOnly = TRUE)
config_file_location <- arguments[1]
state_assignments_file <- arguments[2]
output_file_path <- arguments[3]
isolation_sample_size <- as.numeric(arguments[4])

source(config_file_location)

## =================== ##
##   FILE PROCESSING   ##
## =================== ##

# Skipping first two rows of this file as the metadata is not required
state_assignments <- read.table(state_assignments_file, skip = 2)

# This variable is required for a check in check_downstream to avoid
# unwanted behaviour when going out of range
max_bin_index <- nrow(state_assignments)

## ============================= ##
##   ISOLATION SCORE FUNCTIONS   ##
## ============================= ##

matching_bin_distance <- function(reference_bin_index) {
  reference_state <- state_assignments[reference_bin_index, 1]

  # The maximum possible distance between two bins with the same state is
  # the number of bins between the reference bin index and the two ends of the
  # dataset
  max_distance <-
    max(max_bin_index - reference_bin_index, reference_bin_index - 1)

  for (distance in 1:max_distance) {
    # Check upstream of reference bin index

    # If the comparison bin index is less than 1 it is out of range
    if (reference_bin_index - distance >= 1) {
      comparison_state <- state_assignments[reference_bin_index - distance, 1]
      if (comparison_state == reference_state) {
        # we return distance - 1 as we want the number of bins that
        # separate the two bins with the same assignment (if they
        # are adjacent, this value should be 0, not 1)
        return(distance - 1)
      }
    }

    # Check downstream of reference bin index

    # If the comparison bin index is greater than the max bin index
    # it is out of range
    if (reference_bin_index + distance <= max_bin_index) {
      comparison_state <- state_assignments[reference_bin_index + distance, 1]
      if (comparison_state == reference_state) {
        return(distance - 1)
      }
    }
  }

  # Error catching just in case no matching bin is found
  return(NA)
}

# Isolation score is the average distance between two bins that have
# the same state assignment. A higher value indicates that a state is
# more isolated in the state assignment and therefore less stable.
get_isolation_score <- function(bin_indicies) {
  # Check for if the state is only assigned once in the state assignment
  if (all(is.na(bin_indicies))) {
    return(NA)
  }

  sum_of_distances <- 0
  for (bin_index in bin_indicies){
    sum_of_distances <- sum_of_distances + matching_bin_distance(bin_index)
  }

  average_distance <- (sum_of_distances / length(bin_indicies))

  return(average_distance)
}

## ========================= ##
##   SUBSAMPLING FUNCTIONS   ##
## ========================= ##

subsample_target_bins <- function(target_state, sample_percent) {
  bins_with_target_value <- which(state_assignments[[1]] == target_state)

  # If only one bin has been assigned with the target state the recursion
  # will not halt. This is a flag to prevent this scenario.
  if (length(bins_with_target_value) == 1) {
    return(NA)
  }

  # We want a representative sample for each state, some states will be
  # assigned to a much larger number of bins than others
  sample_size <- ceiling(length(bins_with_target_value) * (sample_percent / 100))

  if (sample_size > length(bins_with_target_value)) {
    sample_size <- length(bins_with_target_value)
  }
  sampled_bin_indices <- sample(bins_with_target_value, sample_size)

  return(sampled_bin_indices)
}

## ======== ##
##   MAIN   ##
## ======== ##

states <- unlist(unique(state_assignments))
samples <- rep(isolation_sample_size, times = length(states))

bin_indices_sample <- mapply(subsample_target_bins, states, samples)

isolation_scores <- unlist(lapply(bin_indices_sample, get_isolation_score))

isolation_scores_output <-
  data.frame(state = states, isolation_score = isolation_scores)

# This is purely so that the output text file is easier to read
sorted_isolation_scores <-
  isolation_scores_output[order(isolation_scores_output$states), ]

## =========== ##
##   OUTPUTS   ##
## =========== ##

setwd(output_file_path)

output_file_name <- paste0("Isolation_Scores_model-",model_size,".txt")

write.table(sorted_isolation_scores,
            output_file_name,
            row.names = FALSE)