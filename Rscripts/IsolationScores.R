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
## $1 -> State assignments file                                   ||
## $2 -> Output file path                                         ||
## $3 -> Sample size for isolation score                          ||
## ============================================================== ##
## OUTPUTS:                                                       ||
## Plot of isolation score for each state                         ||
## Text file containing isolation scores for each state           ||
## ============================================================== ##

## ========== ##
##   SET UP   ##
## ========== ##

rm(list = ls())

setwd("/lustre/projects/Research_Project-MRC190311/scripts/integrative")
source("ChromOptimise/configuration/config.R")

setwd(model_dir)

library(ggplot2)

arguments <- commandArgs(trailingOnly = TRUE)
state_assignments_file <- arguments[1]
output_file_path <- arguments[2]
isolation_sample_size <- as.numeric(arguments[3])

## =================== ##
##   FILE PROCESSING   ##
## =================== ##

# Skipping first two rows of this file as the metadata is not required
state_assignments <- read.table(state_assignments_file, skip = 2)

## ============================= ##
##   ISOLATION SCORE FUNCTIONS   ##
## ============================= ##

check_upstream <- function(bin_index, distance, reference_state) {
  state <- state_assignments[bin_index - distance, 1]
  if (state == reference_state) {
    return(distance - 1)
  }
  return(check_downstream(bin_index, distance, reference_state))
}

check_downstream <- function(bin_index, distance, reference_state) {
  state <- state_assignments[bin_index + distance, 1]
  if (state == reference_state) {
    return(distance - 1)
  }
  return(check_upstream(bin_index, distance + 1, reference_state))
}

# This function utilises recursion of check_(up/down)stream to
# find the closest bin that has the same state assignment to the
# input bin index
matching_bin_distance <- function(bin_index) {
  reference_state <- state_assignments[bin_index, 1]
  distance <- check_upstream(bin_index, 1, reference_state)

  return(distance)
}

# Isolation score is the average distance between two bins that have
# the same state assignment. A higher value indicates that a state is
# more isolated in the state assignment and therefore less stable.
get_isolation_score <- function(bin_indicies) {
  sum_of_distances <- 0
  for (bin_index in bin_indicies){
    sum_of_distances <- sum_of_distances + matching_bin_distance(bin_index)
  }
  sample_size <- length(bin_indicies)
  average_distance <- (sum_of_distances / sample_size)

  return(average_distance)
}

## ========================= ##
##   SUBSAMPLING FUNCTIONS   ##
## ========================= ##

subsample_target_bins <- function(target_state, sample_percent) {
  bins_with_target_value <- which(state_assignments[, 1] == target_state)

  sample_size <- length(bins_with_target_value) * (sample_percent / 100)

  sampled_bin_indices <- sample(bins_with_target_value, sample_size)

  return(sampled_bin_indices)
}

## ======== ##
##   MAIN   ##
## ======== ##

states <- unlist(unique(state_assignments))
samples <- rep(sample_size, times = length(states))

bin_indices_sample <- mapply(subsample_target_bins, states, samples)

isolation_scores <- unlist(lapply(bin_indices_sample, get_isolation_score))

isolation_scores_output <-
  data.frame(states = states, isolation_scores = isolation_scores)

sorted_isolation_scores <-
  isolation_scores_output[order(isolation_scores_output$states), ]

## =========== ##
##   OUTPUTS   ##
## =========== ##

setwd(output_file_path)

isolation_scores_scatter <-
  ggplot(sorted_isolation_scores, aes(x = states, y = isolation_scores)) +
  geom_point() +
  scale_x_continuous(breaks = seq(min(sorted_isolation_scores$states),
                                  max(sorted_isolation_scores$states),
                                  by = 1)) +
  labs(title = "Isolation scores", x = "State number", y = "Isolation score") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(
  "Isolation_Scores.png"
)

write.table(sorted_isolation_scores, "Isolation_Scores.txt", row.names = FALSE)
