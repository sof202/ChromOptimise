## ============================================================== ##
##                                                                ||
##                            PREAMBLE                            ||
##                                                                ||
## ============================================================== ##
## PURPOSE:                                                       ||
## Using the thresholds set in the config.txt file, this script   ||
## determines whether or not the inputted HMM has redundant       ||
## states.                                                        ||
## The criteria for a redundant state is:                         ||
##   (i) There exists another state with similar emission         ||
##       parameters (under Euclidean distance)                    ||
##   (ii) At least one of the following:                          ||
##      (a) The state is highly isolated (under isolation score)  ||
##      (b) The similar states are likely to be flanked by the    ||
##          same states                                           ||
##      (c) At least one of the similar states is flanked on both ||
##          sides by the other state (in the similar state pair)  ||
##                                                                ||
## In the event that (i) and (ii)(b) are triggered, the state     ||
## with the higher isolation score is assumed to be the redundant ||
## state (though this choice is inconsequential)                  ||
## ============================================================== ##
## AUTHOR: Sam Fletcher                                           ||
## CONTACT: s.o.fletcher@exeter.ac.uk                             ||
## CREATED: February 2023                                         ||
## ============================================================== ##
## PREREQUISITES:                                                 ||
## Run 5_batch_CreateIncrementalModels.sh                         ||
## ============================================================== ##
## INPUTS:                                                        ||
## $1 -> Location of configuation file                            ||
## $2 -> Model size                                               ||
## $3 -> Directory containing redundancy criteria files           ||
## ============================================================== ##
## OUTPUTS:                                                       ||
## Boolean response for the presence of redundant states          ||
## List of states that are considered redundant                   ||
## ============================================================== ##

## ========== ##
##   SET UP   ##
## ========== ##

rm(list = ls())

arguments <- commandArgs(trailingOnly = TRUE)
config_file_location <- arguments[1]
model_size <- as.numeric(arguments[2])
output_file_path <- arguments[3]

source(config_file_location)
setwd(output_file_path)

## ============================== ##
##   COLLATE REDUNDANCY METRICS   ##
## ============================== ##

euclidean_distances_file <-
  paste0(output_file_path, "/Euclidean_distances/",
         "Euclidean_distances_model-", model_size, ".txt")

euclidean_distances <- read.table(euclidean_distances_file, header = TRUE)


flanking_states_file <-
  paste0(output_file_path, "/Flanking_states/",
         "Likeliest_flanking_states_model-", model_size, ".txt")

flanking_data <- read.table(flanking_states_file, header = TRUE)


isolation_scores_file <-
  paste0(output_file_path, "/Isolation_scores/",
         "Isolation_Scores_model-", model_size, ".txt")

isolation_scores <- read.table(isolation_scores_file, header = TRUE)

## ============= ##
##   FUNCTIONS   ##
## ============= ##

# Using the duplicated() function doesn't cut it for this script as we
# want to know the pairs of states that have the same flanks, not just which
# states have the same flanks as another state
find_flanks_pairs <- function(flanking_data, for_output = FALSE) {
  same_flanks <-
    data.frame(reference_state = numeric(), comparison_state = numeric())

  for (reference_state in 1:(nrow(flanking_data) - 1)) {
    for (comparison_state in (reference_state + 1):nrow(flanking_data)) {
      reference_state_flanks <- unlist(flanking_data[reference_state, 2:3])
      comparison_state_flanks <- unlist(flanking_data[comparison_state, 2:3])

      if (identical(reference_state_flanks, comparison_state_flanks)) {
        same_flanks[nrow(same_flanks) + 1, ] <-
          c(reference_state, comparison_state)
      }
    }
  }

  if (for_output) {
    same_flanks <-
      cbind(same_flanks, flanking_data[same_flanks$reference_state, 2:3])
  }

  return(same_flanks)
}

# This is for assessing states that satisfy (ii)(c)
find_identical_flanks <- function(flanking_data) {
  identical_flank_rows <- 
    apply(flanking_data, 1, function(row) row[2] == row[3])
  
  # We only need columns 1 and 2 here as columns 2 and 3 have the same data
  # (from the logic above). Adding column 3 adds no extra information.
  identical_flanks <- flanking_data[identical_flank_rows, 1:2]
  
  return(identical_flanks)
}

# This is required as a state might satisfy (i) and (ii)(c), but the order
# of the columns are the incorrect way around.
# Suppose that state 2 is likely to be flanked on both sides by state 1 and 
# states 1 and 2 have similar emission vectors. Then simply checking if the rows 
# of the identical flanks and the similar state pairs are identical is not 
# enough. The rows in this case will be c(2,1) and c(1,2) which are not 
# identical. They will be identical after sorting however.
check_rows_have_same_values <- function(row1, row2) {
  row1_sorted <- sort(row1)
  row2_sorted <- sort(row2)
  
  return(identical(row1_sorted, row2_sorted))
}

# A simple merge does not work for finding states that satisfy (i) and (ii)(c)
# This is because the order of the rows could be the wrong way around
# We could just sort the columns first and use a merge function, but this has
# the possibility of losing the information of "which state is flanking which"
# (as the order of the rows could have flipped).
find_equivalent_rows <- function(identical_flanks, similar_state_pairs) {
  equivalent_rows <- data.frame( state = numeric(), flank = numeric())

  for (row1 in 1:nrow(identical_flanks)) {
    reference_row <- unlist(identical_flanks[row1, ])
    for (row2 in 1:nrow(similar_state_pairs)) {
      comparison_row <- unlist(similar_state_pairs[row2, ])

      if (check_rows_have_same_values(reference_row, comparison_row)) {
        equivalent_rows[, nrow(equivalent_rows) + 1] <- reference_row
      }
    }
  }
}

## ======================================== ##
##   DETERMINE REDUNDANT STATE CANDIDATES   ##
## ======================================== ##

# By candidate, we mean a state that falls under any of the categories set out
# in the preamble. This is an intermediate step required before we can assess
# which states are redundant under the given criterion.

## Similar emission parameters ##
low_euclidean_distances <-
  euclidean_distances[euclidean_distances$euclidean_distance <
                      emissions_threshold, ]

similar_state_pairs <- low_euclidean_distances[, 1:2]

# When paired with the isolation score, we only care about whether or not
# a state exists in a similar state pair (not which one it is a part of)
states_in_similar_pairs <-
  unique(c(low_euclidean_distances[, 1], low_euclidean_distances[, 2]))


## State pairs with the same flanking states ##
same_flank_pairs <- find_flanks_pairs(flanking_data)

same_flank_pairs_for_output <-
  find_flanks_pairs(flanking_data, for_output = TRUE)


## States with identical flanking states ##
identical_flanks <- find_identical_flanks(flanking_data)


## High isolation scores ##
# Isolation score will be NA if the state in question was only assigned
# once. We need to account for these separately.
highly_isolated_states_data <-  isolation_scores[
  isolation_scores$isolation_score > isolation_threshold &
    !is.na(isolation_scores$isolation_score),
]

highly_isolated_states <- highly_isolated_states_data$state

# States that do not have an isolation score are those that were never
# used in the state assignment.
unassigned_states <- setdiff((1:model_size), isolation_scores$state)
isolated_states <- append(highly_isolated_states, unassigned_states)

# States are assigned NA in IsolationScores.R if they only appear once in the
# state assignment.
single_assigned_states <-
  isolation_scores$state[is.na(isolation_scores$isolation_score)]
isolated_states <- append(isolated_states, single_assigned_states)

## ================================== ##
##   REDUNDANT STATE IDENTIFICATION   ##
## ================================== ##

## The below is using the criterion set out in the preamble ##

## Accounts for (i) and (ii)(b) being satisfied ##

# Due to the way these two dataframes have been created, the entry in
# reference_state will be strictly smaller than the one in comparison_state.
# Therefore, simply merging the two dataframes will find all state pairs
# that satisfy (i) and (ii)(b) simultaneously.
redundant_state_candidates <- 
  merge(similar_state_pairs, same_flank_pairs,
        by = c("reference_state", "comparison_state"))[, 1:2]

# This function is in place so that we can choose the state that has a higher
# isolation score (as this state is more likely to be the redundant one).
# Note that we only care if a model HAS redundant states, not which states are
# redundant. We could just as easily pick the first state every time.
if (nrow(redundant_state_candidates) > 0) {
  # If no candidates were found, this code will fail as there are no rows
  # to compare. The above if statement accounts for this scenario.
  redundant_states <- unlist(apply(redundant_state_candidates, 1, function(row) {
    if (isolation_scores$isolation_score[row[1]] >
          isolation_scores$isolation_score[row[2]]) {
      return(row[1])
    }
    return(row[2])
  }))
} else {
  redundant_states <- c()
}

## Accounts for (i) and (ii)(c) being satisfied ##

# If rows in identical_flanks and similar_state_pairs have the same entries
# (either way around), then this signifies that the state in identical_flanks
# satisfies (i), as it is a part of a similar state pair. It also satisfies
# (ii)(c), as it is flanked on both sides by the same state that it is in
# a similar state pair with
identical_flank_similar_state <- 
  find_equivalent_rows(identical_flanks, similar_state_pairs)

# The state column of this dataframe holds the states that are BEING flanked
# on both sides by the state that is similar to them. This makes them
# redundant.
redundant_states <-
  append(redundant_states, identical_flank_similar_state$state)


## Accounts for (i) and (ii)(a) being satisfied ##
redundant_states <-
  append(redundant_states, intersect(isolated_states, states_in_similar_pairs))

# Some redundant states might be added multiple times from our criteria
# we only find the unique ones to help with readability of the output file.
redundant_states <- unique(redundant_states)

## =========== ##
##   OUTPUTS   ##
## =========== ##

# This is pretty gross looking, but one look at the actual output file should
# clear things up

output_file <- paste0("Redundant_states_model-", model_size, ".txt")
setwd(output_file_path)
separator <- "<------------------------------------------------------------>"

write_table_to_output <- function(table) {
  write.table(table, file = output_file, append = TRUE,
              row.names = FALSE, col.names = FALSE)
  write(separator, file = output_file, append = TRUE)
}

write_text_to_output <- function(text) {
  write(text, file = output_file, append = TRUE)
  write(separator, file = output_file, append = TRUE)
}

## Similar emission parameters ##
write_text_to_output("States with similar emission probabilities:\n")

write_text_to_output("|State pair| |Euclidean distance between states|")

write_table_to_output(low_euclidean_distances)

## State pairs wihth the same flanking states ##

write_text_to_output("\nState pairs with the same flanking states:\n")

write_text_to_output("|State pair| |Upstream flank| |Downstream flank|")

write_table_to_output(same_flank_pairs_for_output)

## States with identical flanking states ##

write_text_to_output("\nStates with identical flanking states:\n")

write_text_to_output("|State| |Flanking state|")

write_table_to_output(identical_flanks)

## Highly isolated states ##
write_text_to_output("\nStates with high isolation score:\n")

write_text_to_output("|State| |Isolation Score|")

write_table_to_output(highly_isolated_states_data)

## No/single assignment ##
write_text_to_output("\nStates with no/single assignment:\n")

write_table_to_output(single_assigned_states)

write_table_to_output(unassigned_states)

## Redundant states ##
write("\nDetermined redundant states:", file = output_file, append = TRUE)
if (length(redundant_states) == 0) {
  write("NONE", file = output_file, append = TRUE)
} else {
  write(redundant_states, file = output_file, append = TRUE)
}
