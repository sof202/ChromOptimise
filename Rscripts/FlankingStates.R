## ============================================================== ##
##                                                                ||
##                            PREAMBLE                            ||
##                                                                ||
## ============================================================== ##
## PURPOSE:                                                       ||
## This calculates the most probable states to flank each         ||
## state in a hidden markov model that is not itself. For         ||
## example: suppose a bin is assigned state 1, then we are        ||
## looking at the most likely state assignment for the next bin   ||
## (up and downstream) of this bin that is not assigned 1.        ||
## ============================================================== ##
## AUTHOR: Sam Fletcher                                           ||
## CONTACT: s.o.fletcher@exeter.ac.uk                             ||
## CREATED: February 2022                                         ||
## ============================================================== ##
## OUTPUTS:                                                       ||
## Text file containing most probable flanking states for each    ||
## state in the model.                                            ||
## ============================================================== ##

## ========== ##
##   SET UP   ##
## ========== ##

arguments <- commandArgs(trailingOnly = TRUE)
renv_environment <- arguments[1]
transitions_file <- arguments[2]
output_file_path <- arguments[3]

renv::load(renv_environment)

## ================= ##
##   LOADING FILES   ##
## ================= ##

# Instead of looking at state assignment files, it is easier to just look
# at the transition probability matrix to inspect the most likely upstream
# or downstream flanking state.
transitions_data <- read.table(transitions_file, skip = 1)

# This is required as ChromHMM writes the file with the state numbers in the
# first row and column (which is not wanted for this analysis)
transitions_data <- subset(transitions_data, select = -V1)

## ====================== ##
##   FLANKING FUNCTIONS   ##
## ====================== ##


upstream_flank <- function(transitions_data, state) {
  # Columns give the transition probability of travelling from some state to
  # the selected state
  # Thus they give the probability of each upstream flank
  state_column <- transitions_data[, state]

  # We don't want to factor in the probability of the upstream bin being
  # the same as the selected state, so we remove this value.
  state_column[state] <- NA

  # The state that has the maximum probability of transitioning towards
  # the selected state is our most likely upstream flank
  return(which.max(state_column))
}

downstream_flank <- function(transitions_data, state) {
  # Rows give the transition probability of travelling from our selected state
  # to some other state
  state_row <- transitions_data[state, ]

  state_row[state] <- NA

  # The state that has the maximum probability of being transitioned to
  # from the selected state is our most likely downstream flank
  return(which.max(state_row))
}

## ======== ##
##   MAIN   ##
## ======== ##

model_size <- nrow(transitions_data)


list_of_states <- seq_along(1:model_size)

list_of_upstream_flanks <-
  sapply(
    list_of_states,
    function(state) upstream_flank(transitions_data, state)
  )

list_of_downstream_flanks <-
  sapply(
    list_of_states,
    function(state) downstream_flank(transitions_data, state)
  )

flanking_states_table <-
  data.frame(
    state = list_of_states,
    likeliest_upstream_flank = list_of_upstream_flanks,
    likeliest_downstream_flank = list_of_downstream_flanks
  )

## =========== ##
##   OUTPUTS   ##
## =========== ##

output_file_name <-
  paste0("Likeliest_flanking_states_model-", model_size, ".txt")

write.table(
  flanking_states_table,
  file.path(
    output_file_path,
    output_file_name
  ),
  row.names = FALSE
)
