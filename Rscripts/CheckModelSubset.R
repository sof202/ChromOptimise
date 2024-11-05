## ============================================================== ##
##                                                                ||
##                            PREAMBLE                            ||
##                                                                ||
## ============================================================== ##
## PURPOSE:                                                       ||
## This calculates the Euclidean distance between each pair of    ||
## states in a HMM.                                               ||
## ============================================================== ##
## AUTHOR: Sam Fletcher                                           ||
## CONTACT: s.o.fletcher@exeter.ac.uk                             ||
## CREATED: February 2022                                         ||
## ============================================================== ##
## OUTPUTS:                                                       ||
## Text file containing the Euclidean distance between each       ||
## possible pair of states in the selected model                  ||
## ============================================================== ##

## ========== ##
##   SET UP   ##
## ========== ##

arguments <- commandArgs(trailingOnly = TRUE)
renv_environment <- args[1]
config_file_location <- arguments[2]
smaller_file <- arguments[3]
bigger_file <- arguments[4]
threshold_type <- arguments[5]

source(config_file_location)
renv::load(renv_environment)

## ================= ##
##   LOADING FILES   ##
## ================= ##

# Overlap data has the genome% column and the Base row which we want to remove
# Emissions data can have emissions in different orders between models
if (threshold_type == "overlap") {
  process_data <- function(file) {
    data <- data.table::fread(file, skip = 1)
    data <- data |>
      dplyr::select(-c(V1, V2)) |>
      dplyr::filter(!dplyr::row_number() == nrow(data))

    return(data)
  }
  smaller_data <- process_data(smaller_file)
  bigger_data <- process_data(bigger_file)
} else if (threshold_type == "emissions") {
  process_data <- function(file) {
    data <- data.table::fread(file)
    colnames(data)[1] <- "state"
    data <- data |>
      dplyr::select(order(colnames(data))) |>
      dplyr::select(-state)

    return(data)
  }
  smaller_data <- process_data(smaller_file)
  bigger_data <- process_data(bigger_file)
}

## ================ ##
##   CALCULATIONS   ##
## ================ ##


calc_min_euclidean_distances <- function(smaller_data, bigger_data) {
  euclidean_distance <- function(a, b) {
    return(sqrt(sum((a - b)^2)))
  }
  distances_matrix <- outer(
    seq_len(nrow(smaller_data)),
    seq_len(nrow(bigger_data)),
    Vectorize(function(i, j) {
      euclidean_distance(smaller_data[i, ], bigger_data[j, ])
    })
  )

  min_euclidean_distances <- as.list(apply(distances_matrix, 1, min))

  return(min_euclidean_distances)
}

calc_min_euclidean_distances(smaller_data, bigger_data)

## ======== ##
##   MAIN   ##
## ======== ##

threshold_to_use <- function(threshold_type) {
  is_valid_threshold <- function(threshold) {
    if (is.na(threshold)) {
      error_message <- paste(threshold, "is not set in Config.R.")
      stop(error_message)
    }
  }
  if (threshold_type == "emissions") {
    is_valid_threshold(emissions_threshold)
    return(emissions_threshold)
  } else if (threshold_type == "overlap") {
    is_valid_threshold(overlap_threshold)
    return(overlap_threshold)
  } else {
    stop("Incorrect threshold_type given (must be emissions or overlap)")
  }
}

are_states_a_subset <- function(smaller_data,
                                bigger_data,
                                threshold) {
  min_euclidean_distances <- calc_min_euclidean_distances(
    smaller_data,
    bigger_data
  )

  states_in_both_models <- min_euclidean_distances < threshold
  return(all(states_in_both_models))
}


## ========== ##
##   OUTPUT   ##
## ========== ##

threshold <- threshold_to_use(threshold_type)
cat(are_states_a_subset(smaller_data, bigger_data, threshold))
