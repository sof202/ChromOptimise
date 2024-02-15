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
## PREREQUISITES:                                                 ||
## Run 5_batch_CreateIncrementalModels.sh                         ||
## ============================================================== ##
## INPUTS:                                                        ||
## $1 -> Location of configuation file                            ||
## $2 -> Location of transition matrix file                       ||
## $3 -> Number of states in model                                ||
## $4 -> Directory to place output files into                     ||
## ============================================================== ##
## OUTPUTS:                                                       ||
## Text file containing the Euclidean distance between each       ||
## possible pair of states in the selected model                  ||
## ============================================================== ##

## ========== ##
##   SET UP   ##
## ========== ##

rm(list = ls())

arguments <- commandArgs(trailingOnly = TRUE)
config_file_location <- arguments[1]
emssions_file <- arguments[2]
model_size <- as.numeric(arguments[3])
output_file_path <- arguments[4]

source(config_file_location)
setwd(model_dir)

## ================= ##
##   LOADING FILES   ##
## ================= ##

# This is required as ChromHMM writes the file with the state numbers in the
# first row and column (which is not wanted for this analysis)
emissions_data <- read.table(emissions_file, skip = 1)
emissions_data <- subset(emissions_data, select = -V1)

## ================================ ##
##   EUCLIDEAN DISTANCE FUNCTIONS   ##
## ================================ ##

calculate_euclidean_distances <- function(emissions_data) {
  state_pair_distances <- data.frame(reference_state = numeric(),
                                     comparison_state = numeric(),
                                     euclidean_distance = double())

  for (reference_state_index in 1:(model_size - 1)){
    for (comparison_state_index in (reference_state_index + 1):model_size){
      reference_state <- as.numeric(emissions_data[reference_state_index, ])
      comparison_state <- as.numeric(emissions_data[comparison_state_index, ])

      euclidean_distance <- sqrt(sum((reference_state - comparison_state) ^ 2))

      state_pair_distances[nrow(state_pair_distances) + 1, ] <-
        c(reference_state, comparison_state, euclidean_distance)
    }
  }

  return(state_pair_distances)
}
## ======== ##
##   MAIN   ##
## ======== ##

euclidean_distances_table <- calculate_euclidean_distances(emissions_data)

## =========== ##
##   OUTPUTS   ##
## =========== ##

setwd(output_file_path)

output_file_name <- paste0("Euclidean_distances_model-",model_size,".txt")

write.table(flanking_states_table,
            output_file_name,
            row.names = FALSE)
