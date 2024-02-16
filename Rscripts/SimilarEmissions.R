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
## $3 -> Directory to place output files into                     ||
## $4 -> Boolean for if plots are to be made                      ||
## ============================================================== ##
## OUTPUTS:                                                       ||
## Text file containing the Euclidean distance between each       ||
## possible pair of states in the selected model                  ||
## ============================================================== ##

## ========== ##
##   SET UP   ##
## ========== ##

rm(list = ls())

library(ggplot2)

arguments <- commandArgs(trailingOnly = TRUE)
config_file_location <- arguments[1]
emissions_file <- arguments[2]
output_file_path <- arguments[3]
plotting_flag <- arguments[4]

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

calculate_euclidean_distances <- function(emissions_data, model_size) {
  state_pair_distances <- data.frame(reference_state = numeric(),
                                     comparison_state = numeric(),
                                     euclidean_distance = double())

  for (reference_state_index in 1:(model_size - 1)){
    for (comparison_state_index in (reference_state_index + 1):model_size){
      reference_state <- as.numeric(emissions_data[reference_state_index, ])
      comparison_state <- as.numeric(emissions_data[comparison_state_index, ])

      euclidean_distance <- sqrt(sum((reference_state - comparison_state) ^ 2))

      state_pair_distances[nrow(state_pair_distances) + 1, ] <-
        c(reference_state_index, comparison_state_index, euclidean_distance)
    }
  }

  return(state_pair_distances)
}

## ======== ##
##   MAIN   ##
## ======== ##

model_size <- nrow(emissions_data)

euclidean_distances_table <-
  calculate_euclidean_distances(emissions_data, model_size)

## ========== ##
##   OUTPUT   ##
## ========== ##

setwd(output_file_path)

output_file_name <- paste0("Euclidean_distances_model-", model_size, ".txt")

write.table(flanking_states_table,
            output_file_name,
            row.names = FALSE)

## ============ ##
##   PLOTTING   ##
## ============ ##

# Due to the heatmap, we now need the euclidean distances for every pair
# of states, even for the repeats (distance between 1 and 2, distance between
# 2 and 1).
distances_for_plotting <- function(emissions_data) {
  plotting_data_frame <- data.frame(reference_state = numeric(),
                                    comparison_state = numeric(),
                                    euclidean_distance = numeric())
  
  for (reference_state_index in 1:model_size) {
    for (comparison_state_index in 1:model_size) {
      reference_state <- emissions_data[, reference_state_index]
      comparison_state <- emissions_data[, comparison_state_index]
      
      euclidean_distance <- sqrt(sum((reference_state - comparison_state) ^ 2))
      
      plotting_data_frame[nrow(plotting_data_frame)+1, ] <- 
        c(reference_state_index, comparison_state_index, euclidean_distance)
    }
  }
  return(plotting_data_frame)
}

create_heatmap <- function(emissions_data) {
  eucldiean_distances <- distances_for_plotting(emissions_data)
  
  euclidean_distances_heatmap <- 
    ggplot(eucldiean_distances, aes(reference_state,
                                    comparison_state,
                                    fill = euclidean_distance)) +
    geom_tile() +
    scale_fill_gradient(low = "blue", high = "white")

  return(euclidean_distances_heatmap)
}


create_histogram <- function(emissions_data) {
  eucldiean_distances <- calculate_euclidean_distances(emissions_data)
  
  euclidean_distance_histogram <- 
    ggplot(eucldiean_distances, aes(x = euclidean_distance)) +
    theme_minimal() +
    geom_histogram(binwidth = 0.05, color = "black", fill = "white") +
    labs(title = "Histogram of Euclidean distances",
         x = "Euclidean Distance", y = "Frequency") +
    theme(plot.title = element_text(hjust = 0.5))
  
  return(euclidean_distance_histogram)
}



if (plotting_flag) {
  euclidean_distances_heatmap <- create_heatmap(emissions_data)
  eucldiean_distances_histogram <- create_histogram(emissions_data)

  heatmap_plot_name <-
    paste0("Euclidean_distances_heatmap_model-", model_size, ".png")
  histogram_plot_name <-
    paste0("Euclidean_distances_histogram_model-", model_size, ".png")

  ggsave(
    heatmap_plot_name,
    plot = euclidean_distances_heatmap
  )
  ggsave(
    histogram_plot_name,
    plot = eucldiean_distances_histogram
  )
}
