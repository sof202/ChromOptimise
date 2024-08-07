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
## $1 -> Location of transition matrix file                       ||
## $2 -> Directory to place output files into                     ||
## $3 -> Boolean for if plots are to be made                      ||
## ============================================================== ##
## OUTPUTS:                                                       ||
## Text file containing the Euclidean distance between each       ||
## possible pair of states in the selected model                  ||
## ============================================================== ##

## ========== ##
##   SET UP   ##
## ========== ##

rm(list = ls())

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}

library(ggplot2)

arguments <- commandArgs(trailingOnly = TRUE)
emissions_file <- arguments[1]
output_file_path <- arguments[2]
plotting_flag <- arguments[3]

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
  state_pair_distances <- data.frame(
    reference_state = numeric(),
    comparison_state = numeric(),
    euclidean_distance = double()
  )

  for (reference_state_index in 1:(model_size - 1)) {
    for (comparison_state_index in (reference_state_index + 1):model_size) {
      reference_state <- as.numeric(emissions_data[reference_state_index, ])
      comparison_state <- as.numeric(emissions_data[comparison_state_index, ])

      euclidean_distance <- sqrt(sum((reference_state - comparison_state)^2))

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

write.table(euclidean_distances_table,
  output_file_name,
  row.names = FALSE
)

## ============ ##
##   PLOTTING   ##
## ============ ##

create_heatmap <- function(emissions_data) {
  eucldiean_distances <-
    calculate_euclidean_distances(emissions_data, model_size)

  euclidean_distances_heatmap <-
    ggplot(eucldiean_distances, aes(comparison_state,
      reference_state,
      fill = euclidean_distance
    )) +
    geom_tile() +
    scale_fill_gradient(low = "blue", high = "white") +
    labs(title = "Heatmap of Euclidean distances between state pairs") +
    theme(plot.title = element_text(hjust = 0.5))

  return(euclidean_distances_heatmap)
}


create_histogram <- function(emissions_data) {
  eucldiean_distances <-
    calculate_euclidean_distances(emissions_data, model_size)

  euclidean_distance_histogram <-
    ggplot(eucldiean_distances, aes(x = euclidean_distance)) +
    geom_histogram(binwidth = 0.05, color = "black", fill = "white") +
    labs(
      title = "Histogram of Euclidean distances",
      x = "Euclidean Distance", y = "Frequency"
    ) +
    theme(plot.title = element_text(hjust = 0.5))

  return(euclidean_distance_histogram)
}

if (!exists("plotting_flag")) {
  plotting_flag <- FALSE
}

if (plotting_flag) {
  euclidean_distances_heatmap <- create_heatmap(emissions_data)
  eucldiean_distances_histogram <- create_histogram(emissions_data)

  heatmap_plot_name <-
    paste0("Euclidean_distances_heatmap_model-", model_size, ".png")
  histogram_plot_name <-
    paste0("Euclidean_distances_histogram_model-", model_size, ".png")

  options(bitmapType = "cairo")
  ggsave(
    heatmap_plot_name,
    plot = euclidean_distances_heatmap
  )
  ggsave(
    histogram_plot_name,
    plot = eucldiean_distances_histogram
  )
}
