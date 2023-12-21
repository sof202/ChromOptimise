## ============================================================== ##
##                                                                ||
##                            PREAMBLE                            ||
##                                                                ||
## ============================================================== ##
## PURPOSE:                                                       ||
## This R script is in place to find and plot the Euclidean       ||
## distances between the emission parameters for each pair of     ||
## states in a ChromHMM modelThese scores are then plotted in a   ||
## histogram so that an appropriate threshold can be chosen for   ||
## when two states are 'sufficiently distinct' from one another.  ||
##                                                                ||
## In addition to this, a suggested threshold is given. This      ||
## value is determined by looking for gaps in the histogram.      ||
## ============================================================== ##
## AUTHOR: Sam Fletcher                                           ||
## CONTACT: s.o.fletcher@exeter.ac.uk                             ||
## CREATED: November 2023                                         ||
## ============================================================== ##
## PREREQUISITES:                                                 ||
## Run Generate_Big_Model.sh or ChromHMM's LearnModel Command     ||
## ============================================================== ##
## INPUTS:                                                        ||
## $1 -> Model size                                               ||
## $2 -> Random seed                                              ||
## $3 -> Path to directory containing the model files             ||
## ============================================================== ##
## OUTPUTS:                                                       ||
## histogram displaying the distribution of Euclidean distances   ||
## for the emission matrix that is loaded.                        ||
## ============================================================== ##

## ========== ##
##    SETUP   ##
## ========== ##

rm(list = ls())

setwd("/lustre/projects/Research_Project-MRC190311/scripts")
source("integrative/ChromHMM_OptimumStates/config/config.R")
setwd(big_models_dir)

if (!require("pracma", quietly = TRUE))
  install.packages("pracma")
library("pracma")
library("ggplot2")
library("stringr")

arguments <- commandArgs(trailingOnly = TRUE)
model_size <- as.numeric(arguments[1])
seed <- arguments[2]
input_path <- arguments[3]
bin_size <- 0.025

setwd(input_path)

file_name <- paste0("emissions_", model_size, "_", seed, ".txt")
emission_data <- read.table(file_name, skip = 1)
emission_data <- subset(emission_data, select = -V1)

euclidean_distance <- function(vector_a, vector_b) {
  sqrt(sum((vector_a - vector_b) ^ 2))
}

## ================================== ##
##   SIMILARITY METRIC CALCULATIONS   ##
## ================================== ##

# Calculate the Euclidean distance for each
# disjoint pair of vectors in the emissions text file
euclidean_distance_scores <-
  data.frame(Score = double(), stringsAsFactors = FALSE)

for (reference_state_index in 1:(model_size - 1)){
  for (comparison_state_index in (reference_state_index + 1):model_size){
    reference_state <- as.numeric(emission_data[reference_state_index, ])
    comparison_state <- as.numeric(emission_data[comparison_state_index, ])

    distance_between_vectors <-
      euclidean_distance(reference_state, comparison_state)

    euclidean_distance_scores[nrow(euclidean_distance_scores) + 1, ] <-
      distance_between_vectors
  }
}

## ========================== ##
##    SUGGESTED THRESHOLDS    ##
## ========================== ##

bins_for_distance <-
  cut(euclidean_distance_scores[, 1],
      breaks = seq(min(euclidean_distance_scores[, 1]),
                   max(euclidean_distance_scores[, 1]), by = bin_size))
frequency_counts_for_bins <-
  table(bins_for_distance)
min_count <- min(frequency_counts_for_bins)
# For segment in ggplot
max_count <- max(frequency_counts_for_bins)
frequency_counts_for_bins <- as.matrix(frequency_counts_for_bins)
first_gap_index <- min(which(frequency_counts_for_bins == min_count))
first_gap_interval <- row.names(frequency_counts_for_bins)[first_gap_index]

# Use regular expression to obtain the upper bound of the bin interval
threshold_suggestion <-
  as.numeric(str_extract_all(first_gap_interval, "\\d+\\.\\d+")[[1]][2])
threshold_suggestion_label <-
  paste0("Suggested Threshold Value: ", threshold_suggestion)

## ============ ##
##   Plotting   ##
## ============ ##

euclidean_distance_histogram <- ggplot(euclidean_distance_scores,
                                       aes(x = Score)) +
  theme_minimal() +
  geom_histogram(binwidth = bin_size, color = "black", fill = "white") +
  geom_segment(aes(x = threshold_suggestion,
                   xend = threshold_suggestion, y = 0, yend = max_count),
               linetype = "dotted") +
  labs(title = "Histogram of Euclidean distances",
       x = "Euclidean Distance", y = "Frequency") +
  theme(plot.title = element_text(hjust = 0.5)) +
  annotate("text", x = threshold_suggestion + 12 * bin_size, y = max_count / 5,
           label = threshold_suggestion_label)

## ============== ##
##   SAVE PLOTS   ##
## ============== ##

euclidean_distance_plot_name <-
  paste0("Euclidean_Distances_Histogram.model_size.", model_size, ".pdf")

ggsave(
  euclidean_distance_plot_name,
  plot = euclidean_distance_histogram,
  path = emission_plotting_dir
)
