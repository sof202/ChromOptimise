## ============================================================== ##
##                                                                ||
##                            PREAMBLE                            ||
##                                                                ||
## ============================================================== ##
## PURPOSE:                                                       ||
## This R script is in place to find and plot the cosine          ||
## similarity scores and the Euclidean distances between the      ||
## emission parameters for each pair of states that are           ||
## in a hidden Markov model that has been produced by ChromHMM.   ||
## These scores are then plotted in a histogram so that an        ||
## appropriate threshold can be chosen for when two states        ||
## are 'sufficiently distinct' from one another.                  ||
##                                                                ||
## In addition to this, a suggested threshold is given for the    ||
## Euclidean distances histogram. This value is determined by     ||
## looking for gaps in the histogram.                             ||
## ============================================================== ##
## AUTHOR: Sam Fletcher s.o.fletcher@exeter.ac.uk                 ||
## CREATED: November 2023                                         ||
## ============================================================== ##
## PREREQUISITES:                                                 ||
## Run: Generate_Big_Model.sh                                     ||
## ============================================================== ##
## INPUTS:                                                        ||
## $1 -> Model size                                               ||
## ============================================================== ##
## OUTPUTS:                                                       ||
## 2 histograms displaying the distribution of Euclidean          ||
## distances and the cosine similarity scores for the emission    ||
## matrix that is loaded.                                         ||
## ============================================================== ##

## ========== ##
##    SETUP   ##
## ========== ##

rm(list = ls())

setwd("/lustre/projects/Research_Project-MRC190311/scripts")
source("integrative/blueprint/config/config.R")
setwd(big_models_dir)

if (!require("pracma", quietly = TRUE))
  install.packages("pracma")
library("pracma")
library("ggplot2")
library("stringr")

model_size <- commandArgs(trailingOnly = TRUE)
model_size <- as.numeric(model_size)
bin_size <- 0.025

file_name <- paste0("emissions_", model_size, "_1.txt")
emission_data <- read.table(file_name, skip = 1)
emission_data <- subset(emission_data, select = -V1)

euclidean_distance <- function(vector_a, vector_b) {
  sqrt(sum((vector_a - vector_b) ^ 2))
}


## ================================== ##
##   SIMILARITY METRIC CALCULATIONS   ##
## ================================== ##

# Calculate the cosine similarity score and Euclidean distance for each
# disjoint pair of vectors in the emissions text file
cosine_similarity_scores <-
  data.frame(Score = double(), stringsAsFactors = FALSE)
euclidean_distance_scores <-
  data.frame(Score = double(), stringsAsFactors = FALSE)

for (reference_state_index in 1:(model_size - 1)){
  for (comparison_state_index in (reference_state_index + 1):model_size){
    reference_state <- as.numeric(emission_data[reference_state_index, ])
    comparison_state <- as.numeric(emission_data[comparison_state_index, ])

    dot_product_of_states <- dot(reference_state, comparison_state)
    reference_state_magnitude <- norm(reference_state, type = "2")
    comparison_state_magnitude <- norm(comparison_state, type = "2")

    cosine_similarity_score <-
      dot_product_of_states /
      (reference_state_magnitude * comparison_state_magnitude)
    distance_between_vectors <-
      euclidean_distance(reference_state, comparison_state)

    cosine_similarity_scores[nrow(cosine_similarity_scores) + 1, ] <-
      cosine_similarity_score
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

cosine_similarity_histogram <-
  ggplot(cosine_similarity_scores, aes(x = Score)) +
  theme_minimal() +
  geom_histogram(binwidth = bin_size, color = "black", fill = "white") +
  labs(title = "Histogram of cosine similarity Scores",
       x = "Cosine Similarity", y = "Frequency") +
  theme(plot.title = element_text(hjust = 0.5))


## ============== ##
##   SAVE PLOTS   ##
## ============== ##

euclidean_distance_plot_name <-
  paste0("euclidean_distance_histogram.model_size.", model_size, ".pdf")
cosin_similarity_plot_name <-
  paste0("cosine_similarity_histogram.model_size.", model_size, ".pdf")

ggsave(
  euclidean_distance_plot_name,
  plot = euclidean_distance_histogram,
  path = emissions_plotting_dir
)

ggsave(
  cosin_similarity_plot_name,
  plot = cosine_similarity_histogram,
  path = emissions_plotting_dir
)
