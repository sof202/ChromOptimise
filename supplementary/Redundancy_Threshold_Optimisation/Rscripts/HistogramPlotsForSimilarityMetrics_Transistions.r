## ============================================================== ##
##                                                                ||
##                            PREAMBLE                            ||
##                                                                ||
## ============================================================== ##
## PURPOSE:                                                       ||
## This R script is in place to find and plot the maximum         ||
## transistion probability seen in each column of the             ||
## transition matrix produced as an output from ChromHMM.         ||
## These can then be used to identify the relative assignment of  ||
## each state. If the maximum transition probability is very low  ||
## this implies that the state is unlikely to be assigned to      ||
## genomic bins.                                                  ||
## ============================================================== ##
## AUTHOR: Sam Fletcher s.o.fletcher@exeter.ac.uk                 ||
## CREATED: December 2023                                         ||
## ============================================================== ##
## PREREQUISITES:                                                 ||
## Run: Generate_Big_Model.sh                                     ||
## ============================================================== ##
## INPUTS:                                                        ||
## $1 -> Model size                                               ||
## $2 -> Random seed                                              ||
## $3 -> Path to directory containing the model files             ||
## ============================================================== ##
## OUTPUTS:                                                       ||
## A scatter plot displaying the maxima for the columns seen in   ||
## the transition parameter matrix that is loaded.                ||
## ============================================================== ##

## ========== ##
##    SETUP   ##
## ========== ##

rm(list = ls())

library("ggplot2")

setwd("/lustre/projects/Research_Project-MRC190311/scripts")
source("integrative/blueprint/config/config.R")

arguments <- commandArgs(trailingOnly = TRUE)
model_size <- as.numeric(arguments[1])
seed <- arguments[2]
input_path <- arguments[3]

setwd(input_path)

file_name <- paste0("transitions_", model_size, "_", seed, ".txt")
transition_data <- read.table(file_name, skip = 1)
transition_data <- subset(transition_data, select = -V1)

# Rename the columns to be the state numbers
for (column in 1:model_size){
  names(transition_data)[names(transition_data) == paste0("V", column + 1)] <-
    column
}


## ================ ##
##   CALCULATIONS   ##
## ================ ##

transition_state_to_maxima <- data.frame(State = integer(), Maxima = double())
transition_state_from_maxima <- data.frame(State = integer(), Maxima = double())
baseline_probability <- 1 / model_size

# Find maximum transistion probability TO each state
for (state in 1:model_size){
  maximum_transition_probability <- max(transition_data[, state])
  transition_state_to_maxima[nrow(transition_state_to_maxima) + 1, ] <-
    c(state, maximum_transition_probability)
}


## ============ ##
##   PLOTTING   ##
## ============ ##

transition_to_scatter <-
  ggplot(transition_state_to_maxima, aes(x = State, y = Maxima)) +
  geom_point() +
  geom_hline(yintercept = baseline_probability, color = "red") +
  scale_y_continuous(limits = c(0, NA),
                     breaks = sort(c(seq(0, 1, by = 0.1),
                                     baseline_probability))) +
  scale_x_continuous(breaks = c(seq(0, model_size))) +
  annotate("text", x = model_size - 2, y = baseline_probability / 2,
           label = "Random chance", color = "red") +
  labs(x = "Arriving state of the model",
       title = "Maximum probabilities of transitioning towards states") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

## ================ ##
##   SAVING PLOTS   ##
## ================ ##

transition_to_scatter_name <-
  paste0("Scatter.Plot.Maximum.Probability.Transition.Towards.model_size",
         model_size, ".pdf")

ggsave(
  transition_to_scatter_name,
  plot = transition_to_scatter,
  path = transition_plotting_dir,
  height = 8.3,
  width = 23.4,
  units = "in"
)
