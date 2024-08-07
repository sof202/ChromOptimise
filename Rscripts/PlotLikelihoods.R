## ============================================================== ##
##                                                                ||
##                            PREAMBLE                            ||
##                                                                ||
## ============================================================== ##
## PURPOSE:                                                       ||
## Plots the estimated log  likelihood values obtained from       ||
## ChromHMM's LearnModel command against the number of states     ||
## used in the model.                                             ||
## ============================================================== ##
## AUTHOR: Sam Fletcher                                           ||
## CONTACT: s.o.fletcher@exeter.ac.uk                             ||
## CREATED: November 2023                                         ||
## ============================================================== ##
## PREREQUISITES:                                                 ||
## Run 5_batch_CreateIncrementalModels.sh                         ||
## ============================================================== ##
## INPUTS:                                                        ||
## $1 -> Location of likelihoods file                             ||
## $2 -> Directory to place output files into                     ||
## ============================================================== ##
## OUTPUTS:                                                       ||
## Line plot of estimated log likelihoods against the number      ||
## of states used.                                                ||
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
likelihoods_file <- arguments[1]
output_file_path <- arguments[2]


## ============== ##
##   PROCESSING   ##
## ============== ##

likelihood_data <- read.table(likelihoods_file)

likelihood_data <- subset(likelihood_data, select = c(V5, V7))
names(likelihood_data) <- c("number_of_states", "estimated_log_likelihood")


## ============== ##
##    PLOTTING    ##
## ============== ##

likelihood_plot <- ggplot(
  likelihood_data,
  aes(
    x = number_of_states,
    y = estimated_log_likelihood
  )
)

likelihood_plot +
  geom_smooth(formula = y ~ log(x), color = "blue", se = FALSE, span = 1.2) +
  geom_point(shape = "square", color = "black") +
  labs(x = "Number of States", y = "Estimated Log likelihood") +
  theme_bw()

setwd(output_file_path)
ggsave(
  "LikelihoodPlot.png"
)
