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
## OUTPUTS:                                                       ||
## Line plot of estimated log likelihoods against the number      ||
## of states used.                                                ||
## ============================================================== ##


## ========== ##
##   SET UP   ##
## ========== ##

arguments <- commandArgs(trailingOnly = TRUE)
renv_environment <- arguments[1]
likelihoods_file <- arguments[2]
output_file_path <- arguments[3]

renv::load(renv_environment)
library(ggplot2)

## ============== ##
##   PROCESSING   ##
## ============== ##

likelihood_data <- data.table::fread(likelihoods_file)

likelihood_data <- dplyr::select(likelihood_data, c(V5, V7))
names(likelihood_data) <- c("number_of_states", "estimated_log_likelihood")


## ============== ##
##    PLOTTING    ##
## ============== ##

likelihood_plot <-
  ggplot(
    likelihood_data,
    aes(
      x = number_of_states,
      y = estimated_log_likelihood
    )
  ) +
  geom_smooth(formula = y ~ log(x), se = FALSE, span = 1.2) +
  geom_point() +
  labs(x = "Number of States", y = "Estimated Log likelihood") +
  theme_bw()

options(bitmapType = "cairo")
ggsave(
  file.path(
    output_file_path,
    "LikelihoodPlot.png"
  ),
  likelihood_plot
)
