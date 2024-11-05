## ============================================================== ##
##                                                                ||
##                            PREAMBLE                            ||
##                                                                ||
## ============================================================== ##
## PURPOSE:                                                       ||
## Calculates the Bayesian Information Critereon (BIC) for        ||
## ChromHMMmodels, plots the relative BIC for each model against  ||
## the number of states.                                          ||
## ============================================================== ##
## AUTHOR: Sam Fletcher                                           ||
## CONTACT: s.o.fletcher@exeter.ac.uk                             ||
## CREATED: January 2024                                          ||
## ============================================================== ##
## OUTPUTS:                                                       ||
## Scatter plot of relative BIC against number of states          ||
## Text file containing BIC for each model                        ||
## ============================================================== ##

## ========== ##
##   SET UP   ##
## ========== ##

arguments <- commandArgs(trailingOnly = TRUE)
renv_environment <- args[1]
config_file_location <- arguments[2]
likelihoods_file <- arguments[3]
optimum_states <- as.numeric(arguments[4])
number_of_observations <- as.numeric(arguments[5])
output_file_path <- arguments[6]

source(config_file_location)
renv::load(renv_environment)
library(ggplot2)

## =============== ##
##   IMPORT DATA   ##
## =============== ##

likelihood_data <- data.table::fread(likelihoods_file)

likelihood_data <- dplyr::select(likelihood_data, c(V5, V7))
names(likelihood_data) <- c("number_of_states", "estimated_log_likelihood")

## ================ ##
##   CALCULATIONS   ##
## ================ ##

# BIC = ln(n)k - 2ln(L)
# k -> Number of estimated parameters in the model
#      Which is size of emission matrix + size of transition matrix
#      Which is (#states*#marks)+(#states)^2
# n -> Number of observations (total number of lines in binary files)
# L -> Likelihood function (output of chromHMM already in ln(L) form)
log_observations <- log(number_of_observations)

likelihood_data$parameters <-
  (likelihood_data$number_of_states * number_of_marks) +
  (likelihood_data$number_of_states^2)

likelihood_data$bic <- (log_observations * likelihood_data$parameters) -
  (2 * likelihood_data$estimated_log_likelihood)

min_bic <- min(likelihood_data$bic)
likelihood_data$relative_bic <- likelihood_data$bic / min_bic

## =================== ##
##   WARNING MESSAGE   ##
## =================== ##

optimum_states_bic <- dplyr::filter(
  likelihood_data,
  number_of_states == optimum_states
)[["bic"]]

sub_optimum_states_bic <- dplyr::filter(
  likelihood_data,
  number_of_states < optimum_states
)[["bic"]]

if (any(optimum_states_bic > sub_optimum_states_bic)) {
  warning(
    "WARNING: A less complex model has a smaller BIC ",
    "than the model with the determined optimum number of states (",
    optimum_states,
    ")"
  )
}

## =========== ##
##   OUTPUTS   ##
## =========== ##

relative_bic_scatter <-
  ggplot(likelihood_data, aes(number_of_states, relative_bic)) +
  geom_point() +
  scale_x_continuous(breaks = seq_along(likelihood_data$number_of_states)) +
  labs(
    title = "Bayesian Information Critereon",
    x = "Number of States", y = "BIC (relative to minimum BIC)"
  ) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

write.csv(
  likelihood_data,
  file.path(
    output_file_path,
    "Bayesian_Information_Criterion.csv"
  )
)

options(bitmapType = "cairo")
ggsave(
  file.path(
    output_file_path,
    "Bayesian_Information_Criterion.png"
  ),
  relative_bic_scatter
)
