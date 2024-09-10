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
## PREREQUISITES:                                                 ||
## Run 6_CompareModels.sh                                         ||
## ============================================================== ##
## INPUTS:                                                        ||
## $1 -> Location of configuation file                            ||
## $2 -> Location of likelihoods file                             ||
## $3 -> Total size of binary files                               ||
## $4 -> Directory to place output files into                     ||
## ============================================================== ##
## OUTPUTS:                                                       ||
## Scatter plot of relative BIC against number of states          ||
## Text file containing BIC for each model                        ||
## ============================================================== ##

## ========== ##
##   SET UP   ##
## ========== ##

library(ggplot2)

arguments <- commandArgs(trailingOnly = TRUE)
config_file_location <- arguments[1]
likelihoods_file <- arguments[2]
optimum_states <- as.numeric(arguments[3])
number_of_observations <- as.numeric(arguments[4])
output_file_path <- arguments[5]

source(config_file_location)

## =============== ##
##   IMPORT DATA   ##
## =============== ##

likelihood_data <- read.table(likelihoods_file)

likelihood_data <- subset(likelihood_data, select = c(V5, V7))
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

# number_of_marks is coming from Config.R
likelihood_data$parameters <-
  (likelihood_data$number_of_states * number_of_marks) +
  (likelihood_data$number_of_states^2)

likelihood_data$bic <- (log_observations * likelihood_data$parameters) -
  (2 * likelihood_data$estimated_log_likelihood)

# Relative BIC (to the maximum BIC)
# Lower BIC is usally better, note that this is a heuristic not a metric
min_bic <- min(likelihood_data$bic)
likelihood_data$relative_bic <- likelihood_data$bic / min_bic

## =================== ##
##   WARNING MESSAGE   ##
## =================== ##

# It's not the worst thing, but the user should be warned if a model has a
# larger BIC than the next smallest state

sub_optimum_states <- optimum_states - 1
optimum_states_bic <-
  subset(likelihood_data, number_of_states == optimum_states)[["bic"]]
sub_optimum_states_bic <-
  subset(likelihood_data, number_of_states == sub_optimum_states)[["bic"]]

if (optimum_states_bic > sub_optimum_states_bic) {
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
