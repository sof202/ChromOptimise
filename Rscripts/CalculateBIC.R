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
## $1 -> Bin size                                                 ||
## $2 -> Sample size                                              ||
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

rm(list = ls())

setwd("/lustre/projects/Research_Project-MRC190311/scripts/integrative")
source("ChromOptimise/configuration/config.R")
setwd(likelihood_dir)

library(ggplot2)

arguments <- commandArgs(trailingOnly = TRUE)
bin_size <- arguments[1]
sample_size <- arguments[2]
number_of_observations <- as.numeric(arguments[3])
output_file_path <- arguments[4]

## =============== ##
##   IMPORT DATA   ##
## =============== ##

file_name <- paste0(
  "likelihood.BinSize.", bin_size, ".SampleSize.", sample_size, ".txt"
)
likelihood_data <- read.table(file_name)

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

likelihood_data$parameters <-
  (likelihood_data$number_of_states * number_of_marks) +
  (likelihood_data$number_of_states ^ 2)

likelihood_data$bic <- (log_observations * likelihood_data$parameters) -
  (2 * likelihood_data$estimated_log_likelihood)

# Relative BIC (to the maximum BIC)
# Lower BIC is usally better, note that this is a heuristic not a metric
min_bic <- min(likelihood_data$bic)
likelihood_data$relative_bic <- likelihood_data$bic / min_bic

## =========== ##
##   OUTPUTS   ##
## =========== ##

relative_bic_scatter <-
  ggplot(likelihood_data, aes(number_of_states, relative_bic))

relative_bic_scatter +
  geom_point() +
  scale_x_continuous(breaks = seq(min(likelihood_data$number_of_states),
                                  max(likelihood_data$number_of_states),
                                  by = 1)) +
  labs(title = "Bayesian Information Critereon",
       x = "Number of States", y = "BIC (relative to minimum BIC)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))


setwd(output_file_path)
ggsave(
  "Bayesian_Information_Criterion.png"
)
write.csv(likelihood_data, paste0(output_file_path,
                                  "/Bayesian_Information_Criterion.csv"))
