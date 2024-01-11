## ============================================================== ##
##                                                                ||
##                            PREAMBLE                            ||
##                                                                ||
## ============================================================== ##
## PURPOSE:                                                       ||
## Calculates the Akaike Information Critereon (AIC) for ChromHMM ||
## models, plots the relative AIC for each model against the      ||
## number of states.                                              ||
## ============================================================== ##
## AUTHOR: Jessica Shields                                        ||
## CONTACT: s.o.fletcher@exeter.ac.uk                             ||
## CREATED: October 2022                                          ||
## MODIFIED: Sam Fletcher (January 2024)                          ||
## ============================================================== ##
## PREREQUISITES:                                                 ||
## Run 5_batch_CreateIncrementalModels.sh                         ||
## ============================================================== ##
## INPUTS:                                                        ||
## $1 -> Bin size                                                 ||
## $2 -> Sample size                                              ||
## $3 -> Directory to place output files into                     ||
## ============================================================== ##
## OUTPUTS:                                                       ||
## Scatter plot of relative AIC against number of states          ||
## Text file containing AIC and relative AIC for each model       ||
## ============================================================== ##

## ========== ##
##   SET UP   ##
## ========== ##

rm(list = ls())

setwd("/lustre/projects/Research_Project-MRC190311/scripts/integrative")
source("ChromOptimise/configuration/config.R")
setwd(likelihood_dir)

library(ggplot2)
library(magrittr)

arguments <- commandArgs(trailingOnly = TRUE)
bin_size <- arguments[1]
sample_size <- arguments[2]
output_file_path <- arguments[3]

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

# AIC = 2k - 2ln(L)
# k -> Number of estimated parameters in the model
#      Which is size of emission matrix + size of transition matrix
#      Which is (#states*#marks)+(#states)^2
# L -> Likelihood function (output of chromHMM already in ln(L) form)
likelihood_data$parameters <- (likelihood_data$number_of_states*number_of_marks)+(likelihood_data$number_of_states^2) 
likelihood_data$aic <- (2*likelihood_data$parameters)-(2*likelihood_data$estimated_log_likelihood)

# Relative AIC
# These values are proportional to the probability that each model minimises
# the estimated information loss

min_aic <- min(likelihood_data$aic)
likelihood_data$relative_aic <- exp((min_aic - likelihood_data$aic) / 2)

## =========== ##
##   OUTPUTS   ##
## =========== ##

relative_aic_scatter <- 
  ggplot(likelihood_data, aes(number_of_states, relative_aic)) 

relative_aic_scatter +
  geom_point() +
  theme_bw()

setwd(output_file_path)
ggsave(
  "Akaike_Information_Criterion.png"
)
write.csv(likelihood_data, paste0(output_file_path,"/Akaike_Information_Criterion.csv"))

