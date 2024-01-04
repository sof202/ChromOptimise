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
## Run 6_CompareModels.sh                                         ||
## ============================================================== ##
## INPUTS:                                                        ||
## $1 -> Bin size                                                 ||
## $2 -> Sample size                                              ||
## $3 -> Directory to place output files into                     ||
## ============================================================== ##
## OUTPUTS:                                                       ||
## Line plot of estimated log likelihoods against the number      ||
## of states used.                                                ||
## ============================================================== ##


## ========== ##
##   SET UP   ##
## ========== ##

rm(list = ls())

setwd("/lustre/projects/Research_Project-MRC190311/scripts/integrative")
source("ChromOptimise/configuration/config.R")
setwd(likelihood_dir)

rm(list = ls())

if (!require("dplyr", quietly = TRUE))
  install.packages("dplyr")

library("ggplot2")
library("dplyr")

arguments <- commandArgs(trailingOnly = TRUE)
bin_size <- arguments[1]
sample_size <- arguments[2]
output_file_path <- arguments[3]


## ============== ##
##   PROCESSING   ##
## ============== ##

file_name <- paste0(
  "likelihood.BinSize.", bin_size, ".SampleSize.", sample_size, ".txt"
)
likelihood_data <- read.table(file_name)

likelihood_data <- subset(likelihood_data, select = c(V5, V7))
names(likelihood_data) <- c("Number_Of_States", "Estimated_Log_Likelihood")


## ============== ##
##    PLOTTING    ##
## ============== ##

likelihood_plot <- ggplot(likelihood_data,
                          aes(x = Number_Of_States,
                              y = Estimated_Log_Likelihood))

likelihood_plot +
  geom_smooth(formula = y ~ log(x), color = "blue", se = FALSE, span = 1.2) +
  theme_bw() +
  geom_point(shape = "square", color = "black")

setwd(output_file_path)
ggsave(
  "LikelihoodPlot.png"
)
