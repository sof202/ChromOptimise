## -------------------------------------------------------------- ##
##                                                                ##
##                            PREAMBLE                            ##
##                                                                ##
## -------------------------------------------------------------- ##
##                            PURPOSE:                            ##
##     This R script is in place to find and plot the cosine      ##
##   similarity scores and the Euclidean distances between the    ##
##      emission parameters for each pair of states that are      ##
##  in a hidden Markov model that has been produced by ChromHMM.  ##
##    These scores are then plotted in a histogram so that an     ##
##   appropriate threshold can be chosen for when two states      ##
##         are 'sufficiently distinct' from one another.          ##
##                                                                ##
##  In addition to this, a suggested threshold is given for the   ##
##  Euclidean distances histogram. This value is determined by    ##
##              looking for gaps in the histogram.                ##
## -------------------------------------------------------------- ##
##        AUTHOR: Sam Fletcher s.o.fletcher@exeter.ac.uk          ##
##                     CREATED: November 2023                     ##
## -------------------------------------------------------------- ##
##                         PREREQUISITES:                         ##
##                   Run: Generate_Big_Model.sh                   ##
## -------------------------------------------------------------- ##
##                            INPUTS:                             ##
##                        $1 -> Model size                        ##
## -------------------------------------------------------------- ##
##                            OUTPUTS:                            ##
##      2 histograms displaying the distribution of Euclidean     ##
##   distances and the cosine similarity scores for the emission  ##
##                    matrix that is loaded.                      ##
## -------------------------------------------------------------- ##

## ---------- ##
##    SETUP   ##
## ---------- ##

rm(list=ls())
setwd("/lustre/projects/Research_Project-MRC190311/blueprint/Big_Model_Files")

# Packages
# install.packages("pracma")
library("pracma")
library("ggplot2")
library("stringr")


# Arguments and variables
Model.Size <- commandArgs(trailingOnly = TRUE)
Model.Size <- as.numeric(Model.Size)
Bin.Size <- 0.025

# Open emissions file and extract the emission parameters from this
File.Name <- paste0("emissions_",Model.Size,"_1.txt")
Emission.Data <- read.table(File.Name, skip=1)
Emission.Data <- subset(Emission.Data, select= -V1)

# Define Euclidean distance function
Euclidean.Distance <- function(Vector.a,Vector.b) {
  sqrt(sum((Vector.a-Vector.b)^2))
}


## ---------------------------------- ##
##   SIMILARITY METRIC CALCULATIONS   ##
## ---------------------------------- ##

# Calculate the cosine similarity score and Euclidean distance for each disjoint pair of vectors in the emissions text file
Cosine.Similarity.Scores <- data.frame(Score=double(), stringsAsFactors=FALSE)
Euclidean.Distance.Scores <- data.frame(Score=double(), stringsAsFactors=FALSE)

for (Reference.State.Index in 1:(Model.Size-1)){
  for (Comparison.State.Index in (Reference.State.Index+1):Model.Size){
    Reference.State <- as.numeric(Emission.Data[Reference.State.Index,])
    Comparison.State <- as.numeric(Emission.Data[Comparison.State.Index,])
    
    Dot.Product.Of.States <- dot(Reference.State, Comparison.State)
    Reference.State.Magnitude <- norm(Reference.State, type="2")
    Comparison.State.Magnitude <- norm(Comparison.State, type="2")
    
    Cosine.Similarity.Score <- Dot.Product.Of.States / (Reference.State.Magnitude * Comparison.State.Magnitude)
    Distance.Between.Vectors <- Euclidean.Distance(Reference.State,Comparison.State)
    
    Cosine.Similarity.Scores[nrow(Cosine.Similarity.Scores)+1,] <- Cosine.Similarity.Score
    Euclidean.Distance.Scores[nrow(Euclidean.Distance.Scores)+1,] <- Distance.Between.Vectors
  }
}


## -------------------------- ##
##    SUGGESTED THRESHOLDS    ##
## -------------------------- ##

## Euclidean Distance
# Find smallest bin where the count is 0
Bins.For.Euclidean.Distance <- cut(Euclidean.Distance.Scores[,1], breaks=seq(min(Euclidean.Distance.Scores[,1]),max(Euclidean.Distance.Scores[,1]),by=Bin.Size))
Frequency.Counts.For.Bins.Euclidean.Distance <- table(Bins.For.Euclidean.Distance)
Min.Count.Euclidean.Distance <- min(Frequency.Counts.For.Bins.Euclidean.Distance)
Max.Count.Euclidean.Distance <- max(Frequency.Counts.For.Bins.Euclidean.Distance) #Purely for segment in ggplot
Frequency.Counts.For.Bins.Euclidean.Distance <- as.matrix(Frequency.Counts.For.Bins.Euclidean.Distance)
First.Gap.In.Euclidean.Distances.Index <- min(which(Frequency.Counts.For.Bins.Euclidean.Distance == Min.Count.Euclidean.Distance))
First.Gap.Interval.In.Euclidean.Distances <- row.names(Frequency.Counts.For.Bins.Euclidean.Distance)[First.Gap.In.Euclidean.Distances.Index]

# Use regular expression to obtain the upper bound of the bin interval
Euclidean.Distances.Threshold.Suggestion <- as.numeric(str_extract_all(First.Gap.Interval.In.Euclidean.Distances,"\\d+\\.\\d+")[[1]][2])
Euclidean.Distances.Threshold.Suggestion.Label <- paste0("Suggested Threshold Value: ", Euclidean.Distances.Threshold.Suggestion)

## ------------ ##
##   Plotting   ##
## ------------ ##  

# Euclidean Distance
Euclidean.Distance.Histogram <- ggplot(Euclidean.Distance.Scores, aes(x=Score)) + 
  theme_minimal() +
  geom_histogram(binwidth=Bin.Size, color="black", fill="white") +
  geom_segment(aes(x=Euclidean.Distances.Threshold.Suggestion, xend=Euclidean.Distances.Threshold.Suggestion, y=0, yend=Max.Count.Euclidean.Distance), linetype="dotted") +
  labs(title="Histogram of Euclidean distances", x="Euclidean Distance", y="Frequency") +
  theme(plot.title=element_text(hjust=0.5)) +
  annotate("text",x=Euclidean.Distances.Threshold.Suggestion+12*Bin.Size, y=Max.Count.Euclidean.Distance/5, label=Euclidean.Distances.Threshold.Suggestion.Label)


# Cosine Similarity
Cosine.Similarity.Histogram <- ggplot(Cosine.Similarity.Scores, aes(x=Score)) + 
  theme_minimal() +
  geom_histogram(binwidth=Bin.Size, color="black", fill="white") +
  labs(title="Histogram of cosine similarity Scores", x="Cosine Similarity", y="Frequency") +
  theme(plot.title=element_text(hjust=0.5))


## -------------- ##
##   SAVE PLOTS   ##
## -------------- ##

Euclidean.Distance.Histogram.Plot.Name <- paste0("Euclidean.Distance.Histogram.Model.Size.",Model.Size,".pdf")
Cosine.Similarity.Histogram.Plot.Name <- paste0("Cosine.Similarity.Histogram.Model.Size.",Model.Size,".pdf")

ggsave(
  Euclidean.Distance.Histogram.Plot.Name,
  plot=Euclidean.Distance.Histogram,
  path="/lustre/projects/Research_Project-MRC190311/blueprint/Plots/Emission_Parameter_Big_Model_Plots"
)

ggsave(
  Cosine.Similarity.Histogram.Plot.Name,
  plot=Cosine.Similarity.Histogram,
  path="/lustre/projects/Research_Project-MRC190311/blueprint/Plots/Emission_Parameter_Big_Model_Plots"
)
