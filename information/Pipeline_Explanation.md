# Pipeline Explanation
For a visual representation of the pipeline, please consult the file: `Optimal_States_Analysis_Pipeline.pdf` (found in the same directory as this README file).
\
The process is split into 3 main steps:
- Manipulate the data set using user defined parameters
- Use ChromHMM's pipeline to obtain multiple hidden Markov models on the processed data
- Use statistical analysis to evaluate the model with the optimum number of states

The pipeline expects aligned `.bam` files as its input. The pipeline will then output:
- A scatter plot that showcases how the number of states in the hidden Markov model influences the estimated log likelihood of the model
- A file containing the optimum number of states to use with the binarized data files (when a specific bin size has been used)

Note 1: This pipeline was built with the blueprint data from EGA in mind. Many steps are in place purely to work with these files specifically. Scripts 3 and onwards will work with any .bam files from a data set that can be peak called.
\
Note 2: This pipeline was built using the SLURM workload manager. Scripts will contain errors (scripts of the form `x_batch_xxx.sh` will fail) if not ran using the `sbatch` SLURM command. Please adapt the scripts accordingly if SLURM is not being used.

## Main pipeline shell scripts
### `0_DownloadBluePrint.sh` 
- This is an artefact of the pipeline originally being built for the blueprint data from EGA.
    - The script will download files from a file of file names using the pyega3 python package
        - Login credentials need to be specified in order for the download to be authorised. This will be a file in `.json` format that includes your EGA username and password.
### `1_MoveFilesToSingleDirectory.sh` 
- Again, this is an artefact of the pipeline originally being built for the blueprint data from EGA.
    - The downloaded files from EGA will be in separate directories (as each directory has its own checksum file to validate the downloads). For easier file maniupulation later in the pipeline, these files are collated into folders denoted by the epigenetic mark they contain information on.
### `2_batch_ProcessBamFiles.sh` 
- If the aligned `.bam` files you are working with have already been processed, skip this step.
    - The purpose of this step is to sort, index, filter and remove the duplicates of the aligned `.bam` files.
        - Note that, in the blueprint data, the aligned read files are single end and duplicates have already been marked for duplicate reads using picard.
    - This step also outputs some statistics on the aligned read files using `samtools stats` and `samtools idxstats` before and after processing.
### `3_SubsampleBamFiles.sh` 
- This merges all of the processed `.bam` files from the previous step using `samtools merge` and subsequently subsamples this file using `samtools view -s`
    - This script takes a subsample size as the input from the user.
        - Note that samtools had a known bug in version 0.1.18 where subsampling would fail if the sample size was above 50%, ensure that your version of samtools is up to date.
    - The reason for merging files and subsampling instead of simply subsampling the original files is that it is unlikely that every processed `.bam` file is of the same size. If the `.bam` files are of different sizes, then different subsamples (when subsampling the original files without merging) can have different sizes. 
    - This results in the different samples no longer being directly comparable.
    - Merging the files first avoids this problem.
    - Note that if you are working with a lot of files, this step can use a large amount of memory due to the way samtools merges files.
### `4_BinarizeFiles.sh` 
- This uses ChromHMM's `BinarizeBam` command to binarize the subsampled `.bam` file produced in the previous step. 
    - It is assumed at this point that the above steps have been completed for all of the epigenetic marks that you aim to inspect.
    - This step takes a user input for the bin size, see about for how this will effect your results. If no bin size is given, the default of 200 is used (ChromHMM's default bin size)
 ### `5_batch_CreateIncrementalModels.sh` 
 - This uses ChromHMM's `LearnModel` command to learn multiple hidden Markov models for the inputted data with varying number of states. It also obtains the estimated log likelihood value of each model. 
    - The number of models and the increment to be used between them is user specified.
    - Note that this is using the 'information' initialisation method for the starting parameter set for the model. As such, the number of states in the model cannot exceed the total number of combinations of marks in your dataset.
        - One may not know this number before running the script. However, the maximum number of mark combinations one can have for a binary file will always be capped by 2<sup>k</sup> where k is the total number of marks. If one exceeds the maximum number of states permitted for this initialisation technique, ChromHMM will tell the user the maximum number of states allowed in the associated error message. 
    - This script will always generate a model with 2 states. Inspecting the emission parameters for this simple model is a good way of validating your data.
        - Generally, most of your genomic data will be non-coding regions. The 2 state model will show if this is the case in the overlap files produced.
### `6_OptimalNumberOfStates.sh` 
- This is the statistical analysis step where a plot of estimated log likelihood against number of model states is created and the optimum number of states is determined.
    - This script used two R scripts:
        - PlotLikelihoods.R
            - This creates the above described plot giving the user feedback on how the accuracy of the model increases as one increases the number of states used in the model.
            - The script will also give insights on the estimated maximum log likelihood possible on the data set (if one were to use an abundance of states)
        - OptimumNumberOfStates.R
            - This uses cosine similarity to identify how similar each pair states are in terms of their emission parameters in each model
                - If two states are similar under this metric, then one must be redundant.
            - Redudant states are identified and used to reject higher complexity models


## Additional Shell Scripts
There are a number of additional scripts that are given in JobSubmission/Extra_Scripts that are not directly a part of the main pipeline, but proved useful to constructing the pipeline. Such scripts are outlined here:
### `Generate_Big_Model.sh` 
- This uses the binarized data produced by `4_BinarizeBamFiles.sh` to produce a single very complex model with a high number of user specified states.
    - The model produced will use a random initialisation of emission/transition paramters so that an arbitrary number of states can be used in the model (default initialisation method does not allow for this)
    - These large models were used with blueprint data to obtain sensible threshold values for the similarity metrics used in OptimumNumberOfStates.R
### `CompareModels.sh` 
- This uses ChromHMM's `CompareModels` command to generate comparisons between the models produced in the previous step.
    - Note that `CompareModels` only looks at the emission parameters for each state in the models inputted.
    - This script will compare each model with all of the models that are less complex than it (so a model with 8 states will be compared against only those models that have fewer than 8 states).
### `Generate_Redundancy_Metrics_Plots.sh`
- This runs two R scripts on a specified model (recommended to be a very complex model) so that the user can make informed decisions on the thresholds used in determining redundant states in `6_OptimalNumberOfStates.sh`.
    - HistogramPlotForEuclideanDistances.R -> Creates a histogram for the Euclidean distances between pairs of states in the model provided.
    - ScatterPlotForTransitionMaxima.R -> Creates a scatter plot for the maximum transition probability towards each state in the model provided.
    - The user should input a very complex model as this will result in the two plots possessing two obvious groups, redundant states and useful states.
        - For Euclidean distances -> Redundant state pairs have very low Euclidean distance, other wise they have much higher Euclidean distance (with a noticable gap between these groups).
        - For Transition probabilities -> Redundant state candidates score very low whilst useful state candidates score close to 1 (with a noticable gap between these groups).
- This script is only optional as the user is free to set their own redundancy determining thresholds if they have a good idea for the values from other sources.