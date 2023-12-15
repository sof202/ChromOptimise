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

Note: This pipeline was built with the blueprint data from EGA in mind. Many steps are in place purely to work with these files specifically. Scripts 3 and onwards will work with any data set that can be peak called.
\
Note 2: This pipeline was built using the SLURM Workload Manager. Scripts will likely fail if not ran using the `sbatch` SLURM command. This is because many scripts have paralellisation built into them using SLURM's array feature. Please adapt the scripts accordingly if SLURM is not being used.

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
    - This script takes a subsample size as the input from the user. If no input is recieved, a default of 50 is used.
        - Note that samtools had a known bug in version 0.1.18 where subsampling would fail if the sample size was above 50%, ensure that your version of samtools is up to date.
    - The reason for merging files and subsampling instead of simply subsampling the original files is that it is unlikely that every processed `.bam` file is of the same size. If the `.bam` files are of different sizes, then different subsamples (when subsampling the original files without merging) can have different sizes. 
    - This results in the different samples no longer being directly comparable.
    - Merging the files first avoids this problem.
    - Note that if you are working with a lot of files, this step can use a large amount of memory due to the way samtools merges files.
### `4_BinarizeFiles.sh` 
- This uses ChromHMM's `BinarizeBam` command to binarize the subsampled `.bam` file produced in the previous step. 
    - It is assumed at this point that the above steps have been completed for all of the epigenetic marks that you aim to inspect.
    - This step takes a user input for the bin size, see about for how this will effect your results. If no bin size is given, the default of 200 is used (ChromHMM's default bin size)
    - Note that, due to the merging process completed in the previous step, `BinarizeBam` will have a very large number of (ignored) errors due to the absence of the header in the `.bam` file being binarized.
        - The number of errors will equate to the number of lines in each of the `.bam` files being binarized, which can be hundreds of Gigabytes large.
        - Do not store the error files for this step.
 ### `5_batch_CreateIncrementalModels.sh` 
 - This uses ChromHMM's `LearnModel` command to learn multiple hidden Markov models for the inputted data with varying number of states. It also obtains the estimated log likelihood value of each model.
    - **IMPORTANT**: The bin size (3rd input) used at this stage **must** be the same as what was given in the previous step.  
    - The number of models and the increment to be used between them is user specified.
    - Note that this is using the 'information' initialisation method for the starting parameter set for the model. As such, the number of states in the model cannot exceed the total number of combinations of marks in your dataset.
        - Though as of version 1.24, the number of states does not always reach the theoretical total combinations of marks. It is likely limited by computational complexity and the observed frequency of each mark in the input data.
            - This only seems to happen when the number of marks is greater than or equal to 5 in testing (at 5, the theoretical maximum should be 2<sup>5</sup> but instead it is 27 in testing).
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
- This uses the binarized data produced by `4_BinarizeBamFiles.sh` to produce a signle very complex model with a high number of user specified states.
    - The model produced will use a random initialisation of emission/transition paramters so that an arbitrary number of states can be used in the model (default initialisation method does not allow for this)
    - These large models were used with blueprint data to obtain sensible threshold values for the similarity metrics used in OptimumNumberOfStates.R
### `CompareModels.sh` 
- This uses ChromHMM's `CompareModels` command to generate comparisons between the models produced in the previous step.
    - Note that `CompareModels` only looks at the emission parameters for each state in the models inputted.
    - This script will compare each model with all of the models that are less complex than it (so a model with 8 states will be compared against only those models that have fewer than 8 states).