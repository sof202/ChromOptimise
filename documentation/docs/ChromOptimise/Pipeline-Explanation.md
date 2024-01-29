---
sidebar_position: 3
---

# Pipeline explanation

The main pipeline is split into 3 main steps:

- Manipulate the data set using user defined parameters
- Use ChromHMM's pipeline to obtain multiple hidden Markov models on the processed data
- Use statistical analysis to evaluate the model with the optimum number of states

\
The pipeline expects aligned `.bam` files as its input. One may want to download these files from EGA, this can be done using [0_EGADownloading.sh](#0_egadownloadingsh).
\
The pipeline will then output:

- A scatter plot that showcases how the number of states in the hidden Markov model relates to the estimated log likelihood of the model
- A file containing the optimum number of states to use with the binarized data files (when a specific bin size has been used)

For a visual representation of the pipeline, please consult the file: [`Optimal_States_Analysis_Pipeline.pdf`](https://github.com/sof202/ChromOptimise/blob/main/information/Optimal_States_Analysis_Pipeline.pdf).

## 0_EGADownloading.sh

This is an artefact of the pipeline originally being built for the blueprint data from EGA.
The script will download files from a file of file names using the pyega3 python package.
There are a few prerequisites to use this properly:

Step one: Login credentials need to be saved in `.json` format:

```json
{
    "username":"username@domain.com",
    "password":"password"
}
```

Step Two: Create a conda environment with pyega3 installed:

```shell
conda create -n myenv pyega3
```

Step Three: Change the following paths:

- Path to file of file names
- Path to conda config file (`.../etc/profile.d/conda.sh`)
- Path to login credentials created in step 1
- Path to conda environment created in step 2

Example:

```shell
# Downloads all files listed in FileOfFileNames.txt from EGA
sbatch 0_EGADownloading.sh FileOfFileNames.txt
```

## 1_MoveFilesToSingleDirectory.sh

Again, this is an artefact of the pipeline originally being built for the blueprint data from EGA.

The downloaded files from EGA will be in separate directories (as each directory has its own checksum file to validate the downloads). For easier file manipulation later in the pipeline, these files are collated into folders denoted by the epigenetic mark they contain information on.

:::warning[Important]
This script relies on the `.bam` files containing the epigenetic mark in the file names. Please ensure this is the case for your raw files.
:::

Example:

```shell
# Moves all files containing 'H3K27me3' to: ${RAW_DIR}/H3k27me3
sbatch 1_MoveFilesToSingleDirectory.sh H3K27me3
```

**Note**: If this step is skipped, please ensure that the raw `.bam` files are organised into folders with memorable names (epigenetic mark name is reccomended).

## 2_batch_ProcessBamFiles.sh

If the aligned `.bam` files you are working with have already been processed, skip this step.

The purpose of this step is to sort, index, filter and remove the duplicates of the aligned `.bam` files.
\
Note that, in the blueprint data, the aligned read files are single end and duplicates have already been marked for duplicate reads using picard.
\
Therefore, the processing stage will need to be adapted if your `.bam` files don't have these properties.

This step also outputs some statistics on the aligned read files using `samtools stats` and `samtools idxstats` before and after processing.

Example:

```shell
# Processes all files in ${RAW_DIR}/H3k27me3 using a Phred score threshold
# of 30
sbatch --array=1-4 2_batch_ProcessBamFiles.sh H3K27me3 30
```

## 3_SubsampleBamFiles.sh

This merges all of the processed `.bam` files from the previous step using `samtools merge` and subsequently subsamples this file using `samtools view -s`.
\
**Note**: samtools had a known bug in version 0.1.18 where subsampling would fail if the sample size was above 50%, ensure that your version of samtools satisfies the [software requirements](../intro.md#software-requirements).

The reason for merging files and subsampling instead of simply subsampling the original files is that it is unlikely that every processed `.bam` file is of the same size. If the `.bam` files are of different sizes, then different subsamples with the same number of files can have different sizes.
\
This results in the different samples no longer being directly comparable.

**Note**: If you are working with a lot of files, this step can use a large amount of memory due to the way samtools merges files, see [memory profiling](./Memory-Profiling.md) for more information.

Example:

```shell
# Merges all files in ${PROCESSED_DIR}/H3k27me3 and subsamples this file
# with sample size 75%
sbatch 3_SubsampleBamFiles.sh H3K27me3 75
```

## 4_BinarizeFiles.sh

This uses ChromHMM's `BinarizeBam` command to binarize the subsampled `.bam` file produced in the previous step.
\
It is assumed at this point that the above steps have been completed for all of the epigenetic marks that you aim to use in the hidden Markov models, and that the sample size used for each mark is the same.
\
This step takes a user input for the bin size, this will have a noticable effect on the quality of your models, see [this wiki page](./Factors-that-affect-the-output.md#bin-size) for further information.

**Note**: This file has commented lines at the end. Uncomment these lines to make the script remove binary files associated with mitochondrial DNA.

Example:

```shell
# Takes all files in ${SUBSAMPLE_DIR} that use a sample size of 75% and
# binarizes them using a bin size of 200 and the hg19 assembly
sbatch 4_BinarizeFiles.sh 200 75 hg19
```

## 5_batch_CreateIncrementalModels.sh

This uses ChromHMM's `LearnModel` command to learn multiple hidden Markov models. It also obtains the [estimated log likelihood](./ChromHMM-overview.md#estimated-log-likelihood) value of each model.

**Note**: This script uses the 'information' initialisation method for the starting parameter set for the model. Due to this, the number of states in the model cannot exceed the total number of combinations of marks in your dataset.
\
One may not know this number before running the script. However, the maximum number of mark combinations one can have for a binary file will always be less than or equal to $2^k$, where $k$ is the total number of marks. If one exceeds the maximum number of states permitted for this initialisation technique, ChromHMM will tell the user the maximum number of states allowed in the associated error message.

This script will always generate a model with 2 states. Inspecting the emission parameters for this simple model is a good way of validating your data.
\
Generally, most of your genomic data will be non-coding regions. The 2 state model should reflect this property in the overlap files produced.

Example:

```shell
# Generates 6 models using the hg19 assembly (sizes 2-7) using 
# the binarized files with bin size 200 and sample size 75%
# --array=1-3 is used so that each array element learns 2 models (optimal).
sbatch --array=1-3 5_batch_CreateIncementalModels.sh 6 1 200 75 hg19
```

## 6_OptimalNumberOfStates.sh

This is the statistical analysis step where the optimal number of states to use with the model is determined (subject to many [factors](./Factors-that-affect-the-output.md)).

This script uses four R scripts:

- PlotLikelihoods.R
  - This creates a plot of estimated log likelihood against number of model states. This gives the user feedback on how the accuracy of the model increases as one increases the number of states used in the model.
- OptimumNumberOfStates.R
  - This uses vector similarity (Euclidean distance) to identify how similar each pair states are in terms of their emission parameters in each model.
  - The script also looks at the maximum transition probability towards each state to determine the relative assignment of each state (a low maximum transition probability implies that the [forwards-backwards algorithm](https://github.com/sof202/ChromOptimise/wiki/ChromHMM-overview#estimated-log-likelihood) assigns the state to relatively few genomic bins).
    - If a state scores low enough in both of these metrics (threshold is defined in [config.R](./Configuration-Files-Setup.md#configr)), then the state is classed as redundant.
  - Models with redundant states are rejected in favour of less complex models.
- CalculateAIC.R
  - This script will calculate the [Akaike information critereon](https://en.wikipedia.org/wiki/Akaike_information_criterion) for each of the models created in [5_batch_CreateIncremenatlModels.sh](#5_batch_createincrementalmodelssh).
  - AIC is a heuristic for probabilistic models that attempts to explain how accurate a model is in spite of its complexity.
    - It is given by: $AIC=2k-2ln(L)$.
      - Where $k$ is the number of parameters in the model and $L$ is the estimated likelihood.
    - If two models have the same log likelihood, but one has more states, the model with fewer states will have a lower (and therefore better) AIC.
  - Considering the estimated log likelihood is generally very low for ChromHMM, these values are very large. To combat this, we employ the relative AIC measure: $AIC_{relative} = exp((AIC_{min}-AIC_i)/2)$. This measure is proportional to the probability that the $i^{th}$ model minimizes the (estimated) information loss ([source](https://link.springer.com/book/10.1007/b97636)).
- CalculateBIC.R
  - The Bayesian information critereon is very similar to AIC. The difference here is that this is much more punishing on the number of parameters in the model when the dataset is very large (which it will be in the case of genomic datasets).
  - The Bayesian information critereon is given by: $BIC = ln(n)k - 2ln(L)$
    - Where $n$ is the total number of observations, $k$ is the number of parameters in the model and $L$ is the estimated likelihood.
  - Again, these values are very large and hard to compare. As a result we find the BIC of each model relative to the minimum BIC of all models (this time using a simple ratio as there is no merit to using $exp$ with this heuristic).
    - As with AIC, a smaller value for BIC generally indicates a better model.
  - It is important to note that AIC and BIC are heuristics, not metrics. They only approximate how good a model is relative to its complexity. As such the graphics and .csv files that come from these scripts are in place to provide further information to the user (Perhaps a model with 6 states is better under this metric than a model with 8 states, despite the more complex model having no 'redundant states').

Example:

```shell
# Generates metrics and heuristics for the models currently in ${MODEL_DIR}
# and uses these to determine  the optimal number of states for the binarized
# data set.
sbatch 6_OptimalNumberOfStates.sh
```

## Note

This pipeline was built using the SLURM workload manager. Scripts will result in errors or even fail entirely if not ran using the `sbatch` SLURM command. Please adapt the scripts accordingly if SLURM is not being used. For information on how to do this, please consult [the page on SLURM](./SLURM-Workload-Manager-Information.md).
