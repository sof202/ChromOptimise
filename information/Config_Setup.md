# How to set up the config files for pipeline
You will need two config files for this pipeline to work, config.txt and config.r.
\
The paths are used by each of the R and bash scripts to aid in organisation of the scripts as they avoid repetition.
\
Alternatively, you could set each path manually in all scripts, but that would be a lot of effort.
\
You must ensure that the directories actually exist on your system before running any of the scripts in the pipeline.
\
The scripts do not create their own directories (except sub-directories) as this could lead to unwanted folder positions if the user incorrectly assigns the paths in the config files.
## Templates for config files
Note that the pipeline was completed with blueprint data in mind, if your data is already downloaded, processed, binarized etc. then the associated lines in the config files will not be required.
### config.txt
```
## Data directories
export MAIN_DIR="path/to/main/directory"
export DOWNLOAD_DIR="${MAIN_DIR}/path/to/downloads/"
export RAW_DIR="${MAIN_DIR}/path/to/raw/data/"
export PROCESSED_DIR="${MAIN_DIR}/path/to/processed/data/"
export SUBSAMPLED_DIR="${MAIN_DIR}/path/to/subsampled/data/"
export BINARY_DIR="${MAIN_DIR}/path/to/binary/data/"
export MODEL_DIR="${MAIN_DIR}/path/to/chromHMM/models/"
export OPTIMUM_STATES_DIR="${MAIN_DIR}/path/to/optimum/states/output/"
export COMPARE_DIR="${MAIN_DIR}/path/to/comparison/files/"
export BIG_MODELS_DIR="${MAIN_DIR}/path/to/big/models/"

## Script directories
export SCRIPTS_DIR="/path/to/this/repository"
export RSCRIPTS_DIR="${SCRIPTS_DIR}/Rscripts"
export SUPPLEMENTARY_DIR="${SCRIPTS_DIR}/supplementary"
export LOG_DIR="${SCRIPTS_DIR}/LogFiles"

## ChromHMM file locations
export CHROMHMM_MAIN_DIR="/path/to/ChromHMM/main/directory"
export CHROMHMM_CHROM_SIZES="${CHROMHMM_MAIN_DIR}/path/to/chromosome/sizes"
```
### config.r
```
## Data Directories
main_dir="path/to/main/directory"
model_dir=paste0(main_dir, "path/to/model/files/")
optimum_states_dir=paste0(main_dir, "path/to/optimum/states/output/")
likelihood_dir=paste0(optimum_states_dir, "path/to/likelihood/files/")
compare_dir=paste0(main_dir, "path/to/comparison/files/")
big_models_dir=paste0(main_dir,"path/to/big/model/files/")

## Plotting directories
transition_plotting_dir=paste0(main_dir,"/path/to/plots")
emission_plotting_dir=paste0(main_dir,"/path/to/plots")

## Redundancy parameters
emissions_threshold=VALUE
transitions_threshold=VALUE
```
