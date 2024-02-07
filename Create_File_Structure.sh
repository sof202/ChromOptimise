#!/bin/bash

## =================================================================================##
##                                                                                  ||
##                                     PREAMBLE                                     ||
##                                                                                  ||
## =================================================================================##
## PURPOSE:                                                                         ||
## Creates the file structure based off of the file structure that is given in      ||
## the FilePaths.txt file in /configuration                                         ||
## =================================================================================##
## AUTHOR: Sam Fletcher                                                             ||
## CONTACT: s.o.fletcher@exeter.ac.uk                                               ||
## CREATED: January 2024                                                            ||
## =================================================================================##

configuration_directory=$1

source "${configuration_directory}/FilePaths.txt"


# Create each data directory sequentially
mkdir -p "${DOWNLOAD_DIR}"
mkdir -p "${RAW_DIR}"
mkdir -p "${PROCESSED_DIR}"
mkdir -p "${SUBSAMPLED_DIR}"
mkdir -p "${BINARY_DIR}"
mkdir -p "${MODEL_DIR}"
mkdir -p "${OPTIMUM_STATES_DIR}"
mkdir -p "${COMPARE_DIR}"
mkdir -p "${BIG_MODELS_DIR}"