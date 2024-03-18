## ============================================================== ##
##                                                                ||
##                            PREAMBLE                            ||
##                                                                ||
## ============================================================== ##
## PURPOSE:                                                       ||
## This annotates the SNPs in a given bim file (for a specific    ||
## chromosome) using an altered UCSC BED format. Specifically to  ||
## be used as an annotation file for LDSC                         ||
## ============================================================== ##
## AUTHOR: Sam Fletcher                                           ||
## CONTACT: s.o.fletcher@exeter.ac.uk                             ||
## CREATED: March 2022                                            ||
## ============================================================== ##
## PREREQUISITES:                                                 ||
## Run 5_batch_CreateIncrementalModels.sh                         ||
## Download/create PLINK files (e.g. 1000 genomes)                ||
## ============================================================== ##
## INPUTS:                                                        ||
## $1 -> Size of the model                                        ||
## $2 -> An altered BED file with columns: start position, end    ||
##       position and annotation for genomic region               ||
## $3 -> The first four columns of the reference bim file         ||
##       (CHR, ID, CM, BP)                                        ||
## $4 -> The full file path for the output .annot file            ||
## ============================================================== ##
## OUTPUTS:                                                       ||
## A .annot file that annotates all SNPS in the given bim file    ||
## ============================================================== ##

## ========== ##
##   SET UP   ##
## ========== ##

rm(list = ls())

if (!requireNamespace("doParallel", quietly = TRUE)) {
  install.packages("doParallel")
}

if (!requireNamespace("foreach", quietly = TRUE)) {
  install.packages("foreach")
}

library(foreach)

if (!requireNamespace("data.table", quietly = TRUE)) {
  install.packages("data.table")
}

arguments <- commandArgs(trailingOnly = TRUE)
model_size <- as.numeric(arguments[1])
bed_file <- arguments[2]
baseline_file <- arguments[3]
output_file_name <- arguments[4]

## ================= ##
##   LOADING FILES   ##
## ================= ##

# data.tables are faster, we are reading from /dev/fd (to avoid temporary files)
# which doesn't work with fread. Hence we read with the base R function
bed_file <- data.table::data.table(read.table(bed_file))
baseline_file <-
  data.table::data.table(read.table(baseline_file, header = TRUE))

## ======================= ##
##   PARALLEL PROCESSING   ##
## ======================= ##

setup_cluster <- function() {
  cores <- parallel::detectCores()
  cluster <- parallel::makeCluster(cores)
  doParallel::registerDoParallel(cluster)
  return(cluster)
}

teardown_cluster <- function(cluster) {
  parallel::stopCluster(cluster)
}

## ============================= ##
##   SNP ANNOTATION ASSIGNMENT   ##
## ============================= ##

# This requires the inputted UCSC BED file to be sorted (which is usually the 
# case) as the binary search relies on this to function.
snp_annotation_binary_search <- function(snp_positions, bed_file) {
  intervals <- unlist(c(bed_file[, 1], tail(bed_file[, 2], 1)))
  
  cluster <- setup_cluster()
  # This is the largest slow down of this script due to the size of the
  # bim and bed files. Parallel processing should speed it up
  indices <-
    foreach::foreach(snp_position = snp_positions, .combine = "c") %dopar% {
      findInterval(snp_position, intervals)
      }
  teardown_cluster(cluster)
  
  state_assignments <- unlist(lapply(indices, function(index) {
    bed_file[index, 3]
  }))
  return(state_assignments)
}


update_baseline_file <- function(row, assignment) {
  state_column <- paste0("state_", assignment)
  data.table::set(row, j = state_column, value = 1)
  return(row)
}

write_snp_annotation <- function(bed_file, baseline_file) {
  assignments <- snp_annotation_binary_search(baseline_file$BP, bed_file)
  
  updated_bim_rows <- lapply(1:nrow(baseline_file), function(row) {
    update_baseline_file(baseline_file[row, ], assignments[[row]])
  })
  annotation_file <- data.table::rbindlist(updated_bim_rows)
  return(annotation_file)
}

## ======== ##
##   MAIN   ##
## ======== ##

# ldsc will break if a state contains SNPs on one chromosome but not on another
# So we need to initialise all of our columns first with zeroes
state_columns <- NULL
for (i in 1:model_size) {
  state_columns <- c(state_columns, paste0("state_", i))
}
baseline_file[, state_columns] <- 0

annotation_file <- write_snp_annotation(bed_file, baseline_file)
data.table::fwrite(annotation_file,
                   file = output_file_name,
                   row.names = FALSE,
                   sep = "\t")

