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

library(dplyr)

arguments <- commandArgs(trailingOnly = TRUE)
model_size <- as.numeric(arguments[1])
bed_file <- arguments[2]
reference_bim_file <- arguments[3]
output_file_name <- arguments[4]

## ================= ##
##   LOADING FILES   ##
## ================= ##

bed_file <- data.table::fread(bed_file)
reference_bim_file <- data.table::fread(reference_bim_file)

## ============================= ##
##   SNP ANNOTATION ASSIGNMENT   ##
## ============================= ##

# This requires the inputted UCSC BED file to be sorted (which is usually the 
# case) as the binary search relies on this to function.
snp_annotation_binary_search <- function(snp_positions, bed_file) {
  intervals <- unlist(c(bed_file[, 1], bed_file[nrow(bed_file), 2]))
  indices <- findInterval(snp_positions, intervals)
  result <- ifelse(indices > 0, bed_file[indices, 3], NA)
  return(result)
}

update_bim_file <- function(row, assignment) {
  state_column <- paste0("state_", assignment)
  data.table::set(row, j = state_column, value = 1)
  return(row)
}

write_snp_annotation <- function(bed_file, bim_file) {
  bim_file_rows <- 1:nrow(bim_file)
  assignments <- unlist(lapply(bim_file_rows, function(row) {
    snp_annotation_binary_search(bim_file$BP[row], bed_file)
  }))
  updated_bim_rows <- lapply(bim_file_rows, function(row) {
    update_bim_file(bim_file[row, ], assignments[[row]])
  })
  annotation_file <- data.table::rbindlist(updated_bim_rows)
  return(annotation_file)
}

## ======== ##
##   MAIN   ##
## ======== ##

names(reference_bim_file) <- c("CHR", "SNP", "CM", "BP")

# ldsc expects a particular column order for the annotation file to function
# correctly
column_order <- c("CHR", "BP", "SNP", "CM")
bim_file <- reference_bim_file[, column_order]

# ldsc will break if a state contains SNPs on one chromosome but not on another
# So we need to initialise all of our columns first with zeroes
state_columns <- NULL
for (i in 1:model_size) {
  state_columns <- c(state_columns, paste0("state_", i))
}
bim_file[state_columns] <- 0

annotation_file <- write_snp_annotation(bed_file, bim_file)

write.table(annotation_file, output_file_name, row.names = FALSE, sep = "\t")

