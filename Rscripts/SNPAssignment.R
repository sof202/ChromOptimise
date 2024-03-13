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
## $1 -> An altered BED file with columns: start position, end    ||
##       position and annotation for genomic region               ||
## $2 -> The first four columns of the reference bim file         ||
##       (CHR, ID, CM, BP)                                        ||
## $3 -> The full file path for the output .annot file            ||
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
bed_file <- arguments[1]
reference_bim_file <- arguments[2]
output_file_name <- arguments[3]

## ================= ##
##   LOADING FILES   ##
## ================= ##

bed_file <- read.table(bed_file)
reference_bim_file <- read.table(reference_bim_file)

## ============================= ##
##   SNP ANNOTATION ASSIGNMENT   ##
## ============================= ##

# This requires the inputted UCSC BED file to be sorted (which is usually the 
# case) as the binary search relies on this to function.
snp_annotation_binary_search <- function(snp_position, bed_file) {
  intervals <- c(bed_file[, 1], bed_file[nrow(bed_file), 2])
  index <- findInterval(snp_position, intervals)
  if (index > 0) {
    # state assignments should be in third column of bed file
    return(bed_file[index, 3])
  } else {
    # Should never happen
    return(NA)
  }
}

write_snp_annotation <- function(bed_file, bim_file) {
  assignments <- apply(bim_file, 1, function(row) {
    snp_annotation_binary_search(row["BP"], bed_file)
  })
  
  # Optimisation, reduces the number of times this length function is called
  length_assignments <- length(assignments)

  for (row in 1:length_assignments) {
    assignment <- assignments[[row]]
    column_name <- paste0("state_", assignment)
    bim_file[row, column_name] <- 1 
  }
  bim_file[is.na(bim_file)] <- 0
  
  return(bim_file)
}

## ======== ##
##   MAIN   ##
## ======== ##

names(reference_bim_file) <- c("CHR", "SNP", "CM", "BP")

# ldsc expects a particular column order for the annotation file to function
# correctly
column_order <- c("CHR", "BP", "SNP", "CM")
annotation_file <- write_snp_annotation(bed_file,
                                        reference_bim_file[, column_order])

write.table(annotation_file, output_file_name, row.names = FALSE, sep = "\t")

