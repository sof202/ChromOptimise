## ============================================================== ##
##                                                                ||
##                            PREAMBLE                            ||
##                                                                ||
## ============================================================== ##
## PURPOSE:                                                       ||
## Converts a list of bim files into a bed file                   ||
## ============================================================== ##
## AUTHOR: Sam Fletcher                                           ||
## CONTACT: s.o.fletcher@exeter.ac.uk                             ||
## CREATED: March 2022                                            ||
## ============================================================== ##
## PREREQUISITES:                                                 ||
## Download PLINK files from Zenodo:                              ||
## https://zenodo.org/records/10515792                            ||
## ============================================================== ##
## INPUTS:                                                        ||
## $1 -> bim file location                                        ||
## $2 -> output file location                                     ||
## ============================================================== ##
## OUTPUTS:                                                       ||
## A UCSC bed file with genomic regions corresponding to the SNPS ||
## in the inputted bim file                                       ||
## ============================================================== ##

## ========== ##
##   SET UP   ##
## ========== ##

rm(list = ls())

if (!requireNamespace("data.table", quietly = TRUE)) {
  install.packages("data.table")
}

if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}

arguments <- commandArgs(trailingOnly = TRUE)
bim_file <- arguments[1]
output_file_path <- arguments[2]

## ===================== ##
##   LOADING BIM FILES   ##
## ===================== ##

# The only data that matters in the bim file is the chromosome number and
# the snp's position, the other columns are not used anywhere
bim_file <-
  data.table::data.table(read.table(bim_file)) |>
  dplyr::select(c(1, 4)) |>
  data.table::setnames(c("chromosome", "snp_position"))

## ===================== ##
##   BED FILE CREATION   ##
## ===================== ##

create_bed_file <- function(bim_file) {
  bed_file <- bim_file |>
    dplyr::mutate(
      chromosome = paste0("chr", chromosome),
      start = snp_position - 1,
      end = snp_position
    ) |>
    dplyr::select(c("chromosome", "start", "end"))
  return(bed_file)
}

## =========== ##
##   OUTPUTS   ##
## =========== ##

write_bed_file <- function(bed_file, output_file_path) {
  # Warnings are suppressed here as the 'appending column names to file'
  # message is exactly the behaviour we want
  suppressWarnings(
    write.table(bed_file,
      file = output_file_path,
      sep = "\t",
      row.names = FALSE,
      col.names = FALSE,
      quote = FALSE
    )
  )
  invisible()
}

## ======= ##
##  MAIN   ##
## ======= ##

bed_file <- create_bed_file(bim_file)

invisible(write_bed_file(bed_file, output_file_path))
