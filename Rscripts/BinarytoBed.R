## ============================================================== ##
##                                                                ||
##                            PREAMBLE                            ||
##                                                                ||
## ============================================================== ##
## PURPOSE:                                                       ||
## Converts binary files produced by ChromHMM into BED format     ||
## Warning: This script assumes that missing data is the same as  ||
## the absence of the mark. It also doesn't take into account the ||
## total chromosome size (which shouldn't end up mattering)       ||
## ============================================================== ##
## AUTHOR: Sam Fletcher                                           ||
## CONTACT: s.o.fletcher@exeter.ac.uk                             ||
## CREATED: March 2022                                            ||
## ============================================================== ##
## PREREQUISITES:                                                 ||
## Run 4_BinarizeBamFiles.sh                                      ||
## ============================================================== ##
## INPUTS:                                                        ||
## $1 -> Binary file location (from ChromHMM)                     ||
## $2 -> Bin size used for binarization                           ||
## $3 -> Chromosome number                                        ||
## $4 -> output file location                                     ||
## ============================================================== ##
## OUTPUTS:                                                       ||
## A UCSC bed file with columns: chromosome, start, end, then     ||
## Boolean columns signifying the presence/absence of each mark   ||
## in the binary file                                             ||
## ============================================================== ##

## ========== ##
##   SET UP   ##
## ========== ##

rm(list = ls())

if (!requireNamespace("stringr", quietly = TRUE)) {
  install.packages("stringr")
}

if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}

arguments <- commandArgs(trailingOnly = TRUE)
binary_file <- arguments[1]
bin_size <- as.numeric(arguments[2])
chromosome <- arguments[3]
output_file_path <- arguments[4]

## ======================== ##
##   LOADING BINARY FILES   ##
## ======================== ##

# We skip the first line as the meta data has 2 entries which can result in
# reading errors (if number of marks is not equal to 2).
binary_file <-
  data.table::data.table(read.table(binary_file,
    header = TRUE,
    skip = 1
  ))


## ===================== ##
##   BED FILE CREATION   ##
## ===================== ##

# This isn't quite correct as it will create a  BED file that is larger than
# the chromosome.
create_bed_file <- function(binary_file, chromosome, bin_size) {
  number_of_intervals <- nrow(binary_file)
  bed_file <- binary_file |>
    dplyr::mutate(
      chromosome = chromosome,
      start = seq(
        from = 0,
        to = (number_of_intervals - 1) * bin_size,
        by = bin_size
      ),
      end = seq(
        from = bin_size,
        to = number_of_intervals * bin_size,
        by = bin_size
      ),
      .before = 1
    )
  return(bed_file)
}

## ======= ##
##  MAIN   ##
## ======= ##

bed_file <- create_bed_file(binary_file, chromosome, bin_size)

data.table::fwrite(bed_file,
  file = output_file_path,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)
