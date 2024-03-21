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
## $1 -> A list of binary files (one for each chromosome)         ||
## $2 -> Bin size used for binarization                           ||
## $3 -> output directory                                         ||
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

arguments <- commandArgs(trailingOnly = TRUE)
binary_file_list <- readLines(arguments[1])
bin_size <- as.numeric(arguments[2])
output_directory <- arguments[3]

## ======================== ##
##   LOADING BINARY FILES   ##
## ======================== ##

binary_files <- lapply(binary_file_list, function(file) {
  read.table(gzfile(file), header = TRUE, skip = 1)
})

names(binary_files) <- unlist(lapply(binary_file_list, function(file) {
  stringr::str_extract(basename(file), "(?<=_).*?(?=_)")
}))

## ===================== ##
##   BED FILE CREATION   ##
## ===================== ##

# This isn't quite correct as it will create BED file that is larger than
# the chromosome.
create_intervals <- function(number_of_intervals, bin_size) {
  start <- seq(
    from = 0,
    to = number_of_intervals * bin_size,
    by = bin_size
  )
  end <- seq(
    from = bin_size,
    to = (number_of_intervals + 1) * bin_size,
    by = bin_size
  )
  return(cbind(start, end))
}

create_bed_file <- function(binary_file, chromosome, bin_size) {
  number_of_rows <- nrow(binary_file)
  chromosome_column <- rep(chromosome, times = number_of_rows)
  interval_columns <- create_intervals(number_of_rows, bin_size)
  return(cbind(chromosome_column, interval_columns, binary_file))
}

## =========== ##
##   OUTPUTS   ##
## =========== ##

write_bed_file <- function(bed_file, chromosome, output_directory) {
  full_file_path <- paste0(output_directory, "/binary_", chromosome, ".bed")
  header_line <- paste("ChromOptimise", chromosome, sep = "\t")
  write(header_line, file = full_file_path)
  write.table(bed_file,
    file = full_file_path,
    append = TRUE,
    sep = "\t",
    row.names = FALSE
  )
}

## ======= ##
##  MAIN   ##
## ======= ##

bed_files <- lapply(names(binary_files), function(file) {
  create_bed_file(binary_files[[file]], file, bin_size)
})

names(bed_files) <- names(binary_files)

lapply(names(bed_files), function(file) {
  write_bed_file(bed_files[[file]], file, output_directory)
})
