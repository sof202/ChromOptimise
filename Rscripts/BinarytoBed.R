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

if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
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

# This isn't quite correct as it will create a  BED file that is larger than
# the chromosome.
create_bed_file <- function(binary_file, chromosome, bin_size) {
  number_of_intervals <- dplyr::n(binary_file)
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

## =========== ##
##   OUTPUTS   ##
## =========== ##

write_bed_file <- function(bed_file, chromosome, output_directory) {
  full_file_path <- paste0(output_directory, "/binary-", chromosome, ".bed")
  # Warnings are suppressed here as the 'appending column names to file'
  # message is exactly the behaviour we want
  suppressWarnings(
    write.table(bed_file,
      file = full_file_path,
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

bed_files <- lapply(names(bed_files), function(file) {
  write_bed_file(bed_files[[file]], file, bin_size)
})

bed_files <- lapply(binary_files, function(file) {
  create_bed_file(file, bin_size)
})

names(bed_files) <- names(binary_files)

invisible(lapply(names(bed_files), function(file) {
  write_bed_file(bed_files[[file]], file, output_directory)
}))
