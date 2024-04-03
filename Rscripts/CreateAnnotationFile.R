## ============================================================== ##
##                                                                ||
##                            PREAMBLE                            ||
##                                                                ||
## ============================================================== ##
## PURPOSE:                                                       ||
## Appends columns onto the end of an annotation file             ||
## ============================================================== ##
## AUTHOR: Sam Fletcher                                           ||
## CONTACT: s.o.fletcher@exeter.ac.uk                             ||
## CREATED: March 2022                                            ||
## ============================================================== ##
## PREREQUISITES:                                                 ||
## Download Baseline files from Zenodo:                           ||
## https://zenodo.org/records/10515792                            ||
## ============================================================== ##
## INPUTS:                                                        ||
## $1 -> baseline annotation file                                 ||
## $2 -> file of state assignments                                ||
## $3 -> file of mark assignments                                 ||
## $4 -> size of the model                                        ||
## $5 -> output file location                                     ||
## ============================================================== ##
## OUTPUTS:                                                       ||
## An annotation file to use with ldsc                            ||
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
baseline_annotation <- arguments[1]
state_assignments <- arguments[2]
mark_assignment <- arguments[3]
model_size <- as.numeric(arguments[4])
output_file_path <- arguments[5]

## ================ ##
##   FILE LOADING   ##
## ================ ##

baseline_annotation <- data.table::data.table(
  read.table(baseline_annotation, header = TRUE)
)
state_assignments <- as.numeric(readLines(state_assignments))
mark_assignment <- data.table::data.table(
  read.table(mark_assignment, header = TRUE)
)

## ======== ##
##   MAIN   ##
## ======== ##

generate_state_columns <- function(model_size, state_assignments) {
  binarized_state_columns <- matrix(0,
    nrow = length(state_assignments),
    ncol = model_size
  )
  binarized_state_columns[cbind(
    seq_len(length(state_assignments)),
    state_assignments
  )] <- 1
  binarized_state_columns <- data.table::as.data.table(binarized_state_columns)
  names(binarized_state_columns) <- paste0("Chromoptimise_state_", 1:model_size)

  return(binarized_state_columns)
}

binarized_state_columns <- generate_state_columns(model_size, state_assignments)
output_annotation <-
  cbind(baseline_annotation, binarized_state_columns, mark_assignment)

# These regions end up causing a lot of trouble as for some reason their
# enrichment is massively negative for all traits considered thus far (>80).
blacklisted_annotations <-
  c("MAF_Adj_LLD_AFR", "MAF_Adj_ASMC", "MAF_Adj_Predicted_Allele_Age")
output_annotation <- output_annotation |>
  dplyr::select(-all_of(blacklisted_annotations))

write.table(
  output_annotation,
  file = output_file_path,
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE,
  sep = "\t"
)
