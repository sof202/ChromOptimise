## ============================================================== ##
##                                                                ||
##                            PREAMBLE                            ||
##                                                                ||
## ============================================================== ##
## PURPOSE:                                                       ||
## Produces a series of plots that show the distribution of the   ||
## lenths of contiguous regions assigned the same state by        ||
## ChromHMM. It also spits out some key metrics about these       ||
## distributions                                                  ||
## ============================================================== ##
## AUTHOR: Sam Fletcher                                           ||
## CONTACT: s.o.fletcher@exeter.ac.uk                             ||
## CREATED: September 2024                                        ||
## ============================================================== ##
## PREREQUISITES:                                                 ||
## Learn a model using ChromHMM                                   ||
## ============================================================== ##
## INPUTS:                                                        ||
## $1 -> The dense assignments bed file (chromHMM output)         ||
## $2 -> The number of states used by the model in question       ||
## $3 -> The desired output directory for plots and metrics       ||
## ============================================================== ##
## OUTPUTS:                                                       ||
## Histogram plots for the lengths of contiguous regions with     ||
## the same state assignment                                      ||
## tsv detailing key metrics of distributions                     ||
## ============================================================== ##


## ========== ##
##   SET UP   ##
## ========== ##

library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)

dense_assignment_file <- args[1]
number_of_states <- as.numeric(args[2])
output_directory <- args[3]

## ================ ##
##   LOADING FILE   ##
## ================ ##

dense_assignments <- data.table::fread(dense_assignment_file, skip = 1)
dense_assignments <- dplyr::select(dense_assignments, c(V1, V2, V3, V4))
colnames(dense_assignments) <- c("chr", "start", "end", "state")

## ============= ##
##   FUNCTIONS   ##
## ============= ##

create_list_of_sizes <- function(dense_assignments, state_number) {
  assignment_sizes <- dense_assignments |>
    dplyr::filter(state == !!state_number) |>
    dplyr::mutate(length = end - start) |>
    dplyr::pull(length)
  return(unlist(assignment_sizes))
}


create_histogram <- function(state_number,
                             dense_assignments,
                             output_directory) {
  plot_title <- paste(
    "Size of contiguous regions with state assignment",
    state_number
  )
  sizes <- create_list_of_sizes(dense_assignments, state_number)
  plot <-
    ggplot() +
    aes(sizes) +
    geom_histogram(bins = 50) +
    labs(title = plot_title, x = "length") +
    theme_bw()

  plot_name <- paste0(
    "state_",
    state_number,
    "_contiguous_length_distribution.png"
  )
  ggsave(
    file.path(output_directory, plot_name),
    plot
  )
}


generate_metrics <- function(number_of_states, dense_assignments) {
  region_metrics <- data.table::data.table(
    "state" = integer(),
    "mean" = double(),
    "variance" = double()
  )
  for (state in 1:number_of_states) {
    sizes <- create_list_of_sizes(dense_assignments, state)
    region_metrics <- rbind(
      region_metrics,
      c(state, mean(sizes), var(sizes))
    )
  }
  return(region_metrics)
}

## ======== ##
##   MAIN   ##
## ======== ##

options(bitmapType = "cairo")
lapply(
  1:number_of_states,
  create_histogram,
  dense_assignments,
  output_directory
)

region_metrics <- generate_metrics(number_of_states, dense_assignments)
data.table::fwrite(
  file.path(output_directory, "key_metrics.tsv"),
  region_metrics,
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE
)
