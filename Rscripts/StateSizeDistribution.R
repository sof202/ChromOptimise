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
## $2 -> The bin size used with ChromHMM                          ||
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
bin_size <- as.numeric(args[2])
output_directory <- args[3]
plotting_flag <- args[4]

number_of_states <- gsub(".*_([0-9]+)_.*", "\\1", dense_assignment_file)

## ================ ##
##   LOADING FILE   ##
## ================ ##

dense_assignments <- data.table::fread(dense_assignment_file, skip = 1)
dense_assignments <- dplyr::select(dense_assignments, c(V1, V2, V3, V4))
colnames(dense_assignments) <- c("chr", "start", "end", "state")

## ============= ##
##   FUNCTIONS   ##
## ============= ##

create_list_of_sizes <- function(dense_assignments, state_number, bin_size) {
  assignment_sizes <- dense_assignments |>
    dplyr::filter(state == !!state_number) |>
    dplyr::mutate(length = (end - start) / !!bin_size) |>
    dplyr::pull(length)
  return(unlist(assignment_sizes))
}


create_histogram <- function(state_number,
                             dense_assignments,
                             output_directory,
                             bin_size) {
  plot_title <- paste(
    "Size of contiguous regions with state assignment",
    state_number
  )
  x_axis_label <- paste0(
    "length (in bins of length: ",
    bin_size,
    ")"
  )
  sizes <- create_list_of_sizes(dense_assignments, state_number, bin_size)
  plot <-
    ggplot() +
    aes(sizes) +
    geom_histogram(bins = 50) +
    labs(title = plot_title, x = x_axis_label) +
    theme_bw()

  plot_name <- paste0(
    "state_",
    state_number,
    "_contiguous_length_distribution.png"
  )
  suppressMessages(
    ggsave(
      file.path(output_directory, plot_name),
      plot
    )
  )
}


generate_metrics <- function(number_of_states, dense_assignments, bin_size) {
  region_metrics <- data.table::data.table(
    "state" = integer(),
    "mean" = double(),
    "variance" = double()
  )
  for (state in 1:number_of_states) {
    sizes <- create_list_of_sizes(dense_assignments, state, bin_size)
    region_metrics <- rbind(
      region_metrics,
      list(state, median(sizes))
    )
  }
  return(region_metrics)
}

## ======== ##
##   MAIN   ##
## ======== ##

region_metrics <-
  generate_metrics(number_of_states, dense_assignments, bin_size)

data.table::fwrite(
  region_metrics,
  file.path(output_directory, "key_metrics.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

if (!exists("plotting_flag")) {
  plotting_flag <- FALSE
}

if (plotting_flag) {
  options(bitmapType = "cairo")
  invisible(
    lapply(
      1:number_of_states,
      function(state) {
        create_histogram(state, dense_assignments, output_directory, bin_size)
      }
    )
  )
}
