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
output_directory <- args[2]
model_size <- as.numeric(args[3])
bin_size <- as.numeric(args[4])
plotting_flag <- args[5]

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
    dplyr::mutate(length = (end - start)) |>
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
  # Dividing by bin size to make histograms more interpretable
  sizes <- create_list_of_sizes(dense_assignments, state_number) / bin_size
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


generate_metrics <- function(model_size, dense_assignments) {
  region_metrics <- data.table::data.table(
    "state" = integer(),
    "median" = double()
  )
  for (state in 1:model_size) {
    sizes <- create_list_of_sizes(dense_assignments, state)
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
  generate_metrics(model_size, dense_assignments)

output_file_name <- paste0("Contiguous_state_length_model-", model_size, ".txt")

data.table::fwrite(
  region_metrics,
  file.path(output_directory, output_file_name),
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
      1:model_size,
      function(state) {
        create_histogram(state, dense_assignments, output_directory, bin_size)
      }
    )
  )
}
