## ============================================================== ##
##                                                                ||
##                            PREAMBLE                            ||
##                                                                ||
## ============================================================== ##
## PURPOSE:                                                       ||
## This creates heatmaps of enrichment and bar plots for p-values ||
## that are found inside the .results files obtained from ldsc    ||
## when looking at partitioned heritability                       ||
## ============================================================== ##
## AUTHOR: Sam Fletcher                                           ||
## CONTACT: s.o.fletcher@exeter.ac.uk                             ||
## CREATED: March 2022                                            ||
## ============================================================== ##
## PREREQUISITES:                                                 ||
## Run 5_batch_CreateIncrementalModels.sh                         ||
## Run partitioned heritability analysis with ldsc                ||
## ============================================================== ##
## INPUTS:                                                        ||
## $1 -> A list of all .results files                             ||
## $2 -> The full path for the output plots                       ||
## ============================================================== ##
## OUTPUTS:                                                       ||
## A heatmap for the enrichment seen across GWAS traits and       ||
## ChromHMM states                                                ||
## A selection of bar plots for p-values seen across ChomHMM      ||
## states for each GWAS trait                                     || 
## ============================================================== ##

## ========== ##
##   SET UP   ##
## ========== ##

rm(list = ls())

if (!requireNamespace("tidyr", quietly = TRUE)) {
  install.packages("tidyr")
}

if (!requireNamespace("ggtext", quietly = TRUE)) {
  install.packages("ggtext")
}

library(ggplot2)

arguments <- commandArgs(trailingOnly = TRUE)
results_file_list <- readLines(arguments[1])
output_directory <- arguments[2]

## ================= ##
##   LOADING FILES   ##
## ================= ##

results_files <- lapply(results_file_list, function(file) {
  read.table(file, header = TRUE)
})

# This is required as the rest of the script iterates over list items.
# If only one file is read into this script, the logic will break as
# results_files will be a data.frame 
if (is.data.frame(results_files)) {
  results_files <- list(results_files)
}

names(results_files) <- unlist(lapply(results_file_list, function(file) {
  strsplit(basename(file),"\\.results")
}))

## ======== ##
##   MAIN   ##
## ======== ##

merge_results_files <- function(results_files, target_column) {
  merged_dataframe <- results_files[[1]]$Category
  for (file in 1:length(results_files)) {
    merged_dataframe <- 
      cbind(merged_dataframe, results_files[[file]][target_column])
  }
  colnames(merged_dataframe) <- c("Category", names(results_files))
  return(merged_dataframe)
}

create_enrichment_heatmap <- function(results_files) {
  enrichment_data <- merge_results_files(results_files, "Enrichment") 

  # pivot_longer is so heatmap can be plotted more easily
  enrichment_data <- tidyr::pivot_longer(enrichment_data,
                                         cols = -Category,
                                         names_to = "gwas_trait",
                                         values_to = "Enrichment")
  enrichment_heatmap <- 
    ggplot(enrichment_data,
          aes(gwas_trait, Category, fill = Enrichment)) +
    geom_tile(color = "black") +
    coord_fixed(ratio = 1) +
    scale_fill_gradient(low = "light green", high = "dark green") +
    theme(panel.background = element_rect(fill = "white")) +
    labs(title = "Enrichment of GWAS traits", x = "GWAS trait", y = "State") +
    theme(axis.text.x = element_markdown(angle = 90, vjust = 0.5, hjust = 1))
  return(enrichment_heatmap)
}

create_pvalue_barplots <- function(results_files, p_value_threshold) {
  list_of_pvalue_plots <- list()
  for (file in 1:length(results_files)) {
    data <- results_files[[file]]
    plot_title <- names(results_files)[[file]]
    data$Enrichment_p <- -log10(data$Enrichment_p)

    bar_plot <-
      ggplot(data,
             aes(x = Category, y = Enrichment_p, fill = Enrichment_p)) +
      geom_bar(stat = "identity", color = "black") +
      scale_fill_gradient(low = "light green", high = "dark green") +
      coord_flip() +
      geom_hline(yintercept = p_value_threshold,
                 linetype = "dashed",
                 color = "black") +
      labs(title = plot_title, x = "Enrichment p-value", y = "State") +
      theme(plot.title = element_text(hjust = 0.5))

    list_of_pvalue_plots[[file]] <- bar_plot
  }
  return(list_of_pvalue_plots)
}


## =========== ##
##   OUTPUTS   ##
## =========== ##

enrichment_heatmap <- create_enrichment_heatmap(results_files)

# pvalue threshold is arbitrarily chosen to be 0.01 (which will be
# 2 after -log10(), change this if you wish)
pvalue_barplots <- create_pvalue_barplots(results_files, 2)
names(pvalue_barplots) <- names(results_files)

setwd(output_directory)
# This usually helps to remove errors around being unable
# to start the PNG device
options(bitmapType='cairo')

for (plot in names(pvalue_barplots)) {
  plot_name <- paste0("Enrichment_pvalues_barplot_", plot,".png")
  ggsave(
    plot_name,
    plot = pvalue_barplots[[plot]]
  )
}

ggsave(
  "Enrichment_heatmap.png",
  plot = enrichment_heatmap
)
