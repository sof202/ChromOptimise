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
## $2 -> The complete path for the output plots                   ||
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
  gsub("\\.sumstats.*", "", basename(file))
}))

## =============================== ##
##   DATA MANIPULATION FUNCTIONS   ##
## =============================== ##

merge_results_files <- function(results_files, target_column) {
  merged_dataframe <- results_files[[1]]$Category
  for (file in 1:length(results_files)) {
    merged_dataframe <-
      cbind(merged_dataframe, results_files[[file]][target_column])
  }
  colnames(merged_dataframe) <- c("Category", names(results_files))
  return(merged_dataframe)
}

# LDSC adds an L2_0 after each category, which we don't want
remove_l2_suffix <- function(results) {
  results$Category <- sub("L2_0", "", results$Category)
  return(results)
}

## =============================== ##
##   P-VALUE THRESHOLD FUNCTIONS   ##
## =============================== ##

bonferroni_correction <- function(results_files, pvalue_threshold) {
  number_of_traits <- length(results_files)
  number_of_annotations <- nrow(results_files[[1]])
  number_of_hypotheses <- number_of_traits * number_of_annotations
  bonferroni_threshold <- pvalue_threshold / number_of_hypotheses
  return(-log10(bonferroni_threshold))
}

fdr_correction <- function(pvalues, fdr_threshold) {
  pvalues <- pvalues[!is.na(pvalues)]
  adjusted_pvalues <- p.adjust(pvalues, method = "BH")
  significant_pvalues <-
    subset(adjusted_pvalues, adjusted_pvalues < fdr_threshold)
  critical_value <- max(significant_pvalues)
  return(-log10(critical_value))
}

## ====================== ##
##   PLOTTING FUNCTIONS   ##
## ====================== ##

create_heatmap_data <- function(results_files, complete = FALSE) {
  enrichment_data <- merge_results_files(results_files, "Enrichment")
  enrichment_data <- remove_l2_suffix(enrichment_data)
  if (!complete) {
    state_assignment_rows <- grepl("^ChromOptimise.*", data$Category)
    enrichment_data <- enrichment_data[state_assignment_rows, ]
  }

  # pivot_longer is so heatmap can be plotted more easily
  enrichment_data <- tidyr::pivot_longer(enrichment_data,
    cols = -Category,
    names_to = "gwas_trait",
    values_to = "Enrichment"
  )
  return(enrichment_data)
}

create_enrichment_heatmap <- function(results_files, complete = FALSE) {
  enrichment_data <- create_heatmap_data(results_files, complete)
  enrichment_heatmap <-
    ggplot(enrichment_data, aes(gwas_trait,
      Category,
      fill = Enrichment,
      label = round(Enrichment, 2)
    )) +
    geom_tile(color = "black") +
    scale_fill_gradient(low = "light green", high = "dark green") +
    geom_text() +
    theme(panel.background = element_rect(fill = "white")) +
    labs(title = "Enrichment of GWAS traits", x = "GWAS trait", y = "State") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  return(enrichment_heatmap)
}

create_pvalue_barplots <-
  function(results_files, pvalue_threshold, fdr_threshold, complete = FALSE) {
    list_of_pvalue_plots <- list()

    bonferroni_threshold <-
      bonferroni_correction(results_files, pvalue_threshold)

    for (file in 1:length(results_files)) {
      data <- results_files[[file]]
      # We remove the base row as it is guaranteed to have a NaN p-value
      data <- data[-1, ]
      fdr_threshold <- fdr_correction(data$Enrichment_p, pvalue_threshold)

      data <- remove_l2_suffix(data)
      if (!complete) {
        state_assignment_rows <- grepl("^ChromOptimise.*", data$Category)
        data <- data[state_assignment_rows, ]
      }
      plot_title <- names(results_files)[[file]]

      data$Enrichment_p <- -log10(data$Enrichment_p)

      bar_plot <-
        ggplot(
          data,
          aes(x = Category, y = Enrichment_p, fill = Enrichment_p)
        ) +
        geom_bar(stat = "identity", color = "black") +
        scale_fill_gradient(low = "light green", high = "dark green") +
        coord_flip() +
        geom_hline(
          yintercept = bonferroni_threshold,
          linetype = "dashed",
          color = "black"
        ) +
        geom_hline(
          yintercept = fdr_threshold,
          linetype = "dashed",
          color = "gray"
        ) +
        labs(
          title = plot_title,
          x = "State",
          y = "-log_10(Enrichment p-value)"
        ) +
        theme(plot.title = element_text(hjust = 0.5))

      list_of_pvalue_plots[[file]] <- bar_plot
    }
    return(list_of_pvalue_plots)
  }


## =========== ##
##   OUTPUTS   ##
## =========== ##

# We create plots that conatin all of the categories as well as just the state
# assignments. This is because the complete heatmap and bar plots can be
# difficult to read.
complete_enrichment_heatmap <-
  create_enrichment_heatmap(results_files, complete = TRUE)
state_enrichment_heatmap <-
  create_enrichment_heatmap(results_files)

# pvalue threshold is arbitrarily chosen to be 0.05
complete_pvalue_barplots <-
  create_pvalue_barplots(results_files, 0.05, 0.05, complete = TRUE)
state_pvalue_barplots <-
  create_pvalue_barplots(results_files, 0.05, 0.05)

names(complete_pvalue_barplots) <- names(results_files)
names(state_pvalue_barplots) <- names(results_files)

# This usually helps to remove errors around being unable
# to start the PNG device
options(bitmapType = "cairo")

setwd(output_directory)

for (plot in names(state_pvalue_barplots)) {
  plot_name <-
    paste0("ChromOptimise_Categories/Enrichment_pvalues_", plot, ".png")
  ggsave(
    plot_name,
    plot = state_pvalue_barplots[[plot]],
    limitsize = FALSE,
    height = (nrow(results_files[[1]]) - 47) / 5
  )
}

ggsave(
  "ChromOptimise_Categories/Enrichment_heatmap.png",
  state_enrichment_heatmap,
  limitsize = FALSE,
  width = length(results_files),
  height = (nrow(results_files[[1]]) - 47) / 5
)

write.table(
  create_heatmap_data(results_files),
  "ChromOptimise_Categories/Enrichments.txt",
  quote = FALSE,
  sep = ","
)

for (plot in names(complete_pvalue_barplots)) {
  plot_name <-
    paste0("All_Categories/Enrichment_pvalues_", plot, ".png")
  ggsave(
    plot_name,
    plot = complete_pvalue_barplots[[plot]],
    limitsize = FALSE,
    height = nrow(results_files[[1]]) / 5
  )
}

ggsave(
  "All_Categories/Enrichment_heatmap.png",
  complete_enrichment_heatmap,
  limitsize = FALSE,
  width = length(results_files),
  height = nrow(results_files[[1]]) / 5
)

write.table(
  create_heatmap_data(results_files, complete = TRUE),
  "All_Categories/Enrichments.txt",
  quote = FALSE,
  sep = ","
)
