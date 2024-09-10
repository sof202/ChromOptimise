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

library(ggplot2)

arguments <- commandArgs(trailingOnly = TRUE)
results_file_list <- readLines(arguments[1])
cell_type <- arguments[2]
pvalue_threshold <- as.numeric(arguments[3])
output_directory <- arguments[4]

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
  gsub("\\.results", "", basename(file))
}))

## =============================== ##
##   DATA MANIPULATION FUNCTIONS   ##
## =============================== ##

merge_results_files <- function(results_files, target_column) {
  merged_dataframe <- results_files[[1]]$Category
  for (file in seq_along(results_files)) {
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

pivot_enrichment_data <- function(enrichment_data, column_name) {
  # pivot_longer is so heatmap can be plotted more easily
  enrichment_data <- tidyr::pivot_longer(enrichment_data,
    cols = -Category,
    names_to = "gwas_trait",
    values_to = "value"
  )
  colnames(enrichment_data)[[ncol(enrichment_data)]] <- column_name
  return(enrichment_data)
}

## =============================== ##
##   QUALITY ASSURANCE FUNCTIONS   ##
## =============================== ##

# The number of hypotheses is technically the number of traits multiplied by
# the number of annotations. However, the user is not able to control the
# number of annotations that come from the baseline annotation files. Thus
# we discount these here when calculating bonferroni correction.
get_bonferroni_threshold <- function(results_files, pvalue_threshold) {
  number_of_traits <- length(results_files)
  relevant_categories <-
    sum(grepl(paste0(cell_type, "_*"), results_files[[1]]$Category))
  number_of_hypotheses <- number_of_traits * relevant_categories
  bonferroni_threshold <- pvalue_threshold / number_of_hypotheses
  return(-log10(bonferroni_threshold))
}

get_fdr_threshold <- function(pvalues, pvalue_threshold) {
  pvalues <- pvalues[!is.na(pvalues)]
  adjusted_pvalues <- p.adjust(pvalues, method = "BH")
  significant_pvalues <-
    subset(adjusted_pvalues, adjusted_pvalues < pvalue_threshold)

  # It is possible for significant_pvalues to be empty if the adjusted_pvalues
  # are too great
  if (length(significant_pvalues) == 0) {
    critical_value <- pvalue_threshold
  } else {
    critical_value <- max(significant_pvalues)
  }
  return(-log10(critical_value))
}

negative_enrichment_proportion <- function(results_files) {
  enrichment_data <- merge_results_files(results_files, "Enrichment")
  enrichment_data <- pivot_enrichment_data(enrichment_data, "Enrichment")
  number_of_negative_values <- sum(enrichment_data$Enrichment < 0, na.rm = TRUE)
  return(number_of_negative_values / nrow(enrichment_data))
}

write_poor_enrichment_warning <- function(results_files) {
  negative_enrichment <- negative_enrichment_proportion(results_files)
  if (negative_enrichment > 0.05) {
    lines <- c(
      "WARNING: A significant proportion of enrichments are negative.",
      "The proportion of negative enrichements for this run was:",
      negative_enrichment,
      "Please check the enrichment heatmap to gain more information.",
      "The wiki gives information on how to avoid this."
    )
    writeLines(
      lines,
      file.path(output_directory, "WARNING.txt")
    )
  }
  invisible()
}

## ====================== ##
##   PLOTTING FUNCTIONS   ##
## ====================== ##

create_heatmap_data <- function(results_files, complete = FALSE) {
  enrichment_data <- merge_results_files(results_files, "Enrichment")
  enrichment_data <- remove_l2_suffix(enrichment_data)

  enrichment_p_data <- merge_results_files(results_files, "Enrichment_p")
  enrichment_p_data <- remove_l2_suffix(enrichment_p_data)

  if (!complete) {
    chromoptimise_rows <-
      grepl(paste0(cell_type, "_*"), enrichment_data[["Category"]])
    enrichment_data <- enrichment_data[chromoptimise_rows, ]
    enrichment_p_data <- enrichment_p_data[chromoptimise_rows, ]
  }

  enrichment_data <- pivot_enrichment_data(enrichment_data, "Enrichment")
  enrichment_p_data <- pivot_enrichment_data(enrichment_p_data, "Enrichment_p")

  heatmap_data <-
    cbind(enrichment_data, enrichment_p_data[, ncol(enrichment_p_data)])
  return(heatmap_data)
}

create_enrichment_heatmap <- function(results_files,
                                      pvalue_threshold,
                                      complete = FALSE) {
  heatmap_data <- create_heatmap_data(results_files, complete)

  bonferroni_threshold <-
    get_bonferroni_threshold(results_files, pvalue_threshold)
  fdr_threshold <-
    get_fdr_threshold(heatmap_data$Enrichment_p, pvalue_threshold)

  heatmap_data[is.na(heatmap_data)] <- 1
  heatmap_data <- dplyr::mutate(
    heatmap_data,
    Enrichment_p = -log10(Enrichment_p)
  )

  # Enrichment values greater than 100 have always been within annotations
  # that have lots of negative enrichments in my experience.
  heatmap_data <- dplyr::mutate(
    heatmap_data,
    Enrichment = dplyr::if_else(Enrichment < 0, NA_real_, Enrichment),
    Enrichment = dplyr::if_else(Enrichment > 100, NA_real_, Enrichment)
  )

  # we mainly care about how the enrichment compares to 1. So we scale the
  # colours of our heatmap around this value
  max_enrichment <- max(heatmap_data[["Enrichment"]], na.rm = TRUE)
  colour_palette <- c("yellow", "pink", "red", "grey")
  colour_scale <- c(0, 1, max_enrichment - 0.5, max_enrichment) / max_enrichment

  enrichment_heatmap <-
    ggplot(heatmap_data, aes(
      x = gwas_trait,
      y = Category,
      fill = Enrichment,
    )) +
    geom_tile(color = "black") +
    scale_fill_gradientn(
      colors = colour_palette,
      values = colour_scale,
      na.value = "grey"
    ) +
    geom_text(
      aes(
        label =
          ifelse(Enrichment_p > !!bonferroni_threshold, "**",
            ifelse(Enrichment_p > !!fdr_threshold, "*", "")
          )
      ),
      color = "black"
    ) +
    theme(panel.background = element_rect(fill = "white")) +
    labs(
      title = "Enrichment of GWAS traits",
      x = "GWAS trait",
      y = "Category"
    ) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  return(enrichment_heatmap)
}

create_pvalue_barplots <- function(results_files,
                                   pvalue_threshold,
                                   complete = FALSE) {
  list_of_pvalue_plots <- list()

  bonferroni_threshold <-
    get_bonferroni_threshold(results_files, pvalue_threshold)

  for (file in seq_along(results_files)) {
    barplot_data <- results_files[[file]]
    barplot_data <- remove_l2_suffix(barplot_data)

    fdr_threshold <-
      get_fdr_threshold(barplot_data[["Enrichment_p"]], pvalue_threshold)

    if (!complete) {
      state_assignment_rows <-
        grepl(paste0(cell_type, "_*"), barplot_data[["Category"]])
      barplot_data <- barplot_data[state_assignment_rows, ]
    }

    barplot_data <- dplyr::mutate(
      barplot_data,
      Enrichment_p = -log10(Enrichment_p)
    )

    plot_title <- names(results_files)[[file]]
    bar_plot <-
      ggplot(
        barplot_data,
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
        x = "Category",
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
complete_heatmap <-
  create_enrichment_heatmap(results_files, pvalue_threshold, complete = TRUE)
chromoptimise_only_heatmap <-
  create_enrichment_heatmap(results_files, pvalue_threshold, complete = FALSE)

complete_barplots <-
  create_pvalue_barplots(results_files, pvalue_threshold, complete = TRUE)
chromoptimise_only_barplots <-
  create_pvalue_barplots(results_files, pvalue_threshold, complete = FALSE)

names(complete_barplots) <- names(results_files)
names(chromoptimise_only_barplots) <- names(results_files)

# This usually helps to remove errors around being unable
# to start the PNG device
options(bitmapType = "cairo")

for (plot in names(chromoptimise_only_barplots)) {
  plot_name <-
    paste0("Enrichment_pvalues_", plot, ".png")
  plot_file_path <-
    file.path(output_directory, "ChromOptimise_Categories", plot_name)
  ggsave(
    plot_file_path,
    plot = chromoptimise_only_barplots[[plot]],
    limitsize = FALSE,
    height = (nrow(results_files[[1]]) - 47) / 5
  )
}

ggsave(
  file.path(
    output_directory,
    "ChromOptimise_Categories",
    "Enrichment_heatmap.png"
  ),
  chromoptimise_only_heatmap,
  limitsize = FALSE,
  width = max(length(results_files), 10),
  height = (nrow(results_files[[1]]) - 47) / 5
)

write.table(
  create_heatmap_data(results_files),
  file.path(output_directory, "ChromOptimise_Categories", "Enrichments.csv"),
  quote = FALSE,
  row.names = FALSE,
  sep = ","
)

for (plot in names(complete_barplots)) {
  plot_name <-
    paste0("Enrichment_pvalues_", plot, ".png")
  plot_file_path <-
    file.path(output_directory, "All_Categories", plot_name)
  ggsave(
    plot_file_path,
    plot = complete_barplots[[plot]],
    limitsize = FALSE,
    height = nrow(results_files[[1]]) / 5
  )
}

ggsave(
  file.path(output_directory, "All_Categories", "Enrichment_heatmap.png"),
  complete_heatmap,
  limitsize = FALSE,
  width = length(results_files),
  height = nrow(results_files[[1]]) / 5
)

write.table(
  create_heatmap_data(results_files, complete = TRUE),
  file.path(output_directory, "All_Categories", "Enrichments.csv"),
  quote = FALSE,
  sep = ","
)

write_poor_enrichment_warning(results_files)
