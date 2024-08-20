args <- commandArgs(trailingOnly = TRUE)

transitions_file <- args[1]
output_directory <- args[2]

library(ggplot2)

# Transitions files an annoying header (doesn't line up with
# data). We just remove it as it isn't really required
transitions_data <- data.table::fread(transitions_file, skip = 1)
transitions_data <- dplyr::select(transitions_data, -V1)

number_of_states <- length(transitions_data)

# The main diagonals are pretty much always close to one whereas the other
# entries are closer to 0, resulting in the transitions pngs looking like just
# a diagonal line. We fix this by converting the main diagonal to NAs so we can
# still colour them.
for (i in seq_along(transitions_data)) transitions_data[i, i] <- NA

names(transitions_data) <- paste0("state_", seq_along(transitions_data))
transitions_data[["state_from"]] <- colnames(transitions_data)

heatmap_data <- tidyr::pivot_longer(
  transitions_data,
  cols = -state_from,
  names_to = "state_to"
)

transitions_heatmap <-
  ggplot(heatmap_data, aes(x = state_to, y = state_from, fill = value)) +
  geom_tile() +
  theme_bw() +
  # This is to be consistent with ChromHMM's layout of the transitions heatmap
  scale_y_discrete(limits = rev(levels(factor(heatmap_data[["state_from"]])))) +
  coord_fixed(ratio = 1) +
  scale_fill_gradient(
    low = "white",
    high = "blue",
    na.value = "black"
  ) +
  labs(
    x = "State transitioning towards",
    y = "State transitioning away from",
    title = "Transition Parameters with diagonal removed"
  )

options(bitmapType = "cairo")
file_name <- paste0("transitions_", number_of_states, ".png")
ggsave(
  file.path(output_directory, file_name),
  plot = transitions_heatmap
)
