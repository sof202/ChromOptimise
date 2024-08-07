required_packages <- c("tidyr", "dplyr", "data.table", "ggplot2")

installed_packages <- unlist(lapply(.libPaths(), list.files))

install.packages(setdiff(required_packages, installed_packages))
