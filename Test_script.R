libs <- c("tidyverse", "leaflet", "lubridate", "readxl", "scales", "ggraph", "networkD3", "igraph", "expm")
installed_libs <- libs %in% rownames(installed.packages())
if (any(installed_libs == F)) {
  install.packages(libs[!installed_libs])
}
invisible(lapply(libs, library, character.only = T))

## finds all .R and .r files within a folder and sources them
sourceFolder <- function(folder, recursive = FALSE, ...)
{
  files <- list.files(folder, pattern = "[.][rR]$",
                      full.names = TRUE, recursive = recursive)
  if (!length(files))
    stop(simpleError(sprintf('No R files in folder "%s"', folder)))
  src <- invisible(lapply(files, source, ...))
  message(sprintf('%s files sourced from folder "%s"', length(src), folder))
}
sourceFolder("./R")





