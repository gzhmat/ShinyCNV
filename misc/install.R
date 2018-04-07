list.of.packages = c("rstudioapi", "shiny", "magrittr", "scales", "dplyr", "data.table", "DT", "readr")
new.packages = list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages))
  install.packages(new.packages)
