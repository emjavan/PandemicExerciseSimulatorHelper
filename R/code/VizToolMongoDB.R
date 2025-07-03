# Look into the data wrtiten by the visualization tool to mongoDB



#install.packages("jsonlite")  # if not already installed
library(jsonlite)
library(tidyverse)

# MongoDB export is not a proper json export as each collection is {...}{...}{...}
# so cannot use fromJSON with these exports
days_data_3 <- stream_in(file("../public_input_data/days_2025-07-02_test3.json"))

