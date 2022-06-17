library(tidyverse)

real_data_res <- read_csv("../output/real_data_analysis_results.csv")

real_data_res %>%
 summarize(tdnn_mse = mean(tdnn), dnn_mse = mean(dnn)) %>%
  print()