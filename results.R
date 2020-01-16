library(stats)

results <-
  read.csv("csv/2020-01-16_13-39-54_results.csv", stringsAsFactors = FALSE)

results <- results %>%
  group_by(method) %>%
  summarize(
    avg_accuracy = mean(accuracy, na.rm = TRUE),
    standard_deviation = sd(accuracy, na.rm = TRUE)
  )
