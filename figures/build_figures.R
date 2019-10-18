library(data.table)
library(tidyverse)
library(ggthemes)

generateFigures <- function(input_name, output_complexity_name, output_time_name, title_description) {
  results <- fread(input_name) %>%
    gather("Method", "Value", -RunTime, -n)
  
  objFunCompPlot <- results %>%
    ggplot(aes(x = n, y = Value, group = Method)) +
    geom_line(aes(color = Method)) +
    scale_color_colorblind() +
    theme_bw() +
    labs(title = paste(title_description, "Analytical and Empirical Query Complexity by Input Size"),
         x = "Input Size (n)",
         y = "Queries") +
    scale_x_continuous(breaks = 1:max(results$n))
  
  runTimePlot <- results %>%
    ggplot(aes(x = n, y = RunTime)) +
    geom_line() +
    theme_bw() +
    labs(title = paste(title_description, "Algorithm Run Time by Input Size"),
         x = "Input Size (n)",
         y = "Run Time (seconds)") +
    scale_x_continuous(breaks = 1:max(results$n))
  
  ggsave(output_complexity_name, plot = objFunCompPlot)
  ggsave(output_time_name, plot = runTimePlot)
}
generateFigures("./output_worst_or.csv", "figure_worst_or_complexity.eps", "figure_worst_or_time.eps", "Worst-case Inputs:")
generateFigures("./output_all_or.csv", "figure_all_or_complexity.eps", "figure_all_or_time.eps", "All Inputs:")