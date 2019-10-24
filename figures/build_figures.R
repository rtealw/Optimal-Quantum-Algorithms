library(data.table)
library(tidyverse)
library(ggthemes)

generateFigures <- function(input_name, output_complexity_name, output_time_name, title_description) {
  results <- fread(input_name) %>%
    gather("Method", "Value", -RunTime, -n)
  
  objFunCompPlot <- results %>%
    ggplot(aes(x = n, y = Value, group = Method)) +
    geom_line(aes(color = Method), size = 2)  +
    scale_color_colorblind() +
    theme_bw() +
    labs(title = paste("Numerical Performance for", title_description),
         x = "Input Size (n)",
         y = "Queries") +
    theme(axis.title.y = element_text(angle = 0, hjust = 1, size = 14),
          axis.title.x = element_text(size = 14),
          plot.title = element_text(hjust = 0.5, size = 16),
          axis.text = element_text(size = 12),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 12))
  
  runTimePlot <- results %>%
    ggplot(aes(x = n, y = RunTime)) +
    geom_line(size = 2) +
    theme_bw() +
    labs(title = paste("Run Time for", title_description),
         x = "Input Size (n)",
         y = "Run Time (s)") +
    theme(axis.title.y = element_text(angle = 0, hjust = 1, size = 14),
          axis.title.x = element_text(size = 14),
          plot.title = element_text(hjust = 0.5, size = 16),
          axis.text = element_text(size = 12),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 12))
  
  ggsave(output_complexity_name, plot = objFunCompPlot)
  ggsave(output_time_name, plot = runTimePlot)
}
generateFigures("./output_worst_or.csv", "figure_worst_or_complexity.eps", "figure_worst_or_time.eps", "Worst-case Inputs:")
generateFigures("./output_all_or.csv", "figure_all_or_complexity.eps", "figure_all_or_time.eps", "All Inputs:")
