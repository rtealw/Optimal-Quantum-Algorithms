library(data.table)
library(tidyverse)
library(ggthemes)

library(devtools)
# devtools::install_github("backlin/treesmapstheorems")
library(treesmapstheorems)

generateFigures <- function(input_name, output_complexity_name, output_time_name, title_description) {
  results <- fread(input_name) %>%
    gather("Method", "Value", -RunTime, -n)
  
  custom_theme <- theme(axis.title.y = element_text(color = "black"),
                        axis.title.x = element_text(color = "black"),
                        plot.title = element_text(hjust = 0.5, color = "black"),
                        axis.text.x = element_text(color = "black"),
                        axis.text.y = element_text(color = "black"),
                        legend.title = element_text(color = "black"),
                        legend.text = element_text(color = "black")) +
                  treesmapstheorems::tmt_theme()

  objFunCompPlot <- results %>%
    ggplot(aes(x = n, y = Value, group = Method)) +
    geom_line(aes(color = Method), size = 2)  +
    scale_color_colorblind() +

    labs(title = paste("Numerical Performance for", title_description),
         x = "Input Size (n)",
         y = "Queries") +
    scale_x_continuous(breaks = 1:max(results$n)) +
    custom_theme


  runTimePlot <- results %>%
    ggplot(aes(x = n, y = RunTime)) +
    geom_line(size = 2) +
    labs(title = paste("Run Time for", title_description),
         x = "Input Size (n)",
         y = "Run Time (s)") +
    scale_x_continuous(breaks = 1:max(results$n)) +
    custom_theme

  ggsave(output_complexity_name, plot = objFunCompPlot)
  ggsave(output_time_name, plot = runTimePlot)
}
generateFigures("./output_worst_or.csv", "figure_worst_or_complexity.eps", "figure_worst_or_time.eps", "Worst-case Inputs")
generateFigures("./output_all_or.csv", "figure_all_or_complexity.eps", "figure_all_or_time.eps", "All Inputs")
