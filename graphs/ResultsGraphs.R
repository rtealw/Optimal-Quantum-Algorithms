library(data.table)
library(ggplot2)
library(ggthemes)

results <- fread("ORWorstCaseOutput.csv") %>%
  gather("Method", "Value", -RunTime, -n)

objFunCompPlot <- results %>%
  ggplot(aes(x = n, y = Value, group = Method)) +
  geom_line(aes(color = Method)) +
  scale_color_colorblind() +
  theme_bw() +
  labs(title = "Comparison of Calculated Complexity vs. True Values",
       x = "Input Size n",
       y = "Optimal Quantum Query Complexity") +
  scale_x_continuous(breaks = 1:max(results$n))

runTimePlot <- results %>%
  ggplot(aes(x = n, y = RunTime)) +
  geom_line() +
  theme_bw() +
  labs(title = "Algorithm Runtime in seconds",
       x = "Input Size n",
       y = "Run Time (s)") +
  scale_x_continuous(breaks = 1:max(results$n))

ggsave("ObjFunCompToTrue_OR.eps", plot = objFunCompPlot)
ggsave("WorstCaseRunTimePlot_OR.eps", plot = runTimePlot)

