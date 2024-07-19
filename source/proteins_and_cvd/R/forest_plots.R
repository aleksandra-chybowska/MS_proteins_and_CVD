library(tidyverse)
library(tidytext)

# setwd(
#   paste0("C:/Users/s1654019/Projects/python/proteins/results/_old/",
#       "incremental_parallel_correct/hosp/survival_analysis_manual/")
# )
# 
# cvd_and_deaths = read_csv("Merged_results_deaths_and_CVD.csv")

setwd("~/Projects/python/proteins/results/incremental_parallel/hosp/agesex/")
cvd_and_deaths = read_csv("merged_results_bonf_significant_full.csv")
x = cvd_and_deaths
dict = c("composite_CVD" = "Composite CVD",
         "CVD_death" = "CVD death",
         "death" = "Death",
         "hf" = "Heart Failure")

x$Outcome <- dict[x$event]
x <- x %>%
  group_by(Outcome) %>%
  mutate(Name = reorder_within(name, hr, Outcome))

stacked <- ggplot(x, aes(y = hr, x = Name)) + 
  geom_point(size = 2, position = position_dodge(0.5), color = "#4297d8") +
  geom_errorbar(aes(ymin = lci, ymax = uci), color = "#4297d8",
                position = position_dodge(0.5), width = 0.1) +
  labs(x = "", y = "Hazard Ratio per SD [95% Confidence Interval]", 
       col = "Model covariates") +
  theme_bw() +
  geom_hline(yintercept = 1, linetype = "dotted") +
  theme(axis.text.x = element_text(size = 8, vjust = 0.5), 
        axis.text.y = element_text(size = 6), 
        legend.position = "bottom",
        plot.title = element_text(size = 8), 
        legend.title = element_text(hjust = 0.5)) +
  scale_x_reordered() +
  coord_flip() +
  facet_wrap(~Outcome, scales = "free_y")

print(stacked)
