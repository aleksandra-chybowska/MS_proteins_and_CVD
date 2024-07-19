library(tidyverse)

setwd("~/Projects/python/proteins/results/attenuation_factor/agesex/hosp/")

cvd = read_csv("attenuation_factor_composite_CVD_agesex.csv")
death = read_csv("attenuation_factor_death_agesex.csv")
bonf = 0.05/439


cvd = subset(cvd, p>bonf)
summary = cvd %>% group_by(formula) %>%
  summarise(assocs=n())


unique_proteins = unique(cvd$name)
