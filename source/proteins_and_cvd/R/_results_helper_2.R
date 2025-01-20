library(tidyverse)

setwd("~/Projects/python/proteins/results/attenuation_factor/hosp/agesex/40-69")

cvd = read_csv("attenuation_factor_composite_CVD_agesex.csv")
unique(cvd$id)
death = read_csv("attenuation_factor_death_agesex.csv")
bonf = 0.05/439


cvd = subset(cvd, p>bonf)
summary = cvd %>% group_by(formula) %>%
  summarise(assocs=n())


unique_proteins = unique(cvd$name)
