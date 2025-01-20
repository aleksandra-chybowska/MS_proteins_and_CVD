# Read csvs

setwd("~/Projects/python/proteins/results/attenuation_factor/hosp/agesex/40-69/")
temp = list.files(pattern="attenuation_factor_")
myfiles = lapply(temp, read.csv)
merged = rbindlist(myfiles)
merged = merged[,-1] # remove the ID col
write.csv(merged, "all_outcomes_attenuation_merged.csv", row.names = F)

cvd = merged %>% filter(event == "composite_CVD")
data.frame(unique = unique(cvd$name))

