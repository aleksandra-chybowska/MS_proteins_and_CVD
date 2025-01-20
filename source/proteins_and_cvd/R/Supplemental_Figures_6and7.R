library("ggplot2")
library("tidyverse")
library("data.table")
library("tidytext")

# setwd("C:/Users/s1654019/Projects/python/proteins/results/")
# results = "_old/incremental_parallel_correct/hosp/agesex/merged_results.csv"
# deaths = "_old/incremental_parallel_deaths/hosp/agesex/merged_results.csv"
setwd("~/Projects/python/proteins/results/incremental_parallel/hosp/agesex/40-69/")
ds = read_csv("merged_results.csv")


## Load data
####################################################################
bonf = 0.05/439
# pheno = fread(results)
# pheno = pheno[,c(-1)]
# pheno = pheno %>% filter(event != "CVD_death") # remove CVD_death
# deaths = fread(deaths)
# ds = rbind(pheno, deaths)

# find only significant 
basic_mod = "age+sex+protein"
full_mod = paste0("age+sex+protein+avg_sys+Total_cholesterol+HDL_cholesterol",
                  "+pack_years+rheum_arthritis_Y+diabetes_Y+years+rank+on_pill")

# list of proteins that are significant in full models
sig_proteins = ds %>% filter(formula == full_mod) %>% 
  group_by(event) %>%
  mutate(p_adjusted = p.adjust(p, method = "BH")) %>%
  mutate(is_bonf_significant = p < bonf) %>%
  mutate(is_bh_significant = p_adjusted < 0.05)

dim(sig_proteins %>% filter(event == "composite_CVD") %>% 
      filter(is_bh_significant == "TRUE"))

sig_proteins_005 = sig_proteins %>% filter(p < 0.05)
write.csv(sig_proteins_005, "full_significant_0_05.csv", row.names = F)

sig_proteins_basic = ds %>% filter(formula == basic_mod) %>%
  group_by(event) %>%
  mutate(p_adjusted = p.adjust(p, method = "BH")) %>%
  mutate(is_bonf_significant = p < bonf) %>%
  mutate(is_bh_significant = p_adjusted < 0.05)

table(sig_proteins_basic$event, sig_proteins_basic$p_adjusted < 0.05) # 4 of 47 dont meet coxph
sig_proteins_basic_005 = sig_proteins_basic %>% filter(p < 0.05)
write.csv(sig_proteins_basic_005, "basic_significant_0_05.csv", row.names = F)
