library("ggplot2")
library("tidyverse")
library("data.table")
library("tidytext")

# setwd("C:/Users/s1654019/Projects/python/proteins/results/")
# results = "_old/incremental_parallel_correct/hosp/agesex/merged_results.csv"
# deaths = "_old/incremental_parallel_deaths/hosp/agesex/merged_results.csv"

setwd("~/Projects/python/proteins/results/incremental_parallel/hosp/agesex/40-69")
ds = read_csv("merged_results.csv")
basic = "age+sex+protein"
full = paste0("age+sex+protein+avg_sys+Total_cholesterol+HDL_cholesterol+",
              "pack_years+rheum_arthritis_Y+diabetes_Y+years+rank+on_pill")

### how many proteins significant, altogether?
significant_basic = ds %>% 
  filter(formula == basic) %>%
  filter(p < 0.05/439)

basic_PH = significant_basic %>%
  filter(test_p < 0.05)

significant_basic = ds %>% 
  filter(formula == basic) %>%
  filter(p < 0.05/439) %>%
  filter(test_p > 0.05)

significant_basic_HR = significant_basic %>% arrange(desc(hr))

basic_duplicates = significant_basic %>% 
  mutate(name_duplicated = duplicated(name))


significant_full = ds %>% 
  filter(formula == full) %>%
  filter(p < 0.05/439) 

full_PH = significant_full %>%
  filter(test_p < 0.05)

significant_full = ds %>% 
  filter(formula == full) %>%
  filter(p < 0.05/439) %>%
  filter(test_p > 0.05)

full_duplicates = significant_full %>% 
  mutate(name_duplicated = duplicated(name))
table(full_duplicates$name_duplicated)

duplicates = subset(full_PH, name_duplicated)
count_duplicates = unique(duplicates$name)


significant_sum = ds %>% 
  filter(formula == full) %>%
  filter(p < 0.05/439) %>%
  group_by(event) %>% 
  summarise(sum = n())

## Load data
####################################################################
bonf = 0.05/439
# pheno = fread(results)
# pheno = pheno[,c(-1)]
# pheno = pheno %>% filter(event != "CVD_death") # remove CVD_death
# deaths = fread(deaths)
# ds = rbind(pheno, deaths)

# find only significant 
full_mod = paste0("age+sex+protein+avg_sys+Total_cholesterol+HDL_cholesterol",
                  "+pack_years+rheum_arthritis_Y+diabetes_Y+years+rank+on_pill")

# list of proteins that are significant in full models
composite_CVD = ds %>% filter(event == "composite_CVD") %>% filter(formula == full_mod)
sig_proteins = ds %>% filter(formula == full_mod) %>% filter(p < bonf) 
composite_CVD = sig_proteins %>% filter(event == "composite_CVD")
bad = composite_CVD %>% filter(hr > 1)
bad_grouped = bad %>% mutate(group_name = word(name, 1)) %>%
  group_by(group_name) %>%
  summarise(mean_hr = mean(hr),
            mean_sd = sd(hr),
            max_P = max(p))
bad_grouped_test = bad %>% mutate(group_name = word(name, 1)) 
good = composite_CVD %>% filter(hr < 1)

# HF #
hf = ds %>% filter(event == "hf") %>% filter(formula == full_mod)
sig_proteins = hf %>% filter(p < bonf) 

bad = sig_proteins %>% filter(hr > 1)
good = sig_proteins %>% filter(hr < 1)

good_grouped = sig_proteins %>% mutate(group_name = word(name, 1)) %>%
  group_by(group_name) %>%
  summarise(mean_hr = mean(hr),
            mean_sd = sd(hr),
            max_P = max(p))
bad_grouped_test = bad %>% mutate(group_name = word(name, 1)) 

# CVD_death #
CVD_death = ds %>% filter(event == "CVD_death") %>% filter(formula == full_mod)
sig_proteins = CVD_death %>% filter(p < bonf) 

bad = sig_proteins %>% filter(hr > 1)
good = sig_proteins %>% filter(hr < 1)

# death #
death = ds %>% filter(event == "death") %>% filter(formula == full_mod)
sig_proteins = death %>% filter(p < bonf) 

bad = sig_proteins %>% filter(hr > 1)
good = sig_proteins %>% filter(hr < 1)

good_grouped = sig_proteins %>% mutate(group_name = word(name, 1)) %>%
  group_by(group_name) %>%
  summarise(mean_hr = mean(hr),
            mean_sd = sd(hr),
            max_P = max(p))

bad_grouped_test = bad %>% mutate(group_name = word(name, 1)) 

# CVD_death #

chd = df %>% filter(event == "chd_nos") %>% 
  filter(formula == full_mod) %>% 
  mutate(p_adjusted = p.adjust(p, method = "BH"))
