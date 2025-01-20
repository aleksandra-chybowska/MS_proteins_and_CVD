library(tidyverse)
library(tableone)
library(DataExplorer)

setwd("~/Projects/python/proteins/results/cox/40-69")

cvd = read_csv("cox_hosp_composite_CVD_prepped.csv")
create_report(cvd)
cvd$sex = ifelse(cvd$sex == 'M', 1, 0)
cvd$pack_years = exp(cvd$pack_years) - 1# inverse log+1 transform
cvd$bmi = exp(cvd$bmi)


categorical = c("sex", "rheum_arthritis_Y", "diabetes_Y")
continuous_normal = c("age", "avg_sys", "bmi", "Total_cholesterol")
continuous_non_normal = c("HDL_cholesterol", "tte", "rank", "pack_years", "years")
all = c("age", "sex", "bmi", "rank", "diabetes_Y", 
        "rheum_arthritis_Y", "pack_years", "avg_sys", "Total_cholesterol", 
        "HDL_cholesterol", "years")

tab <- CreateTableOne(vars = all,
                       factorVars = categorical, strata = "event" , data = cvd)
print(tab, nonnormal = continuous_non_normal, formatOptions = list(big.mark = ","))
tab3Mat <- print(tab, nonnormal = continuous_non_normal, 
                 quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
## Save to a CSV file
write.csv(tab3Mat, file = "myTable.csv")

### Ooopsie, we have got males in the group using contraception. Lets fix it.
### Ugh, on pill problematic. I am removing it from this table.

females = subset(cvd, sex == 0)
contingency_table <- table(females$on_pill, females$event)
chi_squared_test <- chisq.test(contingency_table)
chi_squared_test



###############################
#  events
###############################

ds = list(
  composite_CVD = subset(read_csv("cox_hosp_composite_CVD_prepped.csv"), event==1),
  CVD_death = subset(read_csv("cox_hosp_CVD_death_prepped.csv"), event==1),
  death = subset(read_csv("cox_hosp_death_prepped.csv"), event==1),
  chd_nos = subset(read_csv("cox_hosp_chd_nos_prepped.csv"), event == 1),
  myocardial_infarction = subset(read_csv("cox_hosp_myocardial_infarction_prepped.csv"), event==1),
  hf = subset(read_csv("cox_hosp_hf_prepped.csv"), event==1),
  isch_stroke = subset(read_csv("cox_hosp_isch_stroke_prepped.csv"), event==1),
  tia = subset(read_csv("cox_hosp_tia_prepped.csv"), event==1)
)

setwd("~//Projects/python/proteins/results/incremental_parallel/hosp/agesex/40-69")
events = read_csv("merged_results.csv")
basic = "age+sex+protein"
full = paste0("age+sex+protein+avg_sys+Total_cholesterol+HDL_cholesterol+",
              "pack_years+rheum_arthritis_Y+diabetes_Y+years+rank+on_pill")

### how many proteins significant, altogether?

significant_sum = events %>% 
  filter(formula == full) %>%
  filter(p < 0.05/439) %>%
  group_by(event) %>% 
  summarise(sum = n())

significant = events %>% 
  filter(formula == full) %>%
  filter(p < 0.05/439) %>%
  mutate(name_duplicated = duplicated(name))

duplicates = subset(significant, name_duplicated)
count_duplicates = unique(duplicates$name)

summary1 = events %>% 
  group_by(event) %>% 
  summarise(P_nominal_basic   = sum(formula == basic & p<0.05),
            P_nominal_full    = sum(formula == full & p<0.05), 
            ) %>%
  arrange(factor(event, levels = names(ds)))

summary2 = events %>% 
  filter(formula == full) %>%
  group_by(event) %>% 
  mutate(p_adjusted = p.adjust(p, method = "BH")) %>%
  summarise(
    P_FDR_full = sum(formula == full & p_adjusted < 0.05),
    P_bonferroni_full = sum(formula == full & p<(0.05/439))
  ) %>%
  arrange(factor(event, levels = names(ds)))

summary = left_join(summary1, summary2, by="event")

format_name_n_tte_mean = function(set, name_set) {
  return(
    data.frame(
      "event" = name_set,
      "outcome" = paste0(name_set, " (n=", nrow(set),")"),
      "mean_tte" = sprintf("%.1f (%.1f)", mean(set[["tte"]]), sd(set[["tte"]]))
    )
  )
}

second_table = lapply(seq_along(ds), 
                      function(i) {
                        format_name_n_tte_mean(ds[[i]], names(ds)[[i]])
                        })
second_table = bind_rows(second_table)

report_pt2 = left_join(second_table, summary, by="event")

write.csv(report_pt2, "report_pt2.csv", row.names = F)
