library(tidyverse)
library(tableone)
library(DataExplorer)

continous_normal_stats = function(variable, set) {
  mean = mean(set[[variable]])
  sd = sd(set[[variable]])
  ret = c(mean, sd)
  return(ret)
}

continous_non_normal_stats = function(variable, set) {
  median = median(set[[variable]])
  q1 = quantile(set[[variable]], 0.25)
  q3 = quantile(set[[variable]], 0.75)
  ret = c(median, q1, q3)
  return(ret)
}

categorical_stats = function(variable, response, set) {
  condition = paste0(variable, "==", response)
  tmp_set = subset(set, eval(parse(text=condition)))
  n_set = nrow(set)
  n_tmp_set = nrow(tmp_set)
  prop = n_tmp_set * 100/n_set
  ret = c(n_tmp_set, prop)
  return(ret)
}


#categorical = c("sex", "rheum_arthritis_Y", "diabetes_Y", "on_pill", "years")
#continous_normal = c("age", "avg_sys", "bmi", "HDL_cholesterol", "Total_cholesterol")
#continous_not_normal = c("pack_years", "rank")

gather_set_statistics = function(set) {
  set = cvd
  sex = categorical_stats("sex", 1, set)
  ra = categorical_stats("rheum_arthritis_Y", 1, set)
  diabetes = categorical_stats("diabetes_Y", 1, set)
  on_pill = categorical_stats("on_pill", 1, set)

  
  age = continous_normal_stats("age", set)
  avg_sys = continous_normal_stats("avg_sys", set)
  bmi = continous_normal_stats("bmi", set)
  total_cholesterol = continous_normal_stats("Total_cholesterol", set)
  
  hdl = continous_non_normal_stats("HDL_cholesterol", set)
  tte = continous_non_normal_stats("tte", set)
  rank = continous_non_normal_stats("rank", set)
  pack_years = continous_non_normal_stats("pack_years", set)
  years = continous_non_normal_stats("years", set)
  
  formatted_list = data.frame(
    "n" = nrow(set),
    "tte" = sprintf("%.1f\r\n[%.1f, %.1f]", tte[1], tte[2], tte[3]),
    "age" = sprintf("%.1f\r\n(%.1f)", age[1], age[2]),
    "sex"= sprintf("%.f\r\n(%.1f%%)", sex[1], sex[2]),
    "bmi" = sprintf("%.1f\r\n(%.1f)", bmi[1], bmi[2]),
    "simd" = sprintf("%.1f\r\n[%.1f, %.1f]", rank[1], rank[2], rank[3]),
    "education" = sprintf("%.1f\r\n[%.1f, %.1f]", years[1], years[2], years[3]),
    "on_pill" = sprintf("%.f\r\n(%.1f%%)", on_pill[1], on_pill[2]),
    "diabetes" = sprintf("%.f\r\n(%.1f%%)", diabetes[1], diabetes[2]),
    "ra" = sprintf("%.f\r\n(%.1f%%)", ra[1], ra[2]),
    "pack_years" = sprintf("%.1f\r\n[%.1f, %.1f]", pack_years[1], pack_years[2], pack_years[3]),
    "sys_bp" = sprintf("%.1f\r\n(%.1f)", avg_sys[1], avg_sys[2]),
    "total_cholesterol" = sprintf("%.1f\r\n(%.1f)", total_cholesterol[1], total_cholesterol[2]),
    "HDL_cholesterol" = sprintf("%.1f\r\n[%.1f, %.1f]", hdl[1], hdl[2], hdl[3])
  )
  
  return(formatted_list)
}

setwd("C:/Users/s1654019/Projects/python/proteins/results/cox/hosp/prepped/")

cvd = read_csv("cox_hosp_composite_CVD_prepped.csv")
create_report(cvd)
cvd$sex = ifelse(cvd$sex == 'M', 1, 0)
#cvd$years = as.factor(cvd$years)
statistics <- gather_set_statistics(cvd)
statisticsReport = t(bind_rows(statistics))
statisticsReport = as.data.frame(statisticsReport)
write.csv(statisticsReport, "report.csv")

###############################
#  events
###############################

ds = list(
  chd_nos = subset(read_csv("cox_hosp_chd_nos_prepped.csv"), event == 1),
  cvd = subset(read_csv("cox_hosp_composite_CVD_prepped.csv"), event==1),
  csv_death = subset(read_csv("cox_hosp_CVD_death_prepped.csv"), event==1),
  death = subset(read_csv("cox_hosp_death_prepped.csv"), event==1),
  hf = subset(read_csv("cox_hosp_hf_prepped.csv"), event==1),
  stroke = subset(read_csv("cox_hosp_isch_stroke_prepped.csv"), event==1),
  mi = subset(read_csv("cox_hosp_myocardial_infarction_prepped.csv"), event==1),
  tia = subset(read_csv("cox_hosp_tia_prepped.csv"), event==1)
)

get_name_n_tte_mean = function(set, name_set) {
  return(
    data.frame(
      "outcome" = name_set, 
      "n" = nrow(set), 
      "mean_tte" = mean(set[["tte"]]),
      "sd" = sd(set[["tte"]])
    )
  )
}
second_table = lapply(seq_along(ds), function(i) {get_name_n_tte_mean(ds[[i]], names(ds)[[i]])})
second_table = bind_rows(second_table)
write.csv(second_table, "report_pt2.csv")
## in the future, use tableOne (https://cran.r-project.org/web/packages/tableone/vignettes/introduction.html)