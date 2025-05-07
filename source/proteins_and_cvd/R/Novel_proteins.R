library(tidyverse)
library(data.table)

ot_path = "~/Projects/python/proteins/data/disease/open_targets/"
wd_path = "~/Projects/python/proteins/results/incremental_parallel/hosp/agesex/40-69"
pattern = "ot_data"

files = list.files(ot_path, pattern=pattern, full.names = T)
data_ot = files %>% lapply(read.csv)
names(data_ot) = c("chd", "composite_CVD", "death", 
                   "hf", "isch_stroke", "mi", "tia")

setwd(wd_path)
merged = read.csv("merged_results_bonf_significant_full.csv")
ds = merged[c("id", "name", "event", "covar")]
ds = ds %>% separate_rows(id, sep = "\\.")  # Splits IDs at each "."
ds_split = split(ds, ds$event)

matched_files = list(list(data_ot$composite_CVD, ds_split$composite_CVD),
                     list(data_ot$composite_CVD, ds_split$CVD_death),
                     list(data_ot$death, ds_split$death),
                     list(data_ot$hf, ds_split$hf))


process_match <- function(ot_data, ds_data) {
  if (!"UniProt_ID" %in% colnames(ot_data)) {
    stop("Column 'UniProt_ID' not found in Open Targets dataset")
  }
  
  if (!is.null(ds_data)) {
    ds_data$known = ds_data$id %in% ot_data$UniProt_ID
  }
  
  return(ds_data)
}

# Process all matched datasets in fully significant dataset
processed_data = map2(matched_files, seq_along(matched_files), ~{
  process_match(.x[[1]], .x[[2]])
})

all = rbindlist(processed_data)
all = all %>% filter(
  !id %in% c("A0A0A0MRD9", "A0A087WU08", "A0A024R6I7", "A0A087WWT3") # deleted proteins
)

all_groupped = all %>% group_by(event, covar) %>%
  summarise(
    all_known = all(known),  # AND operation: If any FALSE exists, result is FALSE
    unknown = paste(id[!known], collapse = ",")  # Concatenates IDs where known == FALSE
  ) %>%
  ungroup()

# modify Supplemental Table 7 based on these results
to_modify = read.csv("comparison_cox_coxme.csv")

merged_to_modify = merge(to_modify, all_groupped, 
                         by.x = c("id", "event"),
                         by.y=c("covar", "event"), all.x=TRUE) %>% arrange(event, desc(hr))
write.csv(merged_to_modify, "comparison_cox_coxme_with_novelty.csv", 
          row.names = FALSE)

# =================

setwd("../../../../attenuation_factor/hosp/agesex/40-69/")
cvd = read.csv("attenuation_factor_composite_CVD_agesex.csv") 
ds = cvd %>% select("id", "name", "covar", "formula") %>% distinct()
ds = ds %>% separate_rows(id, sep = "\\.")
ds$known = ds$id %in% data_ot$composite_CVD$UniProt_ID
ds = ds %>% filter(
  !id %in% c("A0A087WWT3", "A0A087X0I2", "A0A087X1J7", "A0A0U1RR20") # deleted proteins
)

ds_groupped = ds %>% group_by(covar) %>%
  summarise(
    all_known = all(known),  # AND operation: If any FALSE exists, result is FALSE
    unknown = paste(id[!known], collapse = ",")  # Concatenates IDs where known == FALSE
  ) %>%
  ungroup()

ds_groupped$unknown <- sapply(ds_groupped$unknown, function(x) {
  unique_proteins <- unique(strsplit(x, ",")[[1]])  # Split by comma and get unique values
  paste(unique_proteins, collapse = ",")  # Join unique values back into a string
})


merged_attenuation = merge(cvd, ds_groupped, by="covar", all.x=TRUE) %>% arrange(X)
merged_attenuation = merged_attenuation[-c(1, 2)]
write.csv(merged_attenuation, "Supplemental_Table_8_CVD_attenuation.csv", 
          row.names = FALSE)
