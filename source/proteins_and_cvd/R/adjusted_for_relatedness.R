library(survival)
library(survminer)
library(kinship2)
library(coxme)
library(tidyverse)
library(data.table)
library(gridExtra)

extract_coxme_table <- function (mod){
  beta <- mod$coefficients #$fixed is not needed
  hr <-exp(beta)
  nvar <- length(beta)
  nfrail <- nrow(mod$var) - nvar
  se <- sqrt(diag(mod$var)[nfrail + 1:nvar])
  z<- round(beta/se, 2)
  p<- signif(1 - pchisq((beta/se)^2, 1), 2)
  lci = exp(beta-1*1.96*se)
  uci = exp(beta+1*1.96*se)
  table=data.frame(cbind(beta,hr,se,z,lci,uci,p))
  return(table)
}

# Rank Based Inverse Normalisation of the data
transform <- function(x) {
  transformed <- qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))
  return(transformed)
}

serialize_coxph_summary <- function(mod, protein) {
  coefs = extract_coxme_table(mod)
  hr <-coefs[3,2]
  p <-coefs[3,7]
  lci <-coefs[3,5]
  uci <-coefs[3,6]
  assumptions = cox.zph(mod)
  local=assumptions$table["protein", 'p']
  global=assumptions$table["GLOBAL", 'p']
  ret = data.frame(protein, hr, p, lci, uci, local, global)
  return(ret)
}

setwd("C:/Users/s1654019/Projects/python/proteins")
flag = "hosp" # hosp_gp, hosp, hosp_gp_cons
run = "agesex"
type = "40-69"
bonf = 0.05/439

path = paste0("results/incremental_parallel/", flag, "/", run, "/", type, "/")
plot_path = paste0("plots/incremental_parallel/", flag, "/", run, 
                   "/", type, "/")

dir.create(file.path(plot_path))
ds = read_csv(paste0(path, "merged_results.csv"))

full_mod = paste0("age+sex+protein+avg_sys+Total_cholesterol+HDL_cholesterol",
                  "+pack_years+rheum_arthritis_Y+diabetes_Y+years+rank+on_pill")
sig_proteins = ds %>% filter(formula == full_mod) %>% filter(p < bonf) 

proteins = read_csv('results/cox/40-69/proteins_hosp_all_events_scaled_8343.csv')
cols = c('id', sig_proteins$id)
proteins = proteins[colnames(proteins) %in% cols] ## 38 unique names

ped = read_csv("data/phenotypes/2022-01-17_pedigree.csv")

# Addressing "Value of 'momid' not found in the id list 4091/103027"
missing_mom = data.frame(famid=4091, volid=103027, father=0, mother=0, sex="F", age=NA)
ped = rbind(ped, missing_mom)
kin = with(ped, pedigree(volid, father, mother, sex, famid=famid))
kin_model = kinship(kin)

results = data.frame()

for (i in 1:nrow(sig_proteins)) {
  row = sig_proteins[i, ]
  
  protein_name = row$id
  protein_df = proteins[c("id", protein_name)]
  colnames(protein_df) = c("id", "protein")
  
  event = row$event
  cox_path = paste0("results/cox/", type, "/cox_", flag, "_", event, "_prepped.csv")
  cox = read_csv(cox_path)
  
  cox = merge(cox, protein_df, by="id")
  
  model = coxme(Surv(tte, event) ~
                      age + 
                      sex +
                      protein + 
                      avg_sys + 
                      Total_cholesterol + 
                      HDL_cholesterol +
                      pack_years + 
                      rheum_arthritis_Y + 
                      diabetes_Y + 
                      years + 
                      rank + 
                      on_pill +
                      (1|id), data = cox, varlist = kin_model*2)
  
  
  # same model for diagnostics
  res.cox = coxph(Surv(tte, event) ~
                  age + 
                  sex +
                  protein + 
                  avg_sys + 
                  Total_cholesterol + 
                  HDL_cholesterol +
                  pack_years + 
                  rheum_arthritis_Y + 
                  diabetes_Y + 
                  years + 
                  rank + 
                  on_pill, data = cox)
  
  P1 <- ggcoxdiagnostics(res.cox, type = "deviance",
                   linear.predictions = FALSE, ggtheme = theme_bw())
  P2 <- ggcoxfunctional(Surv(tte, event) ~ age, data = cox)
  P3 <- ggcoxfunctional(Surv(tte, event) ~ protein, data = cox)
  P4 <- ggcoxfunctional(Surv(tte, event) ~ avg_sys, data = cox)
  P5 <- ggcoxfunctional(Surv(tte, event) ~ Total_cholesterol + sqrt(Total_cholesterol) + log(Total_cholesterol), data = cox)
  P6 <- ggcoxfunctional(Surv(tte, event) ~ HDL_cholesterol, data = cox)
  P7 <- ggcoxfunctional(Surv(tte, event) ~ pack_years, data = cox)
  
  test.ph <- cox.zph(res.cox)
  
  pdf(paste0(plot_path, protein_name, "_", event, "_diagnostics.pdf"), paper="a4r")

  print(P1 + ggtitle(paste0("Cox residuals for ", row$name, "\n(", protein_name, ")")))
  
  print(grid.arrange(P2[[1]], P3[[1]], P4[[1]], P5[[1]], P5[[2]], P5[[3]], P6[[1]], P7[[1]],
               top = "Martingale Residuals \nof Null Cox Model", nrow = 4, ncol = 2))
  
  print(ggcoxzph(test.ph, newpage = FALSE, point.size=0.4, font.main = 8,
           font.submain =8,
           font.caption = 8,
           font.x = 8,
           font.y = 8, 
           font.tickslab = 8,
           font.legend = 8))

  dev.off()
  
  summary(model)
  to_add = serialize_coxph_summary(model, protein_name)
  results = rbind(results, to_add)
}

cresults = cbind(sig_proteins, results)

results_missing_added = cresults[c(1:2, 9, 12:17)]
comparison = cresults[c(1:2, 9, 3, 12, 4, 13, 7, 16)]

write.csv(results_missing_added, paste0(path, "significant_coxme.csv"), row.names = F)
write.csv(comparison, paste0(path, "comparison_cox_coxme.csv"), row.names = F)
