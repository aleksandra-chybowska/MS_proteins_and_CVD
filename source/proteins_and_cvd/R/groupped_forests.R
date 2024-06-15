library("ggplot2")
library("tidyverse")
library("data.table")

setwd("C:/Users/s1654019/Projects/python/proteins/results/")
results = "_old/incremental_parallel_correct/hosp/agesex/merged_results.csv"
deaths = "_old/incremental_parallel_deaths/hosp/agesex/merged_results.csv"
# path = "/Cluster_Filespace/Marioni_Group/Ola/Smoking/BayesR_EpiScore/data/"
# neuro_path = "/Cluster_Filespace/Marioni_Group/Ola/Smoking/Neuroimaging/"

transform = function(x) {
  transformed = qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))
  return(transformed)
}

## Load brain data
####################################################################
bonf = 0.05/439
pheno = fread(results)
pheno = pheno[,c(-1)]
pheno = pheno %>% filter(event != "CVD_death") # remove CVD_death
deaths = fread(deaths)

ds = rbind(pheno, deaths)

# find only significant 
full_mod = paste0("age+sex+protein+avg_sys+Total_cholesterol+HDL_cholesterol",
                  "+pack_years+rheum_arthritis_Y+diabetes_Y+years+rank+on_pill")
sig_proteins = ds %>% filter(formula == full_mod) %>% filter(p < bonf) # list of proteins that are significant in full models
ids = sig_proteins$id

plotting = ds %>% filter(formula == full_mod) %>% filter(id %in% ids)

table(plotting$event)
heart = c("chd_nos", "myocardial_infarction", "hf")
brain = c("tia", "isch_stroke")
heart_without_hf = c("chd_nos", "myocardial_infarction", "composite_CVD")
death = c("death", "CVD_death")

## Variable distribution check, remove NAs, scale
####################################################################
plotting = ds %>% filter(formula == full_mod) %>% filter(id %in% ids)
plotting = plotting %>% filter(event %in% heart)
out_cont = data.frame(Outcome=plotting$name, 
                      Predictor=plotting$event,
                      Beta=plotting$hr,
                      P=plotting$p, 
                      LCI=plotting$lci,
                      UCI=plotting$uci,
                      stringsAsFactors=FALSE)

## Prep plot
####################################################################

size_text = 10
My_Theme = theme(
  panel.border = element_rect(colour="black",size=1, fill = NA),
  axis.title.x = element_text(size = size_text), # controls HR label size 
  axis.text.x = element_text(size = size_text),
  axis.text.y = element_text(size = size_text),
  axis.title.y = element_text(size = size_text),
  strip.text = element_text(size = size_text, face = "bold"),
  legend.text=element_text(size=size_text),
  legend.position="bottom",
  legend.background = element_rect(size=0.5, linetype="solid", colour="black"),
  legend.title=element_blank(),
  axis.title=element_text(size=size_text))

stacked = ggplot(out_cont,aes(y=Beta, x=Outcome, group=Predictor, colour=Predictor)) + 
  geom_point(size = 2, position = position_dodge(0.5))+
  geom_errorbar(aes(ymin = LCI, ymax = UCI),
                position = position_dodge(0.5), width = 0.1)+
  ylab("Hazard Ratio per SD")+ 
  xlab ("") +
  theme_bw() +
  geom_hline(yintercept = 1, linetype = "dotted")+
  theme(axis.text.x = element_text(size = 8, vjust = 0.5), axis.text.y = element_text(size = 8), legend.position = "right",
        plot.title = element_text(size = 8))+ theme(legend.title = element_text(hjust = 0.5)) +
  coord_flip() + My_Theme


pdf(paste0("Forest_plot_heart.pdf"), height = 30, width = 10)
stacked
dev.off()
