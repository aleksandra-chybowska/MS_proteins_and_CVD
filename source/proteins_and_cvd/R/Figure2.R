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
full_mod = paste0("age+sex+protein+avg_sys+Total_cholesterol+HDL_cholesterol",
                  "+pack_years+rheum_arthritis_Y+diabetes_Y+years+rank+on_pill")

# list of proteins that are significant in full models
sig_proteins = ds %>% filter(formula == full_mod) %>% filter(p < bonf) 
table(sig_proteins$test_p < 0.05) # 4 of 47 dont meet coxph


heart = c("hf", 
          "composite_CVD")
death = c("death", "CVD_death")
brain = c("tia", "isch_stroke")

sig_heart = sig_proteins %>% filter(event %in% heart)
sig_death = sig_proteins %>% filter(event %in% death)

plotting_heart = ds %>% 
  filter(formula == full_mod) %>% 
  filter(event %in% heart) %>% 
  filter(id %in% sig_heart$id) %>%
  mutate(group = "Cardiovascular Diseases")

# order by HRs for heart failure
ordering_heart = plotting_heart %>%
  filter(event == "composite_CVD") %>%
  select(name, hr) %>%
  rename(ob = hr)

plotting_heart = left_join(plotting_heart, ordering_heart, by="name")

plotting_death = ds %>% 
  filter(formula == full_mod) %>% 
  filter(event %in% death) %>%
  filter(id %in% sig_death$id) %>%
  mutate(group = "Death")

# order by HRs for death
ordering_death = plotting_death %>%
  filter(event == "CVD_death") %>%
  select(name, hr) %>%
  rename(ob = hr)

plotting_death = left_join(plotting_death, ordering_death, by="name")

plotting = rbind(plotting_heart, plotting_death) 
# plotting = plotting %>% 
#   mutate(shape = ifelse(p<bonf, "triangle", "circle"))

plotting = plotting %>% 
  mutate(shape = as.factor(case_when(
    p < bonf & test_p >= 0.05 ~ 17, #"full_triangle",
    p < bonf & test_p < 0.05 ~ 2, # "open_triangle"
    p >= bonf & test_p >= 0.05 ~ 19, # "full_circle",
    p >= bonf & test_p < 0.05 ~ 1 # "open_circle",
  ))
)

## TODO: adjusting for relatedness fix

## TODO: remove items that do not meet cox assumptions
dim(subset(plotting, test_p<0.05)) # 10

#plotting = filter(plotting, plotting$test_p>0.05)

## Prep a plotting ds
####################################################################

dict = c("composite_CVD" = "Composite CVD",
         "CVD_death" = "CVD death",
         "death" = "Death",
         "hf" = "Heart Failure")

out_cont = data.frame(Outcome=plotting$name,
                      Predictor=dict[plotting$event],
                      Beta=plotting$hr,
                      P=plotting$p, 
                      LCI=plotting$lci,
                      UCI=plotting$uci,
                      stringsAsFactors=FALSE,
                      Shape=plotting$shape,
                      Group=plotting$group, 
                      OB = plotting$ob)

out_cont = out_cont %>% 
  mutate(Outcome = reorder_within(Outcome, OB, Group))


# out_cont_heart = out_cont %>% filter(Group == "Cardiovascular Diseases")
# out_cont_death = out_cont %>% filter(Group == "Death")

## Prep plot
####################################################################


size_text = 12
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

color_scale <- c(
  "Heart Failure" = "#c77cff",
  "Composite CVD" = "#f8766d",
  "Death" = "#00bfc4",
  "CVD death" = "#f3b772")

color_palette <- scale_color_manual(values = color_scale)

stacked = ggplot(out_cont, 
                 aes(y=Beta, x=Outcome, group=Predictor, colour=Predictor)) + 
  geom_point(aes(shape=Shape), size = 2, position = position_dodge(0.5))+
  geom_errorbar(aes(ymin = LCI, ymax = UCI),
                position = position_dodge(0.5), width = 0.1)+
  ylab("Hazard Ratio per SD")+ 
  xlab ("") +
  theme_bw() +
  geom_hline(yintercept = 1, linetype = "dotted")+
  coord_flip() + scale_x_reordered() + 
  scale_shape_manual(values = c(1, 2, 17, 19)) +
  My_Theme + 
  facet_wrap(~Group, ncol = 1, scales = "free_y") +
  guides(shape = "none") + color_palette

stacked
# stop()
library(grid)
gt = ggplot_gtable(ggplot_build(stacked))
gt$heights[10] = 0.5*gt$heights[10]
grid.draw(gt)

setwd("~/Desktop")
png(paste0("Stacked_Forest_plot_heart_mapped.png"), height = 880, width=720)
grid.draw(gt)
dev.off()

# 
# color_scale <- c(
#   "hf" = "#c77cff",
#   "composite_CVD" = "#f8766d",
#   "death" = "#00bfc4",
#   "CVD_death" = "#f3b772")
# 
# color_palette <- scale_color_manual(values = color_scale, 
#                                     limits = names(color_scale))
# 
# g <- guides(fill = guide_legend(override.aes = list(color = color_palette,
#                                                     shape = c(16, 16, 16), 
#                                                     size = c(1, 1, 1),
#                                                     alpha = c(1, 1, 1)),))
# 
# stacked_heart = ggplot(out_cont_heart, 
#                  aes(y=Beta, x=Outcome, group=Predictor, colour=Predictor)) + 
#   geom_point(aes(shape=Shape), size = 2, position = position_dodge(0.5), 
#              key_glyph = "point") + g +
#   geom_errorbar(aes(ymin = LCI, ymax = UCI),
#                 position = position_dodge(0.5), width = 0.1)+
#   ylab("Hazard Ratio per SD")+ 
#   xlab ("") +
#   theme_bw() +
#   geom_hline(yintercept = 1, linetype = "dotted")+
#   coord_flip() + My_Theme + 
#   facet_wrap(~Group, ncol = 1, scales = "free_y") +
#   guides(shape = "none")
# 
# stacked_death = ggplot(out_cont_death, 
#                        aes(y=Beta, x=Outcome, group=Predictor, colour=Predictor)) + 
#   geom_point(aes(shape=Shape), size = 2, position = position_dodge(0.5), 
#              key_glyph = "point") + g +
#   geom_errorbar(aes(ymin = LCI, ymax = UCI),
#                 position = position_dodge(0.5), width = 0.1)+
#   ylab("Hazard Ratio per SD")+ 
#   xlab ("") +
#   theme_bw() +
#   geom_hline(yintercept = 1, linetype = "dotted")+
#   coord_flip() + My_Theme + 
#   facet_wrap(~Group, ncol = 1, scales = "free_y") +
#   guides(shape = "none")
# 
# groups = unique(out_cont$Predictor)
# 
# combined_plot = stacked_heart + stacked_death + 
#   plot_layout(heights = c(1, 2), guides = "collect") & 
#   theme(legend.position = "bottom") & scale_color_manual(values = color_scale, 
#                                                          limits = names(color_scale))
# 


