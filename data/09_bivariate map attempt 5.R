#Hanta/Arena stuff...bivariate attempt 5

#load in packages
library(readr)
library(tidyverse)
library(dbplyr)
library(biscale)

## clean environment & plots
rm(list=ls()) 
graphics.off()

#1 hanta data
#hanta_ids <- read_csv("~/Documents/Hanta/Data/Hanta_withHumans.csv") %>% 
 # drop_na(Country)
#country_genomes <- hanta_ids %>% 
 # group_by(Country) %>% 
  #summarize(Number = n()) %>% 
  #arrange(desc(Number))
#colnames(country_genomes) <- c('region', "Number_Hanta")

#1 CoV clade data
setwd("~/Desktop/GitHub/phylofatality/csv files")
species <- read_csv("pf_riskyspecies.csv")
#subset data 
cov_mean<-species %>% 
  filter(host=="bat", var=="mean", virus=="cov")
cov_mean$species<- gsub("_", " ", cov_mean$species)


#2 Arena data
#diversity_genomes <- hanta_ids %>% 
 # group_by(Country) %>% 
  #summarize(Number = length(unique(Species))) %>% 
  #arrange(desc(Number))
#colnames(diversity_genomes) <- c('region', 'Number_species')
#world <- map_data("world")

#2 bat data
setwd("~/Desktop/GitHub/phylofatality/data")
bats=readRDS("bat shp.rds")
bats <- bats[bats$binomial %in% cov_mean$species,] 
world <- map_data("world")

#3 Categorize Hanta genome counts using computed Jenks natural breaks
breaks <- BAMMtools::getJenksBreaks(var = country_genomes$Number_Hanta, k = 4)

genomes_data <- country_genomes %>% 
  mutate(cat_hanta = cut(Number_Hanta, breaks = breaks,
                         include.lowest = T)) 

#4 Categorize Arena genome counts using computed Jenks natural breaks
diversity_breaks <- BAMMtools::getJenksBreaks(var = diversity_genomes$Number_species, k = 4)

#Add Arena data
genomes_data2 <- diversity_genomes %>% 
  mutate(cat_arena = cut(Number_species, breaks = diversity_breaks,
                         include.lowest = T))




#Hanta_arena
hanta_arena <- full_join(genomes_data, genomes_data2, by="region") %>% 
  mutate(Number_Arena = coalesce(Number_species,0),
         Number_Hanta = coalesce(Number_Hanta,0))
colnames(hanta_arena) <- c('region','number_hanta','cat_hanta','number_species','num_spc2','cat_species')

arena_breaks <- BAMMtools::getJenksBreaks(var = hanta_arena$number_arena, k = 7)

data <- bi_class(hanta_arena, x=number_hanta, y=number_species, style = "jenks", dim = 4)


#pick color palette
bivariate_scale <- c("3 - 3" = "#3F2949", 
                     "2 - 3" = "#435786",
                     "1 - 3" = "#4885C1", 
                     "3 - 2" = "#77324C",
                     "2 - 2" = "#806A8A", 
                     "1 - 2" = "#89A1C8",
                     "3 - 1" = "#AE3A4E", 
                     "2 - 1" = "#BC7C8F",
                     "1 - 1" = "#CABED0" )


#CREATE BIVARIATE CLASSES
library(biscale)


map_data <- full_join(world, data, by="region")




hanta_number_of_genomes <- ggplot() +
  geom_polygon(data = map_data, aes(x=long, y=lat, group = group, fill = bi_class),color = "grey60",
               linewidth = 0.05, alpha = 1) +
  bi_scale_fill(pal = "GrPink2", dim = 4, rotate_pal = F, na.value="grey95", flip_axes = F) +
  geom_segment(aes(x=-18, xend=-18, y = -60, yend=90), linetype = "dashed")+
  labs(x=NULL,
       y=NULL,
       fill = NULL)+
  bi_theme(base_family = "Times New Roman")+
  theme(panel.background = element_rect(fill = "white",
                                        color = "grey40"),
        legend.position = "none",
        legend.key.width = unit(3, 'cm'),
        legend.spacing.x = unit(0,'cm'),
        text = element_text(size = 20, family = "Times New Roman"))

breaks_vals <- bi_class_breaks(hanta_arena, x = number_hanta, y = number_species, dim = 4, style = "jenks", clean_levels = T, dig_lab = 4)

legend <- bi_legend(pal = "GrPink2",
                    dim =4,
                    xlab = "No. Hantavirus genomes",
                    ylab = "No. Hantavirus Species",
                    size = 6,
                    breaks = breaks_vals,
                    arrows = F,
                    flip_axes = F,
                    rotate_pal = F)


final_map <- hanta_number_of_genomes +
  annotate("segment", x=-30, xend = -70, y=-62, yend = -62, arrow=arrow(length = unit(0.2, 'cm')))+
  annotate("segment", x=-6, xend = 34, y=-62, yend = -62, arrow=arrow(length = unit(0.2, 'cm')))+
  annotate("text", x=-50, y=-57,label="New World Viruses")+
  annotate("text", x=14, y=-57,label="Old World Viruses")

library(cowplot)

ggdraw() +
  draw_plot(final_map, 0, 0, 1, 1) +
  draw_plot(legend, 0.08, .20, 0.2, 0.2)

ggsave("~/Documents/Hanta/Plots/NumberofGenomesMap2.jpg", final_map, dpi = 1500, height = 10, width = 18.08, units = "in")
#save TIFF

library(ragg)
agg_tiff("~/Documents/Hanta/Plots/NumberofGenomesMap2.tiff", res = 300, height = 10, width = 18.08, units = "in")

dev.off()