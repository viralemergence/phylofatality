## phylofatality
## 07_tables 
## danbeck@ou.edu, carolinecummings@ou.edu
## last update 2/4/2024

## clean environment & plots
rm(list=ls()) 
graphics.off()

## packages
library(dplyr)
library(gt)
library(gtsummary)
library(plyr)
library(tidyr)
library(tidyverse)
library(kableExtra)
library(knitr)

#load in PS data
setwd("~/Desktop/PCM Class/phylofatality/clean/csv files")
data=read.csv("PS_data.csv")

#lambda --> "\u03BB"

#make a table for phylogenetic signal 
data %>%
  summarize(vfamily, dataset, variable, lambda, lambda_lower, lambda_lower_p, 
            lambda_upper, lambda_upper_p) %>%
  arrange(dataset) %>%
  dplyr::rename("p-value (randomness)" = lambda_lower_p, "virus family"=vfamily,
                "\u03BB lower"=lambda_lower, "\u03BB upper"=lambda_upper,
                "p-value (BM)" = lambda_upper_p) %>%
  kable(align = "ccccclcl", digits = 3) %>%
  kable_styling(font_size = 15) %>%
  add_header_above(c("Pagel's \u03BB" = 8))


#make a table for the clades
setwd("~/Desktop/PCM Class/phylofatality/clean/csv files")
data=read.csv("pf_riskyclades.csv")

#filter to risky bats in mammal-wide analysis
bats<- data%>%
  slice(c(1,3,5,8,9,15,18,25,29))

#make table
bats%>%
  summarize(virus, var, taxa)%>%
  arrange(virus)%>%
  dplyr::rename("variable"=var)%>%
  kable(align = "ccl", digits = 3) %>%
  kable_styling(font_size = 15) %>%
  add_header_above(c("Risky bat clades: Mammal-wide analysis" =3))
