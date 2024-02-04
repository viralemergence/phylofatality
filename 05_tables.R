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
