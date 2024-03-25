## phylofatality
## 05a_tables: PS and phylofactor 
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
library(tools)
library(kableExtra)
library(knitr)
library(readr)
library(forcats)

#load in PS data
setwd("~/Desktop/PCM Class/phylofatality/clean/csv files")
data=read.csv("PS_data.csv")

#note: lambda --> "\u03BB"
test<- data %>% 
  #filter(variable=="meanCFR", dataset=="all mammals")
  summarize(vfamily, variable, lambda, lambda_lower, lambda_lower_p, 
            lambda_upper, lambda_upper_p)%>%
  dplyr::rename("p-value (randomness)" = lambda_lower_p,
                "\u03BB lower"=lambda_lower, "\u03BB upper"=lambda_upper,
                "p-value (BM)" = lambda_upper_p)

test %>% kable(test, 
              caption = "Deaths from unintentional injuries in Scotland",
              format = "latex", booktabs = TRUE) %>%
  kable_styling(font_size = 10) %>%
  pack_rows(tab_kable, colnum = 1,
            index = table(fct_inorder(test$vfamily), useNA = "no"))

view(fig_2)

#tables for MeanCFR
test<- data %>% 
  #filter(variable=="meanCFR", dataset=="all mammals")
  summarize(vfamily, variable, lambda, lambda_lower, lambda_lower_p, 
            lambda_upper, lambda_upper_p)%>%
  dplyr::rename("p-value (randomness)" = lambda_lower_p, "virus family"=vfamily,
                "\u03BB lower"=lambda_lower, "\u03BB upper"=lambda_upper,
                "p-value (BM)" = lambda_upper_p) %>%
  kable(align = "ccclcl", digits = 3) %>%
  kable_styling(font_size = 15) %>%
  add_header_above(c("Pagel's \u03BB across mammals" = 6))

data %>% filter(variable=="meanCFR", dataset=="bats only")%>%
  summarize(vfamily, lambda, lambda_lower, lambda_lower_p, 
            lambda_upper, lambda_upper_p) %>%
  dplyr::rename("p-value (randomness)" = lambda_lower_p, "virus family"=vfamily,
                "\u03BB lower"=lambda_lower, "\u03BB upper"=lambda_upper,
                "p-value (BM)" = lambda_upper_p) %>%
  kable(align = "ccclcl", digits = 3) %>%
  kable_styling(font_size = 15) %>%
  add_header_above(c("Pagel's \u03BB only bats" = 6))

##table for Max CFR
data %>% filter(variable=="maxCFR", dataset=="all mammals")%>%
  summarize(vfamily, lambda, lambda_lower, lambda_lower_p, 
            lambda_upper, lambda_upper_p) %>%
  dplyr::rename("p-value (randomness)" = lambda_lower_p, "virus family"=vfamily,
                "\u03BB lower"=lambda_lower, "\u03BB upper"=lambda_upper,
                "p-value (BM)" = lambda_upper_p) %>%
  kable(align = "ccclcl", digits = 3) %>%
  kable_styling(font_size = 15) %>%
  add_header_above(c("Pagel's \u03BB across mammals" = 6))

data %>% filter(variable=="maxCFR", dataset=="bats only")%>%
  summarize(vfamily, lambda, lambda_lower, lambda_lower_p, 
            lambda_upper, lambda_upper_p) %>%
  dplyr::rename("p-value (randomness)" = lambda_lower_p, "virus family"=vfamily,
                "\u03BB lower"=lambda_lower, "\u03BB upper"=lambda_upper,
                "p-value (BM)" = lambda_upper_p) %>%
  kable(align = "ccclcl", digits = 3) %>%
  kable_styling(font_size = 15) %>%
  add_header_above(c("Pagel's \u03BB only bats" = 6))


##table for OT
data %>% filter(variable=="on.frac", dataset=="all mammals")%>%
  summarize(vfamily, lambda, lambda_lower, lambda_lower_p, 
            lambda_upper, lambda_upper_p) %>%
  dplyr::rename("p-value (randomness)" = lambda_lower_p, "virus family"=vfamily,
                "\u03BB lower"=lambda_lower, "\u03BB upper"=lambda_upper,
                "p-value (BM)" = lambda_upper_p) %>%
  kable(align = "ccclcl", digits = 3) %>%
  kable_styling(font_size = 15) %>%
  add_header_above(c("Pagel's \u03BB across mammals" = 6))

data %>% filter(variable=="on.frac", dataset=="bats only")%>%
  summarize(vfamily, lambda, lambda_lower, lambda_lower_p, 
            lambda_upper, lambda_upper_p) %>%
  dplyr::rename("p-value (randomness)" = lambda_lower_p, "virus family"=vfamily,
                "\u03BB lower"=lambda_lower, "\u03BB upper"=lambda_upper,
                "p-value (BM)" = lambda_upper_p) %>%
  kable(align = "ccclcl", digits = 3) %>%
  kable_styling(font_size = 15) %>%
  add_header_above(c("Pagel's \u03BB only bats" = 6))


#tables for phylofactor

#make a table for the clades
setwd("~/Desktop/PCM Class/phylofatality/clean/csv files")
bats=read.csv("pf_riskyclades.csv")

#filter to risky bats in mammal-wide analysis
bats<- bats%>%
  slice(c(1,3,5,8,9,15,18,25,29))
  
# Function to capitalize the first letter of each word
capitalize_first_letter <- function(text) {
  str_to_title(text)
}

# Apply the function to the taxa
bats$taxa <- sapply(bats$taxa, capitalize_first_letter) 

#make table
bats%>%
  summarize(virus, var, taxa)%>%
  arrange(virus)%>%
  dplyr::rename("variable"=var)%>%
  kable(align = "ccl", digits = 3) %>%
  kable_styling(font_size = 15) %>%
  add_header_above(c("Risky bat clades: Mammal-wide analysis" =3))


library("knitr")
library("kableExtra")
library("dplyr")
library("readr")
library("tidyr")
library("forcats")
library(vroom)

orig_ui_deaths = read_csv("https://www.opendata.nhs.scot/dataset/b0135993-3d8a-4f3b-afcf-e01f4d52137c/resource/89807e07-fc5f-4b5e-a077-e4cf59491139/download/ui_deaths_2020.csv")

summary_tab <- orig_ui_deaths %>%
  # Apply filters in same was as in plots so looking at one year, whole of Scotland, and exlucing "All" entries
  filter(Year == "2018", HBR == "S92000003", AgeGroup != "All" & Sex != "All", InjuryType != "Accidental exposure" & InjuryType != "All") %>%
  # Group the table according to injury type, age, and sex
  group_by(InjuryType, AgeGroup, Sex) %>%
  # Create a summary of total number of deaths for neater appearance to the table and demonstration of the figures
  summarise(total_deaths = sum(NumberOfDeaths)) %>%
  # Change the orientation of the table so that age groups become the variable headings
  pivot_wider(names_from = AgeGroup, values_from = total_deaths) %>%
  mutate(Sex = factor(Sex)) %>%
  arrange(Sex) %>%
  select(Sex, everything())

# Create table to present the data]
Fig_2 <- kable(summary_tab[, -1], 
               caption = "Deaths from unintentional injuries in Scotland",
               format = "latex", booktabs = TRUE) %>%
  kable_styling(font_size = 10) %>%
  pack_rows(tab_kable, colnum = 1,
            index = table(fct_inorder(summary_tab$Sex), useNA = "no"))

