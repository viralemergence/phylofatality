
library(tidyverse); library(vroom); library(magrittr)
setwd("C:/Users/cjcar/Desktop/Phylofatality")

vir <- vroom("~/Github/virion/Virion/virion.csv.gz")
vir %<>% filter(HostClass == 'mammalia')

cfr1 <- read_csv("loose_data.csv.txt")
cfr2 <- read_csv("stringent_data.csv.txt")
cfr <- bind_rows(cfr1, cfr2)

cfr %<>% select(SppName_ICTV_MSL2018b, CFR_avg) %>%
  rename(Virus = SppName_ICTV_MSL2018b, CFR = CFR_avg)

cfr %<>% group_by(Virus) %>% 
  summarize(CFR = 0.01*mean(CFR))

cfr %<>% mutate(Virus = str_to_lower(Virus))

vir %<>%
  filter(Virus %in% cfr$Virus) %>%
  select(Host, Virus) %>% 
  distinct() %>% drop_na()

vir %<>% left_join(cfr) %>%
  group_by(Host) %>% 
  summarize(meanCFR = mean(CFR),
            maxCFR = max(CFR),
            virusesWithCFR = n())

write_csv(vir, "CFRBySpecies.csv")
