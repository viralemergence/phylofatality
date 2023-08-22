## phylofatality 01_generate species-level CFR
## danbeck@ou.edu 

## clean environment & plots
rm(list=ls()) 
graphics.off()
gc()

## load packages
library(tidyverse)
library(vroom)
library(magrittr)

## load virion
setwd("/Users/danielbecker/Desktop/GitHub/virion/Virion")
vir=vroom("virion.csv.gz")
vir %<>% filter(HostClass == 'mammalia')

## load cfr
setwd("/Users/danielbecker/Desktop/GitHub/phylofatality")
cfr1=read_csv("loose_data.csv.txt")
cfr2=read_csv("stringent_data.csv.txt")

## categorize
cfr1$cat="loose"
cfr2$cat="stringent"
cfr=bind_rows(cfr1, cfr2)
rm(cfr1,cfr2)

## rename based on CFR average
cfr %<>% select(SppName_ICTV_MSL2018b, CFR_avg, human.trans) %>%
  rename(Virus = SppName_ICTV_MSL2018b, CFR = CFR_avg, onward=human.trans)

## group
cfr %<>% group_by(Virus) %>% 
  summarize(CFR = 0.01*mean(CFR))

## fix with virion naming
cfr %<>% mutate(Virus = str_to_lower(Virus))

## name matching
table(cfr$Virus %in% vir$Virus)
cfr$Virus[!(cfr$Virus %in% vir$Virus)]
rec <- c("colorado tick fever virus" = "colorado tick fever coltivirus",
         "ebolavirus" = "zaire ebolavirus",
         "sealpox virus" = "seal parapoxvirus",
         "severe acute respiratory syndrome-related coronavirus-2" = "severe acute respiratory syndrome-related coronavirus")
cfr %<>% mutate(Virus = recode(Virus, !!!rec))
cfr$Virus[str_detect(cfr$Virus,'middle')] <- "middle east respiratory syndrome-related coronavirus"

## recheck
table(cfr$Virus %in% vir$Virus)

## regroup
cfr %<>% group_by(Virus) %>% 
  summarize(CFR = mean(CFR))

## trim virion to NCBI resolved
vdata=vir[(vir$HostNCBIResolved==T & vir$VirusNCBIResolved==T),]
table(cfr$Virus %in% vir$Virus)

## simplify vdata to cfr viruses
vdata=vdata[vdata$Virus%in%cfr$Virus,]

## remove missing host
vdata=vdata[!is.na(vdata$Host),]

## filter virion
vdata %<>%
  select(Host, Virus, VirusGenus, VirusFamily) %>% 
  distinct() %>% drop_na()

## save alternative
vraw=vdata

## mean/max across all viruses per host for now
vdata %<>% left_join(cfr) %>%
  group_by(Host) %>% 
  summarize(meanCFR = mean(CFR),
            maxCFR = max(CFR),
            virusesWithCFR = n())

# ## group by host and virus family
# vraw %<>% left_join(cfr) %>%
#   group_by(Host,VirusFamily) %>% 
#   summarize(meanCFR = mean(CFR),
#             maxCFR = max(CFR),
#             virusesWithCFR = n())
# 
# ## from long to wide
# vwide1=spread(vraw,VirusFamily,meanCFR)

## export
setwd("/Users/danielbecker/Desktop/GitHub/phylofatality")
write_csv(vdata, "CFRBySpecies.csv")
