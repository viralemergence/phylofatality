## phylofatality 
## 01.1_Compare species-level data between old and new VIRION
## danbeck@ou.edu, carolinecummings@ou.edu 
## last update 5/27/2025

## clean environment & plots
rm(list=ls()) 
graphics.off()
gc()

## load packages
library(tidyverse)
library(vroom)
library(magrittr)
library(ape)
library(plyr)
library(purrr)
library(Hmisc)
library(fastDummies)

## [1] load in VIRIONS
## load new virion
## download "full dataset" https://github.com/viralemergence/virion 
setwd("~/Desktop/GitHub/")
vir=vroom("virion.csv.gz")
vir %<>% filter(HostClass == 'mammalia')

## old virion
setwd("~/Desktop/GitHub/virion/Virion")
virOLD=vroom("virion.csv.gz")
virOLD %<>% filter(HostClass == 'mammalia')

## [2] load in CFR datq 
## load cfr data
setwd("~/Desktop/GitHub/phylofatality/data")
cfr1=read_csv("loose_data.csv.txt")
cfr2=read_csv("stringent_data.csv.txt")

## categorize
cfr1$cat="loose"
cfr2$cat="stringent"
cfr=bind_rows(cfr1, cfr2)
rm(cfr1,cfr2)

## save to use in 02_summary stats script
setwd("~/Desktop/GitHub/phylofatality/csv files")
#write_csv(cfr, "01_cfr.csv")

## rename based on CFR average
## human.trans is on scale 1-4 from 0 ---> endemic
cfr %<>% dplyr::select(SppName_ICTV_MSL2018b, CFR_avg, human.trans, death_burden_since_1950) %>%
  dplyr::rename(Virus = SppName_ICTV_MSL2018b, CFR = CFR_avg, onward=human.trans, db=death_burden_since_1950)

## fix with virion naming
cfr %<>% mutate(Virus = str_to_lower(Virus))

## fix NAs for DB
## rotavirus A is odd, low data availability
cfr$db[is.na(cfr$db)] <- 0
cfr$db<- ifelse(cfr$Virus=="rotavirus a", NA, cfr$db )
raw<- cfr

## [3] name matching for VIRIONS
## check name matching
miss <- setdiff(cfr$Virus,vir$Virus) %>% as.data.frame() #68

## left is old (CFR) and right is new (match virion)
rec <- c("flexal mammarenavirus"="mammarenavirus flexalense",
         "kasokero orthonairovirus"="orthonairovirus kasokeroense",
         "tacaribe mammarenavirus"="mammarenavirus tacaribeense",
         "rio bravo virus"="orthoflavivirus bravoense",
         "bagaza virus"="orthoflavivirus bagazaense",
         "cali mammarenavirus"="mammarenavirus caliense",
         "carnivore amdoparvovirus 1"="amdoparvovirus carnivoran1",
         "mobala mammarenavirus"="mammarenavirus praomyidis",
         "modoc virus"="orthoflavivirus modocense",
         "pestivirus a"="pestivirus bovis",
         "simian immunodeficiency virus"="lentivirus simimdef",
         "thailand orthohantavirus"="orthohantavirus thailandense",
         "tioman pararubulavirus"="pararubulavirus tiomanense",
         "dugbe orthonairovirus"="orthonairovirus dugbeense",
         "colorado tick fever virus" = "colorado tick fever coltivirus",
         "isfahan vesiculovirus"="vesiculovirus isfahan",
         "ilheus virus"="orthoflavivirus ilheusense",
         "japanese encephalitis virus"="orthoflavivirus japonicum",
         "murray valley encephalitis virus"="orthoflavivirus murrayense",
         "saint louis encephalitis virus"="orthoflavivirus louisense",
         "indiana vesiculovirus"="vesiculovirus indiana",
         "alagoas vesiculovirus"="vesiculovirus alagoas",
         "new jersey vesiculovirus"="vesiculovirus newjersey",
         "andes orthohantavirus"="orthohantavirus andesense",
         "argentinian mammarenavirus"="mammarenavirus juninense",
         "australian bat lyssavirus"="lyssavirus australis",
         "banzi virus"= "orthoflavivirus banziense",
         "bayou orthohantavirus"="orthohantavirus bayoui",
         "black creek canal orthohantavirus"="orthohantavirus nigrorivense",
         "brazilian mammarenavirus"="mammarenavirus brazilense",
         "california encephalitis orthobunyavirus"="orthobunyavirus encephalitidis",
         "caraparu orthobunyavirus"="orthobunyavirus caraparuense",
         "chapare mammarenavirus"="mammarenavirus chapareense",
         "colorado tick fever virus" = "colorado tick fever coltivirus",
         "dobrava-belgrade orthohantavirus"="orthohantavirus dobravaense",
         "duvenhage lyssavirus"="lyssavirus duvenhage",
         "european bat 1 lyssavirus"="lyssavirus hamburg",
         "european bat 2 lyssavirus"="lyssavirus helsinki",
         "guanarito mammarenavirus"="mammarenavirus guanaritoense",
         "hantaan orthohantavirus"="orthohantavirus hantanense",
         "hendra henipavirus"="henipavirus hendraense",
         "irkut lyssavirus"="lyssavirus irkut",
         "kokobera virus"="orthoflavivirus kokoberaorum",
         "laguna negra orthohantavirus"="orthohantavirus negraense",
         "louping ill virus"="orthoflavivirus loupingi",
         "lujo mammarenavirus"="mammarenavirus lujoense",
         "lymphocytic choriomeningitis mammarenavirus"="mammarenavirus choriomeningitidis",
         "machupo mammarenavirus"="mammarenavirus machupoense",
         "mammalian 1 orthobornavirus"="orthobornavirus bornaense",
         "menangle pararubulavirus"="pararubulavirus menangleense",
         "omsk hemorrhagic fever virus"="orthoflavivirus omskense",
         "orthohepevirus a"="paslahepevirus balayani",
         "powassan virus"="orthoflavivirus powassanense",
         "puumala orthohantavirus"="orthohantavirus puumalaense",
         "rabies lyssavirus"="lyssavirus rabies",
         "rift valley fever phlebovirus"="phlebovirus riftense",
         "sealpox virus"="grey sealpox virus",
         "severe acute respiratory syndrome-related coronavirus-2"="betacoronavirus pandemicum",
         "sin nombre orthohantavirus"="orthohantavirus sinnombreense",
         "sosuga pararubulavirus"="pararubulavirus sosugaense",
         "suid alphaherpesvirus 1"="varicellovirus suidalpha1",
         "tick-borne encephalitis virus"="orthoflavivirus encephalitidis",
         "tula orthohantavirus"="orthohantavirus tulaense",
         "whitewater arroyo mammarenavirus"="mammarenavirus whitewaterense",
         "ebolavirus"="orthoebolavirus zairense",
         "kyasanur forest disease virus"="orthoflavivirus kyasanurense",
         "shuni orthobunyavirus"="orthobunyavirus shuniense",
         "choclo orthohantavirus"="orthohantavirus chocloense"
)

cfr %<>% mutate(Virus = recode(Virus, !!!rec))

cfr$Virus[str_detect(cfr$Virus,'middle')] <- "middle east respiratory syndrome-related coronavirus"

## recheck
setdiff(cfr$Virus,vir$Virus) #0!

## as data frame
cfr=data.frame(cfr)

## unique onward
honward=cfr[!duplicated(cfr$Virus),]
honward$CFR=NULL

## code strict/loose
cmeans=aggregate(CFR~Virus,cfr,mean,na.rm=T)
cmeans$CFR=cmeans$CFR/100

## merge
cfr=merge(cmeans,honward,by="Virus")
rm(cmeans,honward)

## recode
cfr$type=ifelse(is.na(cfr$onward),"loose","stringent")

## human.trans--> make binary (causes human-human transmission y/n)
cfr$htrans=ifelse(cfr$onward==1,0,1)

## trim virion to NCBI resolved
vdata=vir[(vir$HostNCBIResolved==T & vir$VirusNCBIResolved==T),]
table(cfr$Virus %in% vir$Virus) # 119

## remove humans
vdata=vdata[!vdata$Host=="homo sapiens",]

## simplify vdata to cfr viruses
vdata=vdata[vdata$Virus%in%cfr$Virus,]

## remove missing hosts
vdata=vdata[!is.na(vdata$Host),]

#check (should be zero)
vdata %>% filter(Host=="bos taurus", Virus=="pestivirus bovis", VirusGenus=="flavivirus") %>% print() # good to go

#save for 02_summary statistics script
setwd("~/Desktop/GitHub/phylofatality/csv files")
#write_csv(vdata, "01_vdata.csv")

## redo for old VIRION
## check name matching
cfr<- raw
setdiff(cfr$Virus,virOLD$Virus)
rec <- c("colorado tick fever virus" = "colorado tick fever coltivirus",
         "ebolavirus" = "zaire ebolavirus",
         "sealpox virus" = "seal parapoxvirus",
         "severe acute respiratory syndrome-related coronavirus-2" = "severe acute respiratory syndrome-related coronavirus")
cfr %<>% mutate(Virus = recode(Virus, !!!rec))
cfr$Virus[str_detect(cfr$Virus,'middle')] <- "middle east respiratory syndrome-related coronavirus"

## recheck
setdiff(cfr$Virus,virOLD$Virus)
rm(rec)

## as data frame
cfr=data.frame(cfr)

## unique onward
honward=cfr[!duplicated(cfr$Virus),]
honward$CFR=NULL

## code strict/loose
cmeans=aggregate(CFR~Virus,cfr,mean,na.rm=T)
cmeans$CFR=cmeans$CFR/100

## merge
cfr=merge(cmeans,honward,by="Virus")
rm(cmeans,honward)

## recode
cfr$type=ifelse(is.na(cfr$onward),"loose","stringent")

## human.trans--> make binary (causes human-human transmission y/n)
cfr$htrans=ifelse(cfr$onward==1,0,1)

## trim virion to NCBI resolved
vdataOLD=virOLD[(virOLD$HostNCBIResolved==T & virOLD$VirusNCBIResolved==T),]
table(cfr$Virus %in% virOLD$Virus)

## remove humans
vdataOLD=vdataOLD[!vdataOLD$Host=="homo sapiens",]

## simplify vdata to cfr viruses
vdataOLD=vdataOLD[vdataOLD$Virus%in%cfr$Virus,]

## remove missing hosts
vdataOLD=vdataOLD[!is.na(vdataOLD$Host),]

#fix VIRION 
vdataOLD$VirusGenus <- ifelse(vdataOLD$Host=="bos taurus" & vdataOLD$Virus=="pestivirus a" &
                                vdataOLD$VirusGenus=="flavivirus", "pestivirus", vdataOLD$VirusGenus)

#check (should be zero)
vdataOLD %>% filter(Host=="bos taurus", Virus=="pestivirus a", VirusGenus=="flavivirus") %>% print()

## clean
rm(miss,raw)

## [4] compare the vdata and vdataOLD
vdata$pair<- paste0(vdata$Host,"_",vdata$Virus)
vdataOLD$pair<- paste0(vdataOLD$Host,"_",vdataOLD$Virus)

vdata$comp<- paste0(vdata$Host,"_",vdata$Virus,"_",vdata$Database)
vdataOLD$comp<- paste0(vdataOLD$Host,"_",vdataOLD$Virus,"_",vdataOLD$Database)

# compare
miss<- setdiff(vdata$pair,vdataOLD$pair) # >1400 added
miss<- setdiff(vdataOLD$pair,vdata$pair) # >2100 gone
miss<- setdiff(vdataOLD$comp,vdata$comp)

gone=vdataOLD[vdataOLD$comp%in%miss,]
gone <- unique(gone)
gone <- gone %>% select(comp, everything()) 

cov<- gone %>% filter(VirusFamily=="coronaviridae")

unique(cov$Database)
table(cov$Database) ## genbank and globi satand out
table(gone$Database) ## eid2, genbank, and globi

## merge to make comparisons later as needed
vdata$update<-"Y"
vdata$Release_Date=NULL
vdata$Collection_Date=NULL

vdataOLD$update<-"N"

comp_data<- rbind(vdata,vdataOLD)

## save
setwd("~/Desktop/GitHub/phylofatality/csv files")
#write.csv(comp_data, "01.1_compare_virions.csv")

## peruse things
tab<- table(gone$Virus)  %>% as.data.frame() 
tab <- tab %>% arrange(desc(Freq))
tab<-  table(gone$pair)  %>% as.data.frame() %>% arrange(desc(Freq))
tab<-  table(gone$comp)  %>% as.data.frame() %>% arrange(desc(Freq))
