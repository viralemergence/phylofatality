## phylofatality 
## 02_summary statistics
## danbeck@ou.edu, carolinecummings2018@gmail.com
## last update: 2/11/2025

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
library(dplyr)

## load virion
setwd("~/Desktop/GitHub/virion/Virion")
vir=vroom("virion.csv.gz")
vir %<>% dplyr::filter(HostClass == 'mammalia')

#load in vdata
setwd("~/Desktop/GitHub/phylofatality/csv files")
vdata<- read_csv("01_vdata.csv")

## filter virion
vdata %<>%
  dplyr::select(Host, Virus, VirusGenus, VirusFamily, HostOrder) %>% 
  distinct() %>% drop_na()

## load in host taxonomy
setwd("~/Desktop/GitHub/phylofatality/phylo")
taxa=read.csv('taxonomy_mamPhy_5911species.csv',header=T)
taxa$tip=taxa$Species_Name

## fix tip
taxa$species=sapply(strsplit(taxa$tip,'_'),function(x) paste(x[1],x[2],sep=' '))

## species in data
vdata$species=capitalize(vdata$Host)

## match
miss=setdiff(vdata$species,taxa$species)

## flag
vdata$flag=ifelse(vdata$species%in%miss,1,0)

## fix data names from CLOVER
setwd("~/Desktop/GitHub/clover/clover/clover_0.1_mammalviruses/phylogenies")
tdata=read.csv("mammal_phylo_translations.csv")
tdata=tdata[!duplicated(tdata$Host),]

## merge
tdata$X=NULL
tdata$species=capitalize(tdata$Host)
vdata=merge(vdata,tdata[c("species","Host_Upham")],by="species",all.x=T)

## fix
vdata$species2=ifelse(vdata$flag==1,vdata$Host_Upham,vdata$species)

## clean
vdata$species=vdata$species2
vdata$species2=NULL
vdata$Host_Upham=NULL
vdata$flag=NULL

## reflag
vdata$flag=ifelse(is.na(vdata$species),1,0)
vdata$species=ifelse(vdata$flag==1,capitalize(vdata$Host),vdata$species)

## manual fix
vdata$species=revalue(vdata$species,
                     c("Allochrocebus preussi"="Cercopithecus preussi",
                       "Apodemus chejuensis"="Apodemus agrarius",
                       "Bos taurus x bison bison"="Bos taurus",
                       "Cavia cutleri"="Cavia tschudii",
                       "Cercopithecus doggetti"="Cercopithecus mitis",
                       "Cercopithecus kandti"="Cercopithecus mitis",
                       "Cercopithecus roloway"="Cercopithecus diana",
                       "Cricetomys ansorgei"="Cricetomys gambianus",
                       "Cricetulus griseus"="Cricetulus barabensis",
                       "Dobsonia magna"="Dobsonia moluccensis",
                       "Eothenomys eleusis"="Eothenomys melanogaster",
                       "Equus asinus x caballus"="Equus africanus",
                       "Equus caballus x asinus"="Equus caballus",
                       "Giraffa giraffa"="Giraffa camelopardalis",
                       "Hypsugo pulveratus"="Pipistrellus pulveratus",
                       "Laephotis capensis"="Neoromicia capensis",
                       "Loxodonta cyclotis"="Loxodonta africana",
                       "Macaca brunnescens"="Macaca ochreata",
                       "Macaca speciosa"="Macaca arctoides",
                       "Macronycteris gigas"="Hipposideros gigas",
                       "Microtus obscurus"="Microtus arvalis",
                       "Molossus ater"="Molossus rufus",
                       "Oligoryzomys utiaritensis"="Oligoryzomys nigripes",
                       "Oryzomys texensis"="Oryzomys palustris",
                       "Piliocolobus tholloni"="Procolobus badius",
                       "Rhabdomys dilectus"="Rhabdomys pumilio",
                       "Rhinolophus monoceros"="Rhinolophus pusillus",
                       "Zygodontomys cherriei"="Zygodontomys brevicauda"))

## rematch
miss=setdiff(vdata$species,taxa$species)

## remove missing species
vdata=vdata[!vdata$species%in%miss,]

#summmary stats
n_distinct(vdata$species) #983 unique
bats<- vdata %>% filter(HostOrder=="chiroptera") 
n_distinct(bats$species) #220 unique bats
n_distinct(vdata$Virus) #115
n_distinct(vdata$VirusFamily) #22

#how many mammals in each virus family
vdata%>% select(species, VirusFamily)%>% filter(VirusFamily=="coronaviridae") %>% n_distinct() #98
vdata%>% select(species, VirusFamily)%>% filter(VirusFamily=="flaviviridae") %>% n_distinct() #381
vdata%>% select(species, VirusFamily)%>% filter(VirusFamily=="rhabdoviridae") %>% n_distinct() #298
vdata%>% select(species, VirusFamily)%>%  filter(VirusFamily=="togaviridae") %>% n_distinct() #169
vdata%>% select(species, VirusFamily)%>% filter(VirusFamily=="paramyxoviridae") %>% n_distinct() #46

#how many bats in each virus family
bats%>% select(species, VirusFamily)%>% filter(VirusFamily=="coronaviridae") %>% n_distinct() #35
bats%>% select(species, VirusFamily)%>% filter(VirusFamily=="flaviviridae") %>% n_distinct() #75
bats%>% select(species, VirusFamily)%>% filter(VirusFamily=="rhabdoviridae") %>% n_distinct() #130
bats%>% select(species, VirusFamily)%>%  filter(VirusFamily=="togaviridae") %>% n_distinct() #36
bats%>% select(species, VirusFamily)%>% filter(VirusFamily=="paramyxoviridae") %>% n_distinct() #35

#host-virus associations, what happens when we cut out vector-borne?

#load in cfr data
setwd("~/Desktop/GitHub/phylofatality/csv files")
cfr<- read_csv("01_cfr.csv")

## fix with virion naming
cfr %<>% dplyr::rename(Virus = SppName_ICTV_MSL2018b, CFR = CFR_avg, onward=human.trans, db=death_burden_since_1950)
cfr %<>% mutate(Virus = str_to_lower(Virus))

## check name matching
setdiff(cfr$Virus,vir$Virus)
rec <- c("colorado tick fever virus" = "colorado tick fever coltivirus",
         "ebolavirus" = "zaire ebolavirus",
         "sealpox virus" = "seal parapoxvirus",
         "severe acute respiratory syndrome-related coronavirus-2" = "severe acute respiratory syndrome-related coronavirus")
cfr %<>% mutate(Virus = recode(Virus, !!!rec))
cfr$Virus[str_detect(cfr$Virus,'middle')] <- "middle east respiratory syndrome-related coronavirus"

## recheck
setdiff(cfr$Virus,vir$Virus)

#cut out VB
{
non_vb<- cfr %>% select(Virus, vFamily, IsVectorBorne) %>% filter(IsVectorBorne=="0")
miss= setdiff(vdata$Virus, non_vb$Virus)
vdata=vdata[!vdata$Virus%in%miss,]

#redo stats
#summmary stats
n_distinct(vdata$species) #784 unique
bats<- vdata %>% filter(HostOrder=="chiroptera") 
n_distinct(bats$species) #186 unique bats
n_distinct(vdata$Virus) #73
n_distinct(vdata$VirusFamily) #19

#how many mammals in each virus family
vdata%>% select(species, VirusFamily)%>% filter(VirusFamily=="coronaviridae") %>% n_distinct() #98
vdata%>% select(species, VirusFamily)%>% filter(VirusFamily=="flaviviridae") %>% n_distinct() #79
vdata%>% select(species, VirusFamily)%>% filter(VirusFamily=="rhabdoviridae") %>% n_distinct() #255
vdata%>% select(species, VirusFamily)%>%  filter(VirusFamily=="togaviridae") %>% n_distinct() #0
vdata%>% select(species, VirusFamily)%>% filter(VirusFamily=="paramyxoviridae") %>% n_distinct() #46

#how many bats in each virus family
bats%>% select(species, VirusFamily)%>% filter(VirusFamily=="coronaviridae") %>% n_distinct() #35
bats%>% select(species, VirusFamily)%>% filter(VirusFamily=="flaviviridae") %>% n_distinct() #0
bats%>% select(species, VirusFamily)%>% filter(VirusFamily=="rhabdoviridae") %>% n_distinct() #127
bats%>% select(species, VirusFamily)%>%  filter(VirusFamily=="togaviridae") %>% n_distinct() #0
bats%>% select(species, VirusFamily)%>% filter(VirusFamily=="paramyxoviridae") %>% n_distinct() #35
}

#load in complete vdata from 01_CFR Mean and replace vdata
setwd("~/Desktop/GitHub/phylofatality/csv files")
vdata<- read_csv("01_CFRBySpecies.csv")

#mean
mean(vdata$`meanCFR_all viruses`, na.rm=T) #0.2425732
mean(vdata$`maxCFR_all viruses`, na.rm=T) #0.3939759
mean(vdata$`on.frac_all viruses`, na.rm=T) #0.3529295
mean(vdata$`meanDB_all viruses`, na.rm=T) ## 124,114.3

## median death buden because of skew
hist(vdata$`meanDB_all viruses`)
median(vdata$`meanDB_all viruses`, na.rm=T) ## 1,633

#sd
sd1 <- sd(vdata$`meanCFR_all viruses`, na.rm=T)
sd2 <- sd(vdata$`maxCFR_all viruses`, na.rm=T)
sd3 <- sd(vdata$`on.frac_all viruses`, na.rm=T)
sd4 <- sd(vdata$`meanDB_all viruses`, na.rm=T)

#sample size
ss1<- length(vdata$`meanCFR_all viruses`)
ss2<- length(vdata$`maxCFR_all viruses`)
ss3<- length(vdata$`on.frac_all viruses`)
ss4<- length(vdata$`meanDB_all viruses`)

# se
se1 <- sd1/sqrt(ss1) #0.010
se2 <- sd2/sqrt(ss2) #0.013
se3 <- sd3/sqrt(ss3) #0.013
se4 <- sd4/sqrt(ss4) #12,928.56 ... look into doing CI

#summary stats bats
bdata=vdata[vdata$species%in%bats$species,]

#mean
mean(bdata$`meanCFR_all viruses`, na.rm=T) #0.600
mean(bdata$`maxCFR_all viruses`, na.rm=T) #0.806
mean(bdata$`on.frac_all viruses`, na.rm=T) #0.315
mean(bdata$`meanDB_all viruses`, na.rm=T) # 257,653.2

## median because of skew
median(bdata$`meanDB_all viruses`, na.rm=T) # 91,642.75

#sd
sd1 <- sd(bdata$`meanCFR_all viruses`, na.rm=T)
sd2 <- sd(bdata$`maxCFR_all viruses`, na.rm=T)
sd3 <- sd(bdata$`on.frac_all viruses`, na.rm=T)
sd4 <- sd(bdata$`meanDB_all viruses`, na.rm=T)

#sample size
ss1<- length(bdata$`meanCFR_all viruses`)
ss2<- length(bdata$`maxCFR_all viruses`)
ss3<- length(bdata$`on.frac_all viruses`)
ss4<- length(bdata$`meanDB_all viruses`)

# se
se1 <- sd1/sqrt(ss1) #0.027
se2 <- sd2/sqrt(ss2) #0.024
se3 <- sd3/sqrt(ss3) #0.028
se4 <- sd4/sqrt(ss4) # 42,037.29
