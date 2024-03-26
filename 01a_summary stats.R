## phylofatality 
## 01a_summary statistics
## danbeck@ou.edu, carolinecummings2018@gmail.com
## last update: 3/26/2024

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
#setwd("~/Desktop/virion/Virion")
setwd("~/Desktop/GitHub/virion/Virion")
vir=vroom("virion.csv.gz")
vir %<>% dplyr::filter(HostClass == 'mammalia')

#load in vdata
setwd("~/Desktop/GitHub/phylofatality/data")
vdata<- read_csv("vdata.csv")

## filter virion
vdata %<>%
  dplyr::select(Host, Virus, VirusGenus, VirusFamily, HostOrder) %>% 
  distinct() %>% drop_na()

## load in host taxonomy
#setwd("~/Desktop/phylofatality/phylo")
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
#setwd("~/Desktop/clover/clover/clover_0.1_mammalviruses/phylogenies")
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
vdata%>% filter(VirusFamily=="coronaviridae") %>% n_distinct() #101
vdata%>% filter(VirusFamily=="flaviviridae") %>% n_distinct() #656
vdata%>% filter(VirusFamily=="rhabdoviridae") %>% n_distinct() #394
vdata%>% filter(VirusFamily=="togaviridae") %>% n_distinct() #251
vdata%>% filter(VirusFamily=="paramyxoviridae") %>% n_distinct() #67

#load in complete vdata from 01_CFR Mean and replace vdata
setwd("~/Desktop/GitHub/phylofatality/csv files")
vdata<- read_csv("CFRBySpecies.csv")

#mean
mean(vdata$`meanCFR_all viruses`, na.rm=T) #0.2425706
mean(vdata$`maxCFR_all viruses`, na.rm=T) #0.3939759
mean(vdata$`on.frac_all viruses`, na.rm=T) #0.3529295

#sd
sd1 <- sd(vdata$`meanCFR_all viruses`, na.rm=T)
sd2 <- sd(vdata$`maxCFR_all viruses`, na.rm=T)
sd3 <- sd(vdata$`on.frac_all viruses`, na.rm=T)

#sample size
ss1<- length(vdata$`meanCFR_all viruses`)
ss2<- length(vdata$`maxCFR_all viruses`)
ss3<- length(vdata$`on.frac_all viruses`)

# se
se1 <- sd1/sqrt(ss1) #0.0103
se2 <- sd2/sqrt(ss2) #0.013
se3 <- sd3/sqrt(ss3) #0.013

#summary stats bats
bdata=vdata[vdata$species%in%bats$species,]

#mean
mean(bdata$`meanCFR_all viruses`, na.rm=T) #0.521641
mean(bdata$`maxCFR_all viruses`, na.rm=T) #0.6981767
mean(bdata$`on.frac_all viruses`, na.rm=T) #0.2947801

#sd
sd1 <- sd(bdata$`meanCFR_all viruses`, na.rm=T)
sd2 <- sd(bdata$`maxCFR_all viruses`, na.rm=T)
sd3 <- sd(bdata$`on.frac_all viruses`, na.rm=T)

#sample size
ss1<- length(bdata$`meanCFR_all viruses`)
ss2<- length(bdata$`maxCFR_all viruses`)
ss3<- length(bdata$`on.frac_all viruses`)

# se
se1 <- sd1/sqrt(ss1) #0.03
se2 <- sd2/sqrt(ss2) #0.03
se3 <- sd3/sqrt(ss3) #0.03
