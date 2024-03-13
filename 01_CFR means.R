## phylofatality 
## 01_generate species-level CFR with reconciled mammal taxonomy
## danbeck@ou.edu 
## last update 3/13/24

## clean environxment & plots
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

## load virion
#setwd("~/Desktop/virion/Virion")
setwd("~/Desktop/GitHub/virion/Virion")
vir=vroom("virion.csv.gz")
vir %<>% filter(HostClass == 'mammalia')

## load cfr
#setwd("~/Desktop/phylofatality/data")
setwd("~/Desktop/GitHub/phylofatality/data")
cfr1=read_csv("loose_data.csv.txt")
cfr2=read_csv("stringent_data.csv.txt")

## categorize
cfr1$cat="loose"
cfr2$cat="stringent"
cfr=bind_rows(cfr1, cfr2)
rm(cfr1,cfr2)

## rename based on CFR average
cfr %<>% dplyr::select(SppName_ICTV_MSL2018b, CFR_avg, human.trans) %>%
  dplyr::rename(Virus = SppName_ICTV_MSL2018b, CFR = CFR_avg, onward=human.trans)

## fix with virion naming
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

## human.trans
cfr$htrans=ifelse(cfr$onward==1,0,1)

## trim virion to NCBI resolved
vdata=vir[(vir$HostNCBIResolved==T & vir$VirusNCBIResolved==T),]
table(cfr$Virus %in% vir$Virus)

## remove humans
vdata=vdata[!vdata$Host=="homo sapiens",]

## simplify vdata to cfr viruses
vdata=vdata[vdata$Virus%in%cfr$Virus,]

## remove missing host
vdata=vdata[!is.na(vdata$Host),]

## summarize detection method
table(vdata$DetectionMethod)

## make binary columns for detection
dums=dummy_cols(vdata["DetectionMethod"])

## merge
dums$DetectionMethod=NULL
vdata=data.frame(vdata,dums)
rm(dums)

## unique ID
vdata$pair=paste(vdata$Host,vdata$Virus)

## aggregate detection and filter
vdata=aggregate(cbind(DetectionMethod_Antibodies,
                DetectionMethod_Isolation.Observation,
                DetectionMethod_Not.specified,
                DetectionMethod_PCR.Sequencing)~pair+Host+Virus+VirusGenus+VirusFamily,data=vdata,sum)

## binary
vdata[c("DetectionMethod_Antibodies",
        "DetectionMethod_Isolation.Observation",
        "DetectionMethod_Not.specified",
        "DetectionMethod_PCR.Sequencing")]=ifelse(vdata[c("DetectionMethod_Antibodies",
               "DetectionMethod_Isolation.Observation",
               "DetectionMethod_Not.specified",
               "DetectionMethod_PCR.Sequencing")]>0,1,0)

## how many unique host-virus associations are PCR or isolation?
vdata$evidence=ifelse(vdata$DetectionMethod_PCR.Sequencing==1 | vdata$DetectionMethod_PCR.Sequencing==1,1,0)
table(vdata$evidence)
table(vdata$VirusFamily,vdata$evidence)

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
mis=setdiff(vdata$species,taxa$species)

## remove missing species
vdata=vdata[!vdata$species%in%miss,]

## save data
vraw=vdata

## for each host species, fraction of all viruses that can infect humans
tmp=merge(cfr,vdata,by="Virus")
tmp$vir=1
tmp=aggregate(cbind(vir,htrans)~species,tmp,sum,na.rm=T)
tmp$on.frac=tmp$htrans/tmp$vir

## mean/max across all viruses per host
vdata %<>% left_join(cfr) %>%
  group_by(species) %>% 
  dplyr::summarize(meanCFR = mean(CFR),
                   maxCFR = max(CFR),
                   virusesWithCFR = n())

## fix tmp names
tmp$virusesWithOT=tmp$vir
tmp$vir=NULL

## merge fraction zoonotic
vdata=merge(vdata,tmp,by="species",all.x=T)

## fix names
names(vdata)[2:ncol(vdata)]=paste(names(vdata)[2:ncol(vdata)],"_all viruses",sep="")

## tabulate unique number of hosts per viral family
tmp=merge(cfr,vraw,by="Virus")
vfam_hosts=sort(sapply(unique(tmp$VirusFamily),function(x){
  
  set=tmp[tmp$VirusFamily==x,]
  return(length(unique(set$Host)))
  
}))
vfam_hosts=data.frame(vfam_hosts)
vfam_hosts$VirusFamily=rownames(vfam_hosts)
rownames(vfam_hosts)=NULL
names(vfam_hosts)=c("hosts","VirusFamily")

## cutoff of n=30 or more host species for now
vfam=vfam_hosts[vfam_hosts$hosts>30,]

## function to derive vfam-specific responses
vfam_out=function(x){
  
  ## subset tmp by given virus family
  set=tmp[tmp$VirusFamily==x,]
  sraw=set
  
  ## calculate fraction of viruses that can infect humans
  set$vir=1
  
  ## ifelse for no htrans
  if(all(is.na(set$htrans))){
    
    ## as NA
    set=aggregate(cbind(vir)~species,set,sum,na.rm=T)
    set$htrans=NA
    set$on.frac=0
    set$virusesWithOT=NA
    
  }else{
    
  set=aggregate(cbind(vir,htrans)~species,set,sum,na.rm=T)
  set$on.frac=set$htrans/set$vir
  set$virusesWithOT=set$vir
  }
  
  ## fix tmp names
  set$vir=NULL
  
  ## mean/max across all viruses per host for this virus family
  sraw %<>% left_join(cfr) %>%
    group_by(species) %>% 
    dplyr::summarize(meanCFR = mean(CFR),
                     maxCFR = max(CFR),
                     virusesWithCFR = n())
  
  ## merge fraction zoonotic
  set=merge(sraw,set,by="species",all.x=T)
  
  ## rename
  names(set)[2:ncol(set)]=paste(names(set)[2:ncol(set)],"_",x,sep="")
  
  ## return
  return(set)
  
}

## run function, save as list
vlist=lapply(vfam_hosts$VirusFamily,vfam_out)

## merge all
vset=vlist %>% purrr::reduce(full_join,by="species")

## merge
vdata=merge(vdata,vset,by="species",all=T)

## export
#setwd("~/Desktop/phylofatality")
setwd("~/Desktop/GitHub/phylofatality")
write_csv(vdata,"CFRBySpecies.csv")