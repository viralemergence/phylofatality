## phylofatality 
## 01_generate species-level CFR with reconciled mammal taxonomy
## danbeck@ou.edu, carolinecummings@ou.edu 
## last update 2/11/2025

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

## load virion
## https://github.com/viralemergence/virion/blob/main/Virion/Virion.csv.gz
setwd("~/Desktop/GitHub/virion/Virion")
vir=vroom("virion.csv.gz")
vir %<>% filter(HostClass == 'mammalia')

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
vdata=vir[(vir$HostNCBIResolved==T & vir$VirusNCBIResolved==T),]
table(cfr$Virus %in% vir$Virus)

## remove humans
vdata=vdata[!vdata$Host=="homo sapiens",]

## simplify vdata to cfr viruses
vdata=vdata[vdata$Virus%in%cfr$Virus,]

## remove missing hosts
vdata=vdata[!is.na(vdata$Host),]

#fix VIRION 
vdata$VirusGenus <- ifelse(vdata$Host=="bos taurus" & vdata$Virus=="pestivirus a" &
         vdata$VirusGenus=="flavivirus", "pestivirus", vdata$VirusGenus)

#check (should be zero)
vdata %>% filter(Host=="bos taurus", Virus=="pestivirus a", VirusGenus=="flavivirus") %>% print()

#save for 02_summary statistics script
setwd("~/Desktop/GitHub/phylofatality/csv files")
#write_csv(vdata, "01_vdata.csv")

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
n_distinct(vdata$pair) #2,915 unique host-virus associations

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

## how many unique host-virus associations are PCR or isolation? #1459
vdata$evidence=ifelse(vdata$DetectionMethod_PCR.Sequencing==1 | vdata$DetectionMethod_Isolation.Observation==1,1,0)
table(vdata$evidence) 
table(vdata$VirusFamily,vdata$evidence)

#how many unique host-virus associations for each detection type?
vdata$evidence=ifelse(vdata$DetectionMethod_Isolation.Observation==1,1,0)
table(vdata$evidence) #797 isolation
vdata$evidence=ifelse(vdata$DetectionMethod_PCR.Sequencing==1,1,0)
table(vdata$evidence) #1103 pcr
vdata$evidence=ifelse(vdata$DetectionMethod_Antibodies==1,1,0)
table(vdata$evidence) #1323 antibodies
vdata$evidence=ifelse(vdata$DetectionMethod_Not.specified==1,1,0)
table(vdata$evidence) #2544 none

#more specifically, which are lacking strong evidence?
vdata$evidence=ifelse(vdata$DetectionMethod_Isolation.Observation==1 | 
                        vdata$DetectionMethod_PCR.Sequencing==1 |
                        vdata$DetectionMethod_Antibodies==1,1,0)
table(vdata$evidence) 
#408 are completely unspecified
#2507 are detected by at least 1 detection method

vdata$evidence=ifelse(vdata$DetectionMethod_Isolation.Observation==1 | 
                        vdata$DetectionMethod_PCR.Sequencing==1 |
                        vdata$DetectionMethod_Not.specified==1,1,0)
table(vdata$evidence) 
#159 are only detected by antibodies
#159/2915 = 5.5%

#which are only viral isolation or only PCR or both?
vdata$evidence=ifelse(vdata$DetectionMethod_Not.specified==1 | 
                        vdata$DetectionMethod_PCR.Sequencing==1 |
                        vdata$DetectionMethod_Antibodies==1,1,0)
table(vdata$evidence) 
#30 are only viral isolation

vdata$evidence=ifelse(vdata$DetectionMethod_Isolation.Observation==1 | 
                        vdata$DetectionMethod_Not.specified==1 |
                        vdata$DetectionMethod_Antibodies==1,1,0)
table(vdata$evidence) 
#160 are only pcr

vdata$evidence=ifelse(vdata$DetectionMethod_Isolation.Observation==1 | 
                        vdata$DetectionMethod_PCR.Sequencing==1,1,0)
table(vdata$evidence)

# 1459 are detected by viral isolation and/or pcr
# 441 are viral isolation AND pcr
# 1459/2915 > 50%

## load in host taxonomy
setwd("~/Desktop/GitHub/phylofatality/phylo")
taxa=read.csv('taxonomy_mamPhy_5911species.csv',header=T)
taxa$tip=taxa$Species_Name

## fix tip
taxa$species=sapply(strsplit(taxa$tip,'_'),function(x) paste(x[1],x[2],sep=' '))

## species in data
vdata$species=capitalize(vdata$Host)

## match
miss=setdiff(vdata$species,taxa$species) #79

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
miss=setdiff(vdata$species,taxa$species) #1

## remove missing species
vdata=vdata[!vdata$species%in%miss,] 

## save data
vraw=vdata #2912 unique host-virus obs

# merge and look at db data
tmp=merge(cfr,vdata,by="Virus")
db <- tmp %>% select(species, Virus, db)
db=merge(db,taxa, by="species")
db <- db %>% select(species, Virus, db, ord,fam,gen) %>% arrange(desc(db))
db2=aggregate(db~fam,db,mean,na.rm=T) 

## look at covs
cov <- tmp %>% select(species, Virus, CFR, VirusFamily)
cov <- cov %>% filter(VirusFamily=="coronaviridae") %>% arrange(desc(CFR))

## look at flavis
fla<- tmp %>% select(species, Virus, CFR, VirusFamily)  %>% filter(VirusFamily=="flaviviridae") %>% arrange(desc(CFR))
fla$Host <- fla$species %>% str_to_lower()
fla$species=NULL

fla2 <- vir %>% filter(HostFamily=="vespertilionidae") %>% select(Host, HostFamily)
fla2 <- merge(fla, fla2, by="Host") %>% unique()
nrow(fla2) #25
sum(fla2$Virus=="japanese encephalitis virus") #10 
sum(fla2$Virus=="west nile virus") #6

fla3<- vir %>% filter(Virus=="japanese encephalitis virus" & HostFamily=="vespertilionidae") %>% unique() %>% select(Host, Virus, DetectionMethod, everything())
nrow(fla3) #52
sum(fla3$DetectionMethod=="Antibodies") #12
fla3$combo<- paste0(fla3$DetectionMethod, "_", fla3$Host)
fla4<- unique(fla3$combo) %>% as.data.frame()
nrow(fla4) #27, and 7 out of 27 are antibodies; 4/ 27 are PCR


fla <- vir %>% filter(Virus=="japanese encephalitis virus" & HostFamily=="vespertilionidae")
fla<- vir %>% filter(Virus=="japanese encephalitis virus" & HostOrder=="chiroptera") %>% select(Host, Virus, HostFamily) %>% unique()
fla<- vir %>% filter(Virus=="yellow fever virus" & HostOrder=="chiroptera") %>% select(Host, Virus, HostFamily) %>% unique()

fla$num <- 1

aggregate(num ~ HostFamily, data = fla, sum) %>% print()

## for each host species, fraction of all viruses that can infect humans
tmp$vir=1
tmp=aggregate(cbind(vir,htrans)~species,tmp,sum,na.rm=T) #number total viruses and number of viruses with human-human transmission per host
tmp$on.frac=tmp$htrans/tmp$vir ## percent of viruses a host has that exhibit human-human transmission 

## mean/max across all viruses per host
vdata %<>% left_join(cfr) %>%
  group_by(species) %>% 
  dplyr::summarize(meanCFR = mean(CFR),
                   maxCFR = max(CFR),
                   meanDB=mean(db),
                   virusesWithCFR = n())

## fix tmp names (get total viruses with OT/host)
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
                     meanDB=mean(db),
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
vdata=merge(vdata,vset,by="species",all=T) #983 obs (unique species)

## pubmed citations
library(easyPubMed)

## function
counter=function(name){
  as.numeric(as.character(get_pubmed_ids(gsub(' ','-',name))$Count))
}
citations=c()

## loop through
for(i in 1:length(vdata$species)) {
  citations[i]=counter(vdata$species[i])
  print(i)
}

## compile
cites=data.frame(species=vdata$species,
                 cites=citations)

## save 
setwd("~/Desktop/GitHub/phylofatality/csv files")
#write_csv(cites,"01_cites.csv")
cites <-read_csv("01_cites.csv")

## merge
vdata=merge(vdata,cites,by='species')

## clean
rm(cites,citations,i,counter)

## export
## save for 02_summary statistics script
setwd("~/Desktop/GitHub/phylofatality/csv files")
#write_csv(vdata,"01_CFRBySpecies.csv")
