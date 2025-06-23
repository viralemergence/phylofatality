## phylofatality 
## 02_summary statistics
## danbeck@ou.edu, carolinecummings2018@gmail.com
## last update: 6/20/2025

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
setwd("~/Desktop/GitHub/")
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
miss=setdiff(vdata$species,taxa$species) # 103

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

## recheck
miss=setdiff(vdata$species,taxa$species) # 69

## fix genus modifications
vdata$species=gsub("Neoeptesicus|Cnephaeus","Eptesicus",vdata$species)

## remove sp.
vdata$drop=ifelse(grepl("sp.",vdata$species),1,0)
vdata=vdata[!vdata$drop==1,]
vdata$drop=NULL

## recheck
miss=setdiff(vdata$species,taxa$species) %>% as.data.frame()

## manual fix
{vdata$species=revalue(vdata$species,
                      c("Alexandromys fortis"="Microtus fortis",
                        "Alexandromys maximowiczii"="Microtus maximowiczii",
                        "Alexandromys oeconomus"="Microtus oeconomus",
                        "Apodemus chejuensis"="Apodemus agrarius",
                        "Apodemus ciscaucasicus"="Apodemus uralensis",
                        "Bubalus kerabau"="Bubalus arnee",
                        "Cephalophorus callipygus"="Cephalophus callipygus",
                        "Clethrionomys gapperi"="Myodes gapperi",
                        "Clethrionomys rutilus"="Myodes rutilus", 
                        "Cricetulus griseus"="Cricetulus barabensis",
                        "Dicotyles tajacu"="Pecari tajacu",
                        "Epomophorus pusillus"="Micropteropus pusillus",
                        "Giraffa giraffa"="Giraffa camelopardalis",
                        "Glossophaga mutica"="Glossophaga soricina",
                        "Heteromys salvini"="Liomys salvini",
                        "Kerivoula furva"="Kerivoula titania",
                        "Lophostoma silvicola"="Lophostoma silvicolum",
                        "Lyroderma lyra"="Megaderma lyra",
                        "Microtus obscurus"="Microtus arvalis",
                        "Microtus rossiaemeridionalis"="Microtus arvalis",
                        "Molossus nigricans"="Molossus rufus",
                        "Mops plicatus"="Chaerephon plicatus",
                        "Mops pumilus"="Chaerephon pumilus",
                        "Murina feae"="Murina aurata",
                        "Neogale frenata"="Mustela frenata",
                        "Neogale vison"="Neovison vison",
                        "Nothocricetulus migratorius"="Cricetulus migratorius",
                        "Oligoryzomys costaricensis"="Oligoryzomys fulvescens",
                        "Otaria byronia"="Otaria bryonia",
                        "Pekania pennanti"="Martes pennanti",
                        "Piliocolobus tholloni"="Procolobus badius",
                        "Pteropus medius"="Pteropus giganteus",
                        "Rhinolophus monoceros"="Rhinolophus pusillus",
                        "Stenocranius gregalis"="Microtus gregalis",
                        "Urva edwardsii"="Herpestes edwardsii"))}

## rematch
miss=setdiff(vdata$species,taxa$species) #0

## remove missing species
vdata=vdata[!vdata$species%in%miss,] 

#summmary stats
n_distinct(vdata$species) ## 889
bats<- vdata %>% filter(HostOrder=="chiroptera") 
n_distinct(bats$species) ## 202
n_distinct(vdata$Virus) ## 112
n_distinct(vdata$VirusFamily) # 23

#how many mammals in each virus family
vdata%>% select(species, VirusFamily)%>% filter(VirusFamily=="coronaviridae") %>% n_distinct() # 60
vdata%>% select(species, VirusFamily)%>% filter(VirusFamily=="flaviviridae") %>% n_distinct() # 344
vdata%>% select(species, VirusFamily)%>% filter(VirusFamily=="rhabdoviridae") %>% n_distinct() # 301
vdata%>% select(species, VirusFamily)%>%  filter(VirusFamily=="togaviridae") %>% n_distinct() # 150
vdata%>% select(species, VirusFamily)%>% filter(VirusFamily=="paramyxoviridae") %>% n_distinct() # 43
vdata%>% select(species, VirusFamily)%>% filter(VirusFamily=="poxviridae") %>% n_distinct() ## 122
vdata%>% select(species, VirusFamily)%>% filter(VirusFamily=="arenaviridae") %>% n_distinct() ## 100
tab<- table(vdata$VirusFamily) %>% as.data.frame() %>% arrange(desc(Freq))

#how many bats in each virus family
bats%>% dplyr::select(species, VirusFamily)%>% filter(VirusFamily=="coronaviridae") %>% n_distinct() ## 23
bats%>% select(species, VirusFamily)%>% filter(VirusFamily=="flaviviridae") %>% n_distinct() # 64
bats%>% select(species, VirusFamily)%>% filter(VirusFamily=="rhabdoviridae") %>% n_distinct() # 130
bats%>% select(species, VirusFamily)%>%  filter(VirusFamily=="togaviridae") %>% n_distinct() # 35
bats%>% select(species, VirusFamily)%>% filter(VirusFamily=="paramyxoviridae") %>% n_distinct() # 31
bats%>% select(species, VirusFamily)%>% filter(VirusFamily=="poxviridae") %>% n_distinct() ## 0
bats%>% select(species, VirusFamily)%>% filter(VirusFamily=="arenaviridae") %>% n_distinct() ## 10
tab<- table(bats$VirusFamily) %>% as.data.frame() %>% arrange(desc(Freq))

#host-virus associations, what happens when we cut out vector-borne?

#load in cfr data
setwd("~/Desktop/GitHub/phylofatality/csv files")
cfr<- read_csv("01_cfr.csv")

## fix with virion naming
cfr %<>% dplyr::rename(Virus = SppName_ICTV_MSL2018b, CFR = CFR_avg, onward=human.trans, db=death_burden_since_1950)
cfr %<>% mutate(Virus = str_to_lower(Virus))

## check name matching
setdiff(cfr$Virus,vir$Virus)

## check name matching
miss <- setdiff(cfr$Virus,vir$Virus) %>% as.data.frame() # 108
{
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
         "laguna negra orthohantavirus"="orthohantavirus mamorense",
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
         "choclo orthohantavirus"="orthohantavirus chocloense",
         
         "yaba monkey tumor virus"="yatapoxvirus yabapox",
         "mucambo virus"="alphavirus mucambo",
         "highlands j virus"="alphavirus highlandsj",
         "equine rhinitis a virus"="aphthovirus burrowsi",
         "everglades virus"="alphavirus everglades",
         "madariaga virus"="alphavirus madariaga",
         "palyam virus"="orbivirus palyamense",
         "pixuna virus"="alphavirus pixuna",
         "thottopalayam thottimvirus"="thottimvirus thottapalayamense",
         "una virus"="alphavirus una",
         "vesicular exanthema of swine virus"="vesivirus exanthema",
         "great island virus"= "orbivirus magninsulae",
         "avian orthoavulavirus 1"="orthoavulavirus javaense",
         "influenza a virus"="alphainfluenzavirus influenzae",
         "tonate virus"="alphavirus tonate",
         "west nile virus"="orthoflavivirus nilense",
         "bovine papular stomatitis virus"="parapoxvirus bovinestomatitis",
         "camelpox virus"="orthopoxvirus camelpox",
         "chikungunya virus"="alphavirus chikungunya",
         "cowpox virus"="orthopoxvirus cowpox",
         "dengue virus"="orthoflavivirus denguei",
         "orthohantavirus negraense"="orthohantavirus mamorense",
         "lassa mammarenavirus"="mammarenavirus lassaense",
         "macacine alphaherpesvirus 1"="simplexvirus macacinealpha1",
         "marburg marburgvirus"="orthomarburgvirus marburgense",
         "monkeypox virus"="orthopoxvirus monkeypox",
         "nipah henipavirus"="henipavirus nipahense",
         "orf virus"="parapoxvirus orf",
         "primate t-lymphotropic virus 1"="deltaretrovirus pritlym1",
         "primate t-lymphotropic virus 2"="deltaretrovirus pritlym2",
         "primate t-lymphotropic virus 3"="deltaretrovirus pritlym3",
         "pseudocowpox virus"="parapoxvirus pseudocowpox",
         "ross river virus"= "alphavirus rossriver",
         "rotavirus a"="rotavirus alphagastroenteritidis",
         "seoul orthohantavirus"="orthohantavirus seoulense",
         "severe acute respiratory syndrome-related coronavirus"="betacoronavirus pandemicum",
         "tanapox virus"="yatapoxvirus tanapox",
         "yellow fever virus"="orthoflavivirus flavi",
         "zika virus"="orthoflavivirus zikaense",
         "changuinola virus"="orbivirus changuinolaense",
         "nelson bay orthoreovirus"="orthoreovirus nelsonense")}

cfr %<>% mutate(Virus = recode(Virus, !!!rec))

cfr$Virus[str_detect(cfr$Virus,'middle')] <- "middle east respiratory syndrome-related coronavirus"

## recheck
miss<- setdiff(cfr$Virus,vir$Virus) ## middle east respiratory syndrome-related coronavirus is the only missing (NA in VIRION)
## remove missing species
cfr=cfr[!cfr$Virus%in%miss,] 

#cut out VB
{
non_vb<- cfr %>% select(Virus, vFamily, IsVectorBorne) %>% filter(IsVectorBorne=="0")
miss= setdiff(vdata$Virus, non_vb$Virus)
vdata=vdata[!vdata$Virus%in%miss,]

#redo stats
#summmary stats
n_distinct(vdata$species) #701
bats<- vdata %>% filter(HostOrder=="chiroptera") 
n_distinct(bats$species) #171
n_distinct(vdata$Virus) #71
n_distinct(vdata$VirusFamily) #20

#how many mammals in each virus family
vdata%>% select(species, VirusFamily)%>% filter(VirusFamily=="coronaviridae") %>% n_distinct() # 06
vdata%>% select(species, VirusFamily)%>% filter(VirusFamily=="flaviviridae") %>% n_distinct() #75
vdata%>% select(species, VirusFamily)%>% filter(VirusFamily=="rhabdoviridae") %>% n_distinct() #262
vdata%>% select(species, VirusFamily)%>%  filter(VirusFamily=="togaviridae") %>% n_distinct() #0
vdata%>% select(species, VirusFamily)%>% filter(VirusFamily=="paramyxoviridae") %>% n_distinct() #43

#how many bats in each virus family
bats%>% select(species, VirusFamily)%>% filter(VirusFamily=="coronaviridae") %>% n_distinct() #23
bats%>% select(species, VirusFamily)%>% filter(VirusFamily=="flaviviridae") %>% n_distinct() #0
bats%>% select(species, VirusFamily)%>% filter(VirusFamily=="rhabdoviridae") %>% n_distinct() #127
bats%>% select(species, VirusFamily)%>%  filter(VirusFamily=="togaviridae") %>% n_distinct() #0
bats%>% select(species, VirusFamily)%>% filter(VirusFamily=="paramyxoviridae") %>% n_distinct() #31
}

#load in complete vdata from 01_CFR Mean and replace vdata
setwd("~/Desktop/GitHub/phylofatality/csv files")
vdata<- read_csv("01_CFRBySpecies.csv")

#mean
mean(vdata$`meanCFR_all viruses`, na.rm=T) ## 0.26
mean(vdata$`maxCFR_all viruses`, na.rm=T)  ## 0.42
mean(vdata$`on.frac_all viruses`, na.rm=T) ## 0.31
mean(vdata$`meanDB_all viruses`, na.rm=T) ## 107,866.2 
min(vdata$`meanDB_all viruses`, na.rm=T)  ## 0
max(vdata$`meanDB_all viruses`, na.rm=T) ## 2,581,976

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
se1 <- sd1/sqrt(ss1) ## 0.01
se2 <- sd2/sqrt(ss2) ## 0.01
se3 <- sd3/sqrt(ss3)  ## 0.01
se4 <- sd4/sqrt(ss4)  ## 11993.78

#summary stats bats
bdata=vdata[vdata$species%in%bats$species,]

#mean
mean(bdata$`meanCFR_all viruses`, na.rm=T) ## 0.64
mean(bdata$`maxCFR_all viruses`, na.rm=T) ## 0.84
mean(bdata$`on.frac_all viruses`, na.rm=T) ## 0.27
mean(bdata$`meanDB_all viruses`, na.rm=T) # 310,732.6
min(bdata$`meanDB_all viruses`, na.rm=T) # 0
max(bdata$`meanDB_all viruses`, na.rm=T) ## 2,581,976

## median because of skew
median(bdata$`meanDB_all viruses`, na.rm=T) # 137,393.8

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
se1 <- sd1/sqrt(ss1) #0.03
se2 <- sd2/sqrt(ss2) #0.02
se3 <- sd3/sqrt(ss3) #0.03
se4 <- sd4/sqrt(ss4) # 50,298.55
