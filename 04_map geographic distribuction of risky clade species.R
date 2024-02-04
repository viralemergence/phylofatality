## phylofatality
## 06_map geographic distribuction of risky clade species
## danbeck@ou.edu, carolinecummings@ou.edu
## last update 2/4/2023

## clean environment & plots
rm(list=ls()) 
graphics.off()

## packages
library(devtools)
library(dplyr)
library(ggplot2)
library(magrittr)
library(plyr)
library(sp)
library(tidyverse)

#load in species data from all risky clades
setwd("~/Desktop/PCM Class/phylofatality/clean/csv files")
data=read.csv("pf_riskyspecies.csv")

## load in mammal shapefile pruned to bats (spatial data)
setwd("~/Desktop/PCM Class/phylofatality/clean")
bats=readRDS("bat shp.rds")

## tip
bats$tip=gsub(" ","_",bats$binomial)
data$tip=data$species

#filter data to include species in risky bat clades for mean cfr
rawdata=data
all_mam_mean= data %>% filter(virus=="all", host=="mammal", var=="mean", factor=="1" | factor=="3")
cov_mam_mean= data %>% filter(virus=="cov", host=="mammal", var=="mean", factor=="1")
fla_mam_mean= data %>% filter(virus=="fla", host=="mammal", var=="mean", factor=="2" | factor=="4")

data<- bind_rows(all_mam_mean, cov_mam_mean, fla_mam_mean )

## check missing
miss=setdiff(data$tip,bats$tip) #105 missing

## fix bats
#sanity check name matching
bats_names<- as.data.frame(bats$tip)

#revalue=c("old tip"= "new tip")
bats$tip=revalue(bats$tip,
                 c("Austronomus_australis"="Tadarida_australis",
                   "Austronomus_kuboriensis"="Tadarida_kuboriensis",
                   "Baeodon_alleni"="Rhogeessa_alleni",
                   "Baeodon_gracilis"="Rhogeessa_gracilis",
                   "Boneia_bidens"="Rousettus_bidens",
                   "Chaerephon_jobimena"="Tadarida_jobimena",
                   #"Coelops_robinsoni"="Coelops_hirsutus",
                   "Dermanura_azteca"="Dermanura_aztecus",
                   "Dermanura_cinerea"="Dermanura_cinereus",
                   "Dermanura_glauca"="Dermanura_glaucus",
                   "Dermanura_gnoma"="Dermanura_gnomus",
                   "Dermanura_rosenbergi"="Dermanura_rosenbergii",
                   "Dermanura_tolteca"="Dermanura_toltecus",
                   #"Dermanura_watsoni"="Dermanura_incomitatus",
                   "Diclidurus_isabella"="Diclidurus_isabellus",
                   "Gardnerycteris_crenulatum"="Mimon_crenulatum",
                   "Gardnerycteris_koepckeae"="Mimon_koepckeae",
                   "Hipposideros_pendleburyi"="Hipposideros_pendelburyi",
                   #"Hipposideros_pomona"="Paracoelops_megalotis",
                   "Hypsugo_affinis"="Falsistrellus_affinis",
                   "Hypsugo_alaschanicus"="Pipistrellus_alaschanicus",
                   "Hypsugo_anthonyi"="Pipistrellus_anthonyi",
                   "Hypsugo_arabicus"="Pipistrellus_arabicus",
                   "Hypsugo_ariel"="Pipistrellus_ariel",
                   "Hypsugo_cadornae"="Pipistrellus_cadornae",
                   "Hypsugo_eisentrauti"="Pipistrellus_eisentrauti",
                   "Hypsugo_joffrei"="Pipistrellus_joffrei",
                   "Hypsugo_kitcheneri"="Pipistrellus_kitcheneri",
                   "Hypsugo_lophurus"="Pipistrellus_lophurus",
                   "Hypsugo_macrotis"="Pipistrellus_macrotis",
                   "Hypsugo_musciculus"="Pipistrellus_musciculus",
                   "Hypsugo_pulveratus"="Pipistrellus_pulveratus",
                   "Hypsugo_savii"="Pipistrellus_savii",
                   "Hypsugo_vordermanni"="Pipistrellus_vordermanni",
                   "Lissonycteris_angolensis"="Myonycteris_angolensis",
                   "Lonchophylla_cadenai"="Hsunycteris_cadenai",
                   "Lonchophylla_pattoni"="Hsunycteris_pattoni",
                   #"Lonchophylla_thomasi"="Hsunycteris_thomasi",
                   "Lophostoma_occidentalis"="Lophostoma_aequatorialis",
                   #"Lophostoma_carrikeri"="Lophostoma_yasuni",
                   "Lyroderma_lyra"="Megaderma_lyra",
                   "Macronycteris_commersoni"="Hipposideros_commersoni",
                   "Macronycteris_gigas"="Hipposideros_gigas",
                   "Macronycteris_thomensis"="Hipposideros_thomensis",
                   "Macronycteris_vittatus"="Hipposideros_vittatus",
                   "Micronomus_norfolkensis"="Mormopterus_norfolkensis",
                   #"Miniopterus_schreibersii"="Miniopterus_fuliginosus",
                   #"Molossus_coibensis"="Molossus_barnesi"
                   "Mormopterus_kalinowskii"="Nyctinomops_kalinowskii",
                   "Murina_feae"="Murina_cineracea",
                   "Murina_lorelieae"="Murina_loreliae", 
                   #"Murina_harrisoni"="Murina_tiensa",
                   #"Murina_peninsularis"="Murina_cyclotis",
                   #"Myotis_albescens"="Nycticeius_aenobarbus",
                   #"Myotis_formosus"="Myotis_flavus",
                   #"Myotis_riparius"="Myotis_handleyi",
                   #"Myotis_simus"="Myotis_midastactus",
                   #"Natalus_mexicanus"="Natalus_lanatus",
                   #"Natalus_stramineus"="Natalus_saturatus",
                   "Nyctophilus_major"="Nyctophilus_timoriensis",
                   "Ozimops_loriae"="Mormopterus_loriae",
                   "Ozimops_planiceps"="Mormopterus_planiceps",
                   "Paremballonura_atrata"="Emballonura_atrata",
                   "Paremballonura_tiavato"="Emballonura_tiavato",
                   "Perimyotis_subflavus"="Pipistrellus_subflavus",
                   "Pipistrellus_anchietae"="Hypsugo_anchietae",
                   #"Pipistrellus_kuhlii"="Pipistrellus_deserti",
                   #"Pteropus_chrysoproctus"="Pteropus_argentatus",
                   #"Pteropus_pelewensis"="Pteropus_yapensis",
                   #"Rhinolophus_borneensis"="Rhinolophus_chaseni",
                   "Rhyneptesicus_nasutus"="Eptesicus_nasutus",
                   #"Scotonycteris_ophiodon"="Casinycteris_campomaanensis",
                   "Scotonycteris_ophiodon"="Casinycteris_ophiodon",
                   "Thainycteris_aureocollaris"="Arielulus_aureocollaris",
                   "Triaenops_rufus"="Triaenops_menamena",
                   "Vampyriscus_bidens"="Vampyressa_bidens",
                   "Vampyriscus_brocki"="Vampyressa_brocki",
                   "Vampyriscus_nymphaea"="Vampyressa_nymphaea"))
## check missing
miss=setdiff(data$tip,bats$tip) #48 missing
miss

#haven't found:
#Anoura_canishina

#drop
#Asellia_arabica - split from Asellia tridens
#Casinycteris_campomaanensis - no geo range info IUCN, relatively new species discovered (2014)
#Coelops_hirsutus - renamed Coelops_robinsoni
#Dermanura_incomitatus - molecularly the same as watsoni, and watsoni is on IUCN
#Desmodus_draculae - extinct
#Eumops_chiribaya - species is a relatively new discovery (2014)
#Glischropus_aquilus - species is a relatively new discovery (2015)
#Harpiocephalus_mordax - little info, also considered synonym of H harpia
#Hsunycteris_thomasi- basonym of Lonchophylla_thomasi
#Lonchophylla_inexpectata" - species is a relatively new discovery (2015)
#Lophostoma_yasuni-synonym of L carrikeri
#Miniopterus_mossambicus- species is a relatively new discovery (2013)
#Miniopterus_oceanensis - was a subspecies of M. schreibersii
#Miniopterus_fuliginosus - was a subspecies of M. schreibersii
#Molossus_barnesi - synonym of M. coibensis
#Murina_chrysochaetes - relatively new discovery, data deficient
#Murina guilleni- relatively new discovery, data deficient
#Murina_jaintiana - relatively new discovery, data deficient
#Murina pluvialis - relatively new discovery, data deficient
#Murina_tiensa - synonym of M. harrisoni
#Myotis formosus flavus - subspecies M. formosus
#Myotis_handleyi- relatively new, similar to M. riparius
#Myotis_midastactus- relatively new, similar to M. simus
#Myotis_phanluongi - relatively new discovery (2008)
#Natalus_lanatus - synonym of N. mexicanus
#Natalus stramineus saturatus - subspecies of N. stramineus
#Nycticeius_aenobarbus - distinct from Myotis albescens but often grouped with it
#Paracoelops_megalotis - now known as "Hipposideros_pomona" 
#Phoniscus_aerosa - data deficient
#Pipistrellus kuhlii deserti- subspecies
#Pipistrellus_sturdeei- extinct
#Pipistrellus_tenuis - missing in bat dataset, data deficient on IUCN
#Platyrrhinus_guianensis- relatively new discovery (2014)
#Pteropus_argentatus - synonym of Pteropus chrysoproctus
#Pteropus_brunneus - extinct
#Pteropus pilosus- extinct
#Pteropus_yapensis- synonym of Pteropus_pelewensis
#Rhinolophus borneensis chaseni - subspecies
#Rhinolophus francisi- relatively new discovery (2015)
#Rhinolophus_luctoides- preciously classified w R. luctus
#Rhinolophus thailandensis- relatively new (2009)
#Triaenops_parvus- now distinct, was a synonym of Triaenops_persicus
#Triaenops_rufus - old name, now Triaenops_menamena
#Uroderma_bakeri- missing in the bats dataset
#Vampyressa_elisabethae - missing in the bats dataset
#Vampyressa_sinchi - missing in the bats dataset


## drop missing
data=data[!data$tip%in%miss,]

## trim
bats=bats[bats$tip%in%data$tip,]

## save id
bats$id=rownames(bats@data)

#rgeos for R 4.1.2
#packageurl <- "https://cran.r-project.org/src/contrib/Archive/rgeos/rgeos_0.3-1.tar.gz"
#install.packages(packageurl, repos=NULL, type="source")

#simplify (with help from ChatGPT)
library(rgeos)
tol = 0.2

prob_geos = c()  #store geometries with problems in a vector

##loop through and simplify
lset=list()
for (i in 1:length(unique(bats$tip))) {
  
  # subset run
  set = bats[bats$tip == unique(bats$tip)[i],]
  
  tryCatch({ #run the code and catch errors that might occur
    # save data from the subset
    sdata = set@data
    
    # simplify
    shp = gSimplify(set, tol, topologyPreserve = TRUE)
    
    # fortify
    shp = data.frame(fortify(shp, region = "ID"))
    
    # merge with sdata by id
    sdata = sdata[c("id", "tip", "binomial")]
    shp = merge(shp, sdata, by = "id", all.x = TRUE)
    
    # save
    lset[[i]] = shp
    
    # print
    print(paste(i, "in", length(unique(bats$tip))))
  }, error = function(e) {
    # capture the index of the problematic geometry
    prob_geos = c(prob_geos, i)
    cat("Error in geometry", i, ": ", conditionMessage(e), "\n")
  })
}
# Print indices of problematic geometries
cat("Indices of problematic geometries:", prob_geos, "\n")

## convert to data
bset=do.call(rbind.data.frame,lset)

## clean
rm(bats)

## merge with data
bats=merge(bset,data,by="tip",all.x=T)

#convert factor variable to be a factor
bats$factor <- factor(bats$factor)

#save 
setwd("~/Desktop/PCM Class/phylofatality/clean/csv files")
#write.csv(bats, "bats_georanges.csv")

#load in bat geo range data
setwd("~/Desktop/PCM Class/phylofatality/clean/csv files")
bats=read.csv("bats_georanges.csv")


#save data into separate virus variables
fla<- bats%>% filter(virus=="fla")
all<- bats %>%filter(virus=="all")
cov<- bats %>%filter(virus=="cov")


## get world map
library(ggalt)
require(proj4)
library(ggthemes)
library(viridis)
library(mapproj)
wdata=map_data("world")
#wdata=wdata[-which(wdata$region=='Antarctica'),]

bmaps<-ggplot() +
  ## base layer
  geom_polygon(data=wdata, aes(x=long, y=lat, group=group),
               fill="grey90", colour="grey90", linewidth=0.2) +
  
  ## add shapefiles
  geom_polygon(data=all, aes(x=long, y=lat, group=paste(tip, group),
               fill=factor), alpha=0.25) +
  scale_fill_manual(values=c("purple2", "orange2","turquoise2", "magenta2"))+
  
  guides(fill = FALSE) +
  theme_void() +
  coord_map("gilbert", xlim = c(-180, 180))+
  ggtitle("Geographic range of risky bat hosts MeanCFR all viruses")+
  theme(plot.title = element_text(hjust = 0.5))

#save
print(bmaps)
setwd("~/Desktop/PCM Class/phylofatality/clean/figs")
#ggsave("map_allviruses.jpg", bmaps, device = "jpeg", width = 6, height = 6, units = "in")



#coronaviridae
bmaps<-ggplot() +
  ## base layer
  geom_polygon(data=wdata, aes(x=long, y=lat, group=group),
               fill="grey90", colour="grey90", linewidth=0.2) +
  
  ## add shapefiles
  ##for data=, "all" to see clades distributed across viruses overall, "cov" 
  ##to see cov clades, and "fla" to see fla clades
  geom_polygon(data=cov, aes(x=long, y=lat, group=paste(tip, group),
                             fill=factor), alpha=0.25) +
  scale_fill_manual(values=c("green4", "orange2","turquoise2", "magenta2"))+
  
  guides(fill = FALSE) +
  theme_void() +
  coord_map("gilbert", xlim = c(-180, 180))+
  ggtitle(expression("Geographic range of risky bat hosts MeanCFR-"~italic("Coronaviridae"))) +
  theme(plot.title = element_text(hjust = 0.5))
  
#save
print(bmaps)
setwd("~/Desktop/PCM Class/phylofatality/clean/figs")
#ggsave("map_cov.jpg", bmaps, device = "jpeg", width = 6, height = 6, units = "in")



#flaviviridae
bmaps<-ggplot() +
  ## base layer
  geom_polygon(data=wdata, aes(x=long, y=lat, group=group),
               fill="grey90", colour="grey90", linewidth=0.2) +
  
  ## add shapefiles
  geom_polygon(data=fla, aes(x=long, y=lat, group=paste(tip, group),
                             fill=factor), alpha=0.25) +
  scale_fill_manual(values=c("turquoise2", "magenta2","purple2", "orange2"))+
  
  guides(fill = FALSE) +
  theme_void() +
  coord_map("gilbert", xlim = c(-180, 180))+
  ggtitle(expression("Geographic range of risky bat hosts MeanCFR-"~italic("Flaviviridae"))) +
  theme(plot.title = element_text(hjust = 0.5))

#save
print(bmaps)
setwd("~/Desktop/PCM Class/phylofatality/clean/figs")
#ggsave("map_flav.jpg", bmaps, device = "jpeg", width = 6, height = 6, units = "in")

