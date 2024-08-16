## phylofatality
## 04_bivariate map
## danbeck@ou.edu, carolinecummings@ou.edu, Colin Carlson
## last update 8/16/2024

#troubleshooting:
#https://cran.r-project.org/web/packages/terra/terra.pdf 
#https://rfunctions.blogspot.com/2015/03/bivariate-maps-bivariatemap-function.html 

## clean environment & plots
rm(list=ls()) 
graphics.off()

#load packages
library(classInt)
library(dismo)
library(fasterize)
library(maps)
library(raster)
library(rgdal)
library(sp)
library(sf)
library(terra)
library(tidyverse)
library(XML)


#1 Create a color matrix (nquantiles=10)
{
colmat<-function(nquantiles=10, upperleft=rgb(0,150,235, maxColorValue=255), 
                 upperright=rgb(130,0,80, maxColorValue=255), bottomleft="grey", 
                 bottomright=rgb(255,230,15, maxColorValue=255), xlab="x label", 
                 ylab="y label"){
  my.data<-seq(0,1,.01)
  my.class<-classIntervals(my.data,n=nquantiles,style="quantile")
  my.pal.1<-findColours(my.class,c(upperleft,bottomleft))
  my.pal.2<-findColours(my.class,c(upperright, bottomright))
  col.matrix<-matrix(nrow = 101, ncol = 101, NA)
  for(i in 1:101){
    my.col<-c(paste(my.pal.1[i]),paste(my.pal.2[i]))
    col.matrix[102-i,]<-findColours(my.class,my.col)}
  plot(c(1,1),pch=19,col=my.pal.1, cex=0.5,xlim=c(0,1),ylim=c(0,1),frame.plot=F, 
       xlab=xlab, ylab=ylab,cex.lab=1.3)
  for(i in 1:101){
    col.temp<-col.matrix[i-1,]
    points(my.data,rep((i-1)/100,101),pch=15,col=col.temp, cex=1)}
  seqs<-seq(0,100,(100/nquantiles))
  seqs[1]<-1
  col.matrix<-col.matrix[c(seqs), c(seqs)]}

col.matrix<-colmat(nquantiles=10, upperright="#be64ac", upperleft="#3b4994",
                   bottomleft="#e8e8e8", bottomright="#5ac8c8", 
                   ylab = "Human footprint/pop density", xlab = "Bat host richness")

#color matrix option 2
#col.matrix<-colmat(nquantiles=10, upperleft="goldenrod1", upperright="deeppink4",
#                   bottomleft="#e8e8e8", bottomright="deepskyblue2", 
#                   ylab = "Human footprint", xlab = "Coronavirus hosts")

#color matrix option 3
#col.matrix<-colmat(nquantiles=10, upperright="goldenrod1", upperleft="deeppink4",
#                   bottomleft="#e8e8e8", bottomright="deepskyblue2", 
#                   ylab = "Human footprint", xlab = "Coronavirus hosts")

}

#2 Generate the bivariate map (nquantiles=10)
{
bivariate.map<-function(rasterx, rastery, colormatrix=col.matrix, nquantiles=10){
  quanmean<-getValues(rasterx)
  temp<-data.frame(quanmean, quantile=rep(NA, length(quanmean)))
  brks<-with(temp, quantile(unique(temp),na.rm=TRUE, probs = c(seq(0,1,1/nquantiles))))
  r1<-within(temp, quantile <- cut(quanmean, breaks = brks, labels = 2:length(brks),include.lowest = TRUE))
  quantr<-data.frame(r1[,2]) 
  quanvar<-getValues(rastery)
  temp<-data.frame(quanvar, quantile=rep(NA, length(quanvar)))
  brks<-with(temp, quantile(unique(temp),na.rm=TRUE, probs = c(seq(0,1,1/nquantiles))))
  r2<-within(temp, quantile <- cut(quanvar, breaks = brks, labels = 2:length(brks),include.lowest = TRUE))
  quantr2<-data.frame(r2[,2])
  as.numeric.factor<-function(x) {as.numeric(levels(x))[x]}
  col.matrix2<-colormatrix
  cn<-unique(colormatrix)
  for(i in 1:length(col.matrix2)){
    ifelse(is.na(col.matrix2[i]),col.matrix2[i]<-1,col.matrix2[i]<-which(col.matrix2[i]==cn)[1])}
  cols<-numeric(length(quantr[,1]))
  for(i in 1:length(quantr[,1])){
    a<-as.numeric.factor(quantr[i,1])
    b<-as.numeric.factor(quantr2[i,1])
    cols[i]<-as.numeric(col.matrix2[b,a])}
  r<-rasterx
  r[1:length(r)]<-cols
  return(r)}
}

#3 Create raster #1: geographic ranges of bats (CoV risky clade)
#first, load in risky species name data
setwd("~/Desktop/GitHub/phylofatality/csv files")
clade <- read_csv("pf_riskyspecies.csv")

#subset data to All viruses-MeanCFR and CoV-MeanCFR
all<-clade %>% 
  filter(host=="mammal", var=="mean", virus=="all", factor=="1" | factor=="3")
cov<-clade %>% filter(host=="mammal", var=="mean", virus=="cov", factor=="1")

#load in IUCN geographic range data
iucn <- st_read("~/Desktop/GitHub/phylofatality/data/IUCN/MAMMALS.shp")
iucn<- iucn %>% filter(order_=="CHIROPTERA")
iucn$binomial<- gsub(" ", "_", iucn$binomial)

#unique species in iucn dataset
names<- unique(iucn$binomial) %>% as.data.frame()
names$names<- names$.
names$.=NULL
names$names<- names[order(names$names),]

#name reconciliation
miss=setdiff(all$species, names$names) #87 missing
miss<- as.data.frame(miss)

#revalue=c("old tip"= "new tip")
{
iucn$binomial= plyr::revalue(iucn$binomial,
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


#drop reasoning
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
}
miss=setdiff(all$species, iucn$binomial) #45 missing

#create a blank raster and increase the resolution
setwd("~/Desktop/GitHub/phylofatality/data/alt_2-5m_bil")
r <- raster("alt.bil")
r <- disaggregate((r)*0,2)

#subset geographic range data to CoV risky species and create a RasterLayer
iucn_all<- iucn[iucn$binomial %in% all$species,]
iucn_cov<- iucn[iucn$binomial %in% cov$species,] 

map_all <- fasterize(iucn_all, r, fun="sum")
map_cov <- fasterize(iucn_cov, r, fun="sum")

# This adds zeros for the continental area so they aren't invisible
fix <- function(x) {sum(x,r,na.rm=TRUE)+r} 
map_all<- fix(map_all)
map_cov<- fix(map_cov)

#clean up environment
rm(names, miss)

#4 Generate heatmap of CoV clade visualization
{
#library(rasterVis)
#library(RColorBrewer)

#mycolors <- colorRampPalette(rev(brewer.pal(10,"Spectral")))(21)
#mycolors[1] <- "#C0C0C0"

#rasterVis::levelplot(map,  
#                     col.regions = mycolors,
#                     #at = seq(0, 15, 1),
#                     alpha = 0.5, 
#                     scales=list(alternating=FALSE),
#                     par.strip.text=list(cex=0),
#                     xlab = NULL, ylab = NULL,
#                     maxpixels = 5e6)
}

#5 Create Raster #2: Load in human footprint data and hum population density data
setwd("~/Desktop/GitHub/phylofatality/data/footprint/")
footprint <-raster('~/Desktop/GitHub/phylofatality/data/footprint/wildareas-v3-2009-human-footprint.tif')
humpop <-raster('~/Desktop/GitHub/phylofatality/data/footprint/gpw_v4_population_density_rev11_2020_1_deg.tif')

#Make sure the projections are the same
#Terra package is faster, so first convert both map and footprint RasterLayer --> SpatRaster
footprint <- rast(footprint)
humpop <- rast(humpop)
map_all <- rast(map_all)
map_cov <- rast(map_cov)

footprint <- project(footprint, map_all)
humpop <- project(humpop, map_all)

#then convert back to RasterLayer
footprint <- raster(footprint)
humpop <- raster(humpop)
map_all <- raster(map_all)
map_cov <- raster(map_cov)

#get rid of outliers and ensure everything is aligned
footprint[footprint>50] <- NA
humpop[humpop<0] <- NA
#foot <- footprint + map*0
#foot1 <- footprint + map*NA

#check extent, resolution, and projections of rasters (need to be identical)
{
print(extent(footprint))
print(extent(humpop))
print(extent(map_all))
print(extent(map_cov))

print(res(footprint))
print(res(footprint))
print(res(map_all))
print(res(map_cov))

print(crs(footprint))
print(crs(map_all))
}
#more experiments
#footprint1 <- aggregate(footprint, fact=1000, fun=mean)
#map_all1 <- aggregate(map_all, fact=1000, fun=mean)
#footprint2<- rast(footprint1)
#map_all2<- rast(map_all1)

##test out layers
#footprint <- aggregate(footprint, fact=100, fun=mean)
#humpop <- aggregate(humpop, fact=100, fun=mean)
#map_all <- aggregate(map_all, fact=100, fun=mean)
#map_cov <- aggregate(map_cov, fact=100, fun=mean) #takes 23 minutes to run to here

#check raster layers individually
my.colors = colorRampPalette(c("#be64ac","lightblue", "yellow","orangered", "red"))
plot(map_all,frame.plot=F,axes=F,box=F,add=F,legend.width=1,legend.shrink=1,col=my.colors(255)) 

#6 Map the bivariate map
bivmap_all_foot<-bivariate.map(footprint, map_all, colormatrix=col.matrix, nquantiles=10)
bivmap_all_hum<-bivariate.map(humpop, map_all, colormatrix=col.matrix, nquantiles=10)
bivmap_cov_foot<-bivariate.map(footprint, map_cov, colormatrix=col.matrix, nquantiles=10)
bivmap_cov_hum<-bivariate.map(humpop, map_cov, colormatrix=col.matrix, nquantiles=10)

all_foot<- terra:: plot(bivmap_all_foot, frame.plot = TRUE, axes = F, box = T, 
                        add = F, legend = F, col = as.vector(col.matrix), asp = 1)
#map(interior = T, add = T)

all_hum<- terra:: plot(bivmap_all_hum, frame.plot = TRUE, axes = F, box = T, 
                       add = F, legend = F, col = as.vector(col.matrix), asp = 1)
#map(interior = T, add = T)

cov_foot<- terra:: plot(bivmap_cov_foot, frame.plot = TRUE, axes = F, box = T, 
                        add = F, legend = F, col = as.vector(col.matrix), asp = 1)
#map(interior = T, add = T)

cov_hum<- terra:: plot(bivmap_cov_hum, frame.plot = TRUE, axes = F, box = T, 
                       add = F, legend = F, col = as.vector(col.matrix), asp = 1)
#map(interior = T, add = T)
