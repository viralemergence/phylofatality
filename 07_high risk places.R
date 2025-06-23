## phylofatality
## 07_high risk places
## carolinecummings@ou.edu
## last update 6/20/2025

#troubleshooting:
#https://rfunctions.blogspot.com/2015/03/bivariate-maps-bivariatemap-function.html 

## clean environment & plots
rm(list=ls()) 
graphics.off()
gc()

#load packages
library(classInt)
library(dismo)
library(fasterize)
library(ggpubr)
library(maps)
library(raster)
#library(rgdal) ##
library(sp)
library(sf)
library(terra)
library(tidyverse)
library(XML) 

#1 Create a color matrix 
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
         xlab=xlab, ylab=ylab,cex.lab=1.3, xaxt="n", yaxt="n")
    par(mgp = c(1, 0.5, 0)) 
    for(i in 1:101){
      col.temp<-col.matrix[i-1,]
      points(my.data,rep((i-1)/100,101),pch=15,col=col.temp, cex=1)}
    seqs<-seq(0,100,(100/nquantiles))
    seqs[1]<-1
    col.matrix<-col.matrix[c(seqs), c(seqs)]}
  
  #color matrix
  col.matrix<-colmat(nquantiles=10, upperleft="goldenrod1", upperright="deeppink4",
                     bottomleft="#e8e8e8", bottomright="deepskyblue2", 
                     ylab = "anthropogenic footprint", xlab = "bat host richness")}

#2 Generate the bivariate map 
{
  bivariate.map<-function(rasterx, rastery, colormatrix=col.matrix, nquantiles=10){
    quanmean<- getValues(rasterx)
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
clade <- read.csv("05_pf_riskyspecies_20250609.csv")

#subset data 
all_mean<-clade %>% 
  filter(host=="mammal", var=="mean", virus=="all", factor=="1" | factor=="3")
all_db<- clade %>%
  filter(host=="mammal", var=="db", virus=="all", factor=="1")

fla_mean<- clade %>% 
  filter(host=="mammal", var=="mean", virus=="fla", factor=="3")
fla_db<- clade %>%
  filter(host=="mammal", var=="db", virus=="fla", factor=="2")
fla_ot<- clade %>%
  filter(host=="mammal", var=="ot", virus=="fla", factor=="2")

tog_mean <- clade %>% 
  filter(host=="mammal", var=="mean", virus=="tog", factor=="1")

#load in IUCN geographic range data
iucn <- st_read("~/Desktop/GitHub/footprint/IUCN/MAMMALS.shp")
iucn<- iucn %>% filter(order_=="CHIROPTERA")
iucn$binomial<- gsub(" ", "_", iucn$binomial)

#unique species in iucn dataset
names<- unique(iucn$binomial) %>% as.data.frame()
names$names<- names$.
names$.=NULL
names <- names %>% arrange(names)

## combine all bat names and compare to iucn to make sure we have the most covered
all_mean_species <- all_mean %>% dplyr::select(species) %>% unique()
all_db_species <- all_db %>% dplyr::select(species) %>% unique()
fla_db_species <- fla_db %>% dplyr::select(species) %>% unique()
fla_mean_species <- fla_mean %>% dplyr::select(species) %>% unique()
fla_ot_species <- fla_ot %>% dplyr::select(species) %>% unique()
tog_mean_species <- tog_mean %>% dplyr::select(species) %>% unique()

all.spec<- rbind(all_db_species,all_mean_species,
                 fla_db_species,fla_mean_species,fla_ot_species,
                 tog_mean_species)

all.spec <- all.spec %>% arrange(species) %>% unique() ##1268

#name reconciliation
miss=setdiff(all.spec$species, names$names) #105
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
miss=setdiff(all.spec$species, iucn$binomial) ##48

#create a blank raster and increase the resolution
#data downloaded from here: https://geodata.ucdavis.edu/climate/worldclim/1_4/grid/cur/
setwd("~/Desktop/GitHub/footprint/alt_2-5m_bil")
r <- raster("alt.bil")
r <- disaggregate((r)*0,2)

#subset geographic range data to risky species, etc. and create a RasterLayer
iucn_allmean<- iucn[iucn$binomial %in% all_mean$species,]
iucn_alldb<- iucn[iucn$binomial %in% all_db$species,] 
iucn_flamean <- iucn[iucn$binomial %in% fla_mean$species,] 
iucn_fladb <- iucn[iucn$binomial %in% fla_db$species,] 
iucn_flaot<- iucn[iucn$binomial %in% fla_ot$species,] 
iucn_togmean<- iucn[iucn$binomial %in% tog_mean$species,] 

map_allmean <- fasterize(iucn_allmean, r, fun="sum")
map_alldb<- fasterize(iucn_alldb, r, fun="sum")
map_flamean<- fasterize(iucn_flamean, r, fun="sum")
map_fladb<- fasterize(iucn_fladb, r, fun="sum")
map_flaot<- fasterize(iucn_flaot, r, fun="sum")
map_togmean<- fasterize(iucn_togmean, r, fun="sum")

#clean up environment
rm(names,miss,r,clade,iucn)

#5 Create Raster #2: Load in human footprint data and human population density data
#footprint data downloaded from here: https://sedac.ciesin.columbia.edu/data/set/wildareas-v3-2009-human-footprint 
setwd("~/Desktop/GitHub/footprint/footprint/")
footprint <-raster('~/Desktop/GitHub/footprint/footprint/wildareas-v3-2009-human-footprint.tif')

#Make sure the projections are the same (i.e., maps overlap correctly)
#Terra package is faster, so first convert both map and footprint RasterLayer --> SpatRaster using rast()
footprint <- rast(footprint)
map_allmean<- rast(map_allmean)
map_alldb<- rast(map_alldb)
map_flamean<- rast(map_flamean)
map_fladb<- rast(map_fladb)
map_flaot<- rast(map_flaot)
map_togmean<- rast(map_togmean)

footprint <- project(footprint, map_allmean)
terramapallmean<- map_allmean ## save for later

#then convert back to RasterLayer using raster()
footprint <- raster(footprint)
map_allmean<- raster(map_allmean)
map_alldb<- raster(map_alldb)
map_flamean<- raster(map_flamean)
map_fladb<- raster(map_fladb)
map_flaot<- raster(map_flaot)
map_togmean<- raster(map_togmean)

#get rid of outliers and ensure everything is aligned
footprint[footprint>50] <- NA  ## plot(footprint)

#check extent, resolution, and projections of rasters (need to be identical)
{
  print(extent(footprint))
  print(extent(map_allmean))
  
  print(res(footprint))
  print(res(map_allmean))
  
  print(crs(footprint))
  print(crs(map_allmean))
}

## test out layers
footprint <- aggregate(footprint, fact=10, fun=mean)
map_allmean <- aggregate(map_allmean, fact=10, fun=mean)
map_alldb <- aggregate(map_alldb, fact=10, fun=mean)
map_flamean <- aggregate(map_flamean, fact=10, fun=mean)
map_fladb <- aggregate(map_fladb, fact=10, fun=mean)
map_flaot <- aggregate(map_flaot, fact=10, fun=mean)
map_togmean <- aggregate(map_togmean, fact=10, fun=mean)

#6 Map the bivariate map
## bivar map
bivmap_allmean<-bivariate.map(map_allmean, footprint, colormatrix=col.matrix, nquantiles=10)
bivmap_alldb<-bivariate.map(map_alldb, footprint, colormatrix=col.matrix, nquantiles=10)
bivmap_flamean<-bivariate.map(map_flamean, footprint, colormatrix=col.matrix, nquantiles=10)
bivmap_fladb<-bivariate.map(map_fladb, footprint, colormatrix=col.matrix, nquantiles=10)
bivmap_flaot<-bivariate.map(map_flaot, footprint, colormatrix=col.matrix, nquantiles=10)
bivmap_togmean<-bivariate.map(map_togmean, footprint, colormatrix=col.matrix, nquantiles=10)

## tinker
# Get matrix blocks of interest
col.matrix
new.mat<- col.matrix
new.mat[1:8,]<- NA
new.mat[,1:8] <- NA

## add country borders
library(rnaturalearth)
library(rnaturalearthdata)

## load in country borders
countries <- ne_countries(scale = "medium", returnclass = "sf")

## make sure the projection matches our map's projection
crsallmean<- crs(terramapallmean)
countries <- st_transform(countries, crs = crsallmean)

# confirm they are the same
st_crs(countries) == st_crs(map_allmean)

# now plot
hrp_allmean<- terra:: plot(bivmap_allmean, frame.plot = TRUE, axes = F, box = T, 
                           add = F, legend = F, col = as.vector(new.mat), asp = 1)
plot(st_geometry(countries), add = TRUE, col = NA, border = "black", lwd = 0.5)

hrp_alldb<- terra:: plot(bivmap_alldb, frame.plot = TRUE, axes = F, box = T, 
                         add = F, legend = F, col = as.vector(new.mat), asp = 1)
plot(st_geometry(countries), add = TRUE, col = NA, border = "black", lwd = 0.5)

hrp_flamean<- terra:: plot(bivmap_flamean, frame.plot = TRUE, axes = F, box = T, 
                       add = F, legend = F, col = as.vector(new.mat), asp = 1)
plot(st_geometry(countries), add = TRUE, col = NA, border = "black", lwd = 0.5)

hrp_fladb<- terra:: plot(bivmap_fladb, frame.plot = TRUE, axes = F, box = T, 
                     add = F, legend = F, col = as.vector(new.mat), asp = 1)
plot(st_geometry(countries), add = TRUE, col = NA, border = "black", lwd = 0.5)

hrp_flaot<- terra:: plot(bivmap_flaot, frame.plot = TRUE, axes = F, box = T, 
                     add = F, legend = F, col = as.vector(new.mat), asp = 1)
plot(st_geometry(countries), add = TRUE, col = NA, border = "black", lwd = 0.5)

hrp_togmean<- terra:: plot(bivmap_togmean, frame.plot = TRUE, axes = F, box = T, 
                       add = F, legend = F, col = as.vector(new.mat), asp = 1)
plot(st_geometry(countries), add = TRUE, col = NA, border = "black", lwd = 0.5)


## add country borders over it
#plot(st_geometry(countries), add = TRUE, col = NA, border = "black", lwd = 0.5)
## save as PDF

## make a list of all the bivmaps you want to loop through
all.risk.maps <- list(
  allmean = bivmap_allmean,
  alldb = bivmap_alldb,
  flamean = bivmap_flamean,
  fladb = bivmap_fladb,
  flaot = bivmap_flaot,
  togmean = bivmap_togmean)

# save empty list where results will go
results <- list()

for (name in names(all.risk.maps)) {
  x <- all.risk.maps[[name]]
  
# Convert each raster cell to a point, and assign it a value (riskiness)
point <- as.data.frame(rasterToPoints(x))

# so for every coordinate (x,y) there is an assigned "risk" value
colnames(point) <- c("x", "y", "value")

# Filter only colored (non-NA) cells
#save the colors we care about in a vector called "color vector"
color_vector <- as.vector(new.mat)

#now in the full matrix, we give ever matrix hex color a value, and most of the 
## matrix hex colors are NA, since we got rid of them above
color_lookup <- data.frame(value = seq_along(color_vector), hex = color_vector)
color_lookup <- color_lookup[!is.na(color_lookup$hex), ]
print(color_lookup)

## now we know the only colors we care about correspond to these "values"
point <- point %>% filter(value==97 |value==98 |value==98 
                          |value==108 |value==109 |value==110
                          |value==119 |value==120 |value==121)

# Convert to spatial points
## now we take those points in the dataframe and transform them back to a spatial object
## and we specify the crs we use
point.sf <- st_as_sf(point, coords = c("x", "y"), crs = crs(countries))

# turn countries into spatial object, and make sure its CRS is the same as the points.sf CRS
countries.sf <- st_as_sf(countries)
countries.sf <- st_transform(countries.sf, crs = st_crs(point.sf))

## locate which country each point falls in
## kinda like the bivariate map situation above
joined <- st_join(point.sf, countries.sf, join = st_intersects)

# Now joined contains country names for each raster cell with color that we care about
places <- table(joined$name_long) %>% as.data.frame()

print(places)

## country counts organized
places <- places %>% arrange(desc(Freq))
names(places) <- c("country","freq")

## save
results[[name]] <- places


}

final <- bind_rows(results, .id = "source")


final$source <- factor(final$source, levels=c("allmean",
                                              "flamean",
                                              "togmean",
                                              "flaot",
                                              "alldb",
                                              "fladb" ))

final <- final %>% arrange(source, desc(freq))

## save
setwd('~/Desktop/GitHub/phylofatality/csv files')
write.csv(final, "07_high risk places.csv")

