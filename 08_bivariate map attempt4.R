## phylofatality
## 08_bivariate map
## danbeck@ou.edu, carolinecummings@ou.edu, Colin Carlson
## last update 7/29/2024

#troubleshooting:
#https://cran.r-project.org/web/packages/terra/terra.pdf 
#https://rfunctions.blogspot.com/2015/03/bivariate-maps-bivariatemap-function.html 

## clean environment & plots
rm(list=ls()) 
graphics.off()

#load packages
library(classInt)
library(fasterize)
library(maps)
library(tidyverse)
library(raster)
library(rgdal)
library(dismo)
library(XML)
library(sp)
library(sf)
library(terra)

#1 Create a color matrix (nquantiles=10)
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

col.matrix<-colmat(nquantiles=10, upperleft="#be64ac", upperright="#3b4994",
                   bottomleft="#e8e8e8", bottomright="#5ac8c8", 
                   ylab = "Human footprint", xlab = "Coronavirus hosts")


#2 Generate the bivariate map (nquantiles=10)
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

#3 Create raster #1: geographic ranges of bats (CoV risky clade)
#first, load in risky species data
setwd("~/Desktop/GitHub/phylofatality/csv files")
clade <- read_csv("pf_riskyspecies.csv")

#subset data (94 species in CoV risky bat clade)
clade<-clade %>% 
  filter(host=="bat", var=="mean", virus=="cov")
clade$species<- gsub("_", " ", clade$species)

#load in IUCN geographic range data
#reduce resolution of shapefile
iucn <- st_read("~/Desktop/GitHub/phylofatality/data/IUCN")

#create a blank raster
setwd("~/Desktop/GitHub/phylofatality/data/alt_2-5m_bil")
r <- raster("alt.bil")
#r <- disaggregate((r)*0,2) ##Colin did this, but it makes it run a long time
r<- aggregate((r)*0,100)

#subset geographic range data to CoV risky species and create a RasterLayer
iucn<- iucn[iucn$binomial %in% clade$species,] 
map <- fasterize(iucn, r, fun="sum")

fix <- function(x) {sum(x,r,na.rm=TRUE)+r} # This adds zeros for the continental area
map<- fix(map)

#crop the map and clean it up
map<- map %>% 
  crop(c(-170, 170,-90,90)) %>% 
  raster::trim()

#clean up environment
rm(r,clade,iucn)

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

#5 Create Raster #2: Load in human footprint data
setwd("~/Desktop/GitHub/phylofatality/data/footprint/")
footprint <-raster('~/Desktop/GitHub/phylofatality/data/footprint/wildareas-v3-2009-human-footprint.tif')

#Make sure the projections are the same
#Terra package is faster, so first convert both map and footprint RasterLayer --> SpatRaster
footprint <- rast(footprint)
map <- rast(map)
footprint <- terra::project(footprint, map)

#then convert back to RasterLayer
footprint <- raster(footprint)
map <- raster(map)

#these might not be necessary, get rid of outliers and ensure everything is aligned I think
footprint[footprint>50] <- 0
footprint <- footprint + map*0

#check extent, resolution, and projections of rasters (need to be identical)
print(extent(footprint))
print(extent(map))
print(res(footprint))
print(res(map))
print(crs(footprint))
print(crs(map))

#check raster layers individually
my.colors = colorRampPalette(c("#be64ac","lightblue", "yellow","orangered", "red"))
plot(footprint,frame.plot=F,axes=F,box=F,add=F,legend.width=1,legend.shrink=1,col=my.colors(255)) 
plot(map,frame.plot=F,axes=F,box=F,add=F,legend.width=1,legend.shrink=1,col=my.colors(255)) 


#6 Map the bivariate map
bivmap<-bivariate.map(footprint, map, colormatrix=col.matrix, nquantiles=10)

terra:: plot(bivmap, frame.plot = TRUE, axes = F, box = T, add = F, legend = F, col = as.vector(col.matrix), asp = 1)
map(interior = T, add = T)
