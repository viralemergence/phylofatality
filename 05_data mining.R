## phylofatality
## 05_data mining (species extraction from risky clades)
## danbeck@ou.edu, carolinecummings2018@gmail.com, Cole Brookson
## last update 6/10/2025

## clean environment & plots
rm(list=ls()) 
graphics.off()
gc()

## packages
library(dplyr)
library(plyr)
library(tidyr)
library(tidyverse)

## load in clade virulence data
setwd("~/Desktop/GitHub/phylofatality/csv files")
#data=read.csv("04_pf_allclades.csv")
data=read.csv("04_pf_allclades_20250609.csv")

## load in host taxonomy
setwd("~/Desktop/GitHub/phylofatality/phylo")
taxa=read.csv('taxonomy_mamPhy_5911species.csv',header=T)
taxa$tip=taxa$Species_Name

#label risky clades
data$risk=ifelse(data$clade>data$other, "risky", "non-risky")
rawdata=data

#filter to risk clades
data<- data%>% filter(risk=="risky")
data$risk=NULL

#separate variable and virus
data=data %>% separate(ID, c('var', 'virus'))

#add all viruses to replace NAs
data[is.na(data)]<- "all"

#separate var into host and var
data=data %>% separate(var, c('host', 'var'), sep=1)

#rename hosts
data$host=ifelse(data$host=="c", "mammal", data$host)
data$host=ifelse(data$host=="d", "mammal", data$host)
data$host=ifelse(data$host=="b", "bat", data$host)


#reorder and clean up table
data=data %>% dplyr::select(virus, host, var, everything())
#data$var=revalue(data$var,c("means"= "mean"))
data[which(data$var == "means"), "var"] <- "mean"
data[which(data$var == "bmean"), "var"] <- "db"
data[which(data$var == "dbmean"), "var"] <- "db"

#save data of risky clades
rawdata_risk<-data

#pull out species in each clade
data$species=data$taxa
data=data %>% dplyr::select(species, everything())
data=data %>% separate_rows(species, sep = ", ")

## Cole helped match the fams and gens to species
`%notin%` <- Negate(`%in%`)
data_fam <- data %>% 
  # pick out the family one by only keeping the ones that are all uppercase
  dplyr::filter(grepl("^[[:upper:]]+$", species)) %>% 
  # rename for easier joining
  dplyr::rename(fam = species)
data_gen <- data %>% 
  dplyr::filter(species %notin% data_fam$fam) %>% 
  dplyr::rename(gen = species)

# now join one at a time - use genera preferentially
joined_gen <- dplyr::full_join(
  taxa,
  data_gen,
  by = "gen"
)

# now filter the ones that got the data from genera and don't need family
which_taxa_for_fam <- joined_gen %>% 
  dplyr::filter_at(vars(virus, host, var, taxa, tips), 
                   all_vars(is.na(.)))

taxa_fam <- taxa %>% 
  dplyr::filter(fam %in% which_taxa_for_fam$fam)

joined_fam <- dplyr::full_join(
  taxa_fam, 
  data_fam,
  by = "fam"
) 

# now remove the ones from joined_gen who had those NA's 
joined_gen <- joined_gen %>% 
  dplyr::filter_at(vars(virus, host, var, taxa, tips), 
                   all_vars(!is.na(.)))

# remove the ones in joined_fam already accounted for in joined_gen
joined_fam <- joined_fam %>% 
  dplyr::filter(gen %notin% joined_gen$gen)

# bind together
all_joined <- rbind(joined_fam, joined_gen)

#remove NAs
species=all_joined[!is.na(all_joined$taxa),]
raw<- species

#clean up table
species=species %>% dplyr::select(Species_Name, virus, host, var, factor, tips, node,
                           clade.y, other, taxa)
species$species=species$Species_Name
species$clade=species$clade.y
species$Species_Name=NULL
species$clade.y=NULL
species=species %>% dplyr::select(species, virus, host, var, factor, tips, node,
                           clade, other, taxa)
#sanity check and save
species %>% n_distinct() ## 19,359
setwd("~/Desktop/GitHub/phylofatality/csv files")
#write.csv(species,"05_pf_riskyspecies.csv")
write.csv(species,"05_pf_riskyspecies_20250609.csv")

## pull out the risky bat clades to see what they look like
ord <- taxa %>% select(Species_Name, ord)
ord$species <- ord$Species_Name
ord$Species_Name=NULL

bats<- merge(ord,raw)

bats <- bats %>% filter(ord=="CHIROPTERA")
bats$ID<- paste0(bats$virus,"_",bats$host,"_",bats$var,"_",bats$factor)
bats<- unique(bats$ID)
bats <- bats %>% as.data.frame()
