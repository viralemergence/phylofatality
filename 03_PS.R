## phylofatality
## 03_PS
## danbeck@ou.edu, carolinecummings2018@gmail.com
## last update 6/19/2025

## clean environment & plots
rm(list=ls()) 
graphics.off()
gc()

## packages
library(ape)
library(caper)
library(dplyr)
#library(plyr)
library(ggtree)
library(ggplot2)
library(data.table)
library(treeio)
library(Hmisc)
library(phylofactor)
library(parallel)
library(emmeans)
library(phytools)
library(egg)
library(devtools)
library(stringr)
#library(ggpubfigs)

## load in virulence data
setwd("~/Desktop/GitHub/phylofatality/csv files")
data=read.csv("01_CFRbySpecies.csv")

## load Upham phylogeny
setwd("~/Desktop/GitHub/phylofatality/phylo")
tree=read.nexus('MamPhy_fullPosterior_BDvr_Completed_5911sp_topoCons_NDexp_MCC_v2_target.tre')

## load in taxonomy
taxa=read.csv('taxonomy_mamPhy_5911species.csv',header=T)
taxa$tip=taxa$Species_Name

## fix tip
tree$tip.label=sapply(strsplit(tree$tip.label,'_'),function(x) paste(x[1],x[2],sep=' '))
taxa$species=sapply(strsplit(taxa$tip,'_'),function(x) paste(x[1],x[2],sep=' '))

## merge data and taxa
data=merge(data,taxa[c("species","tiplabel","gen","fam","ord","clade")],by="species",all.x=T)
rm(taxa)

## trim tree
tree=keep.tip(tree,data$species)

## save data all mammals, all viruses
data$label=data$species
data$Species=data$species
# all mammals, 5 viruses
#coronaviridae
data_cor <- data
data_cor <- dplyr::select(data_cor, species, label, Species, contains("coronaviridae"))
#flaviviridae
data_fla <- data
data_fla <- dplyr::select(data_fla, species, label, Species, contains("flaviviridae"))
#rhabdoviridae
data_rha <- data
data_rha <- dplyr::select(data_rha, species, label, Species, contains("rhabdoviridae"))
#togaviridae
data_tog <- data
data_tog <- dplyr::select(data_tog, species, label, Species, contains("togaviridae"))
#paramyxoviridae
data_par <- data
data_par <- dplyr::select(data_par, species, label, Species, contains("paramyxoviridae"))

##all mammals, all viruses, no onward transmission
data$ntrans_all.viruses=data$virusesWithOT_all.viruses-data$htrans_all.viruses
##all mammals, 5 viruses
#coronaviridae
data_cor$ntrans_cor=data_cor$virusesWithOT_coronaviridae-data_cor$htrans_coronaviridae
#flaviviridiae
data_fla$ntrans_fla=data_fla$virusesWithOT_flaviviridae-data_fla$htrans_flaviviridae
#rhabdoviridae
data_rha$ntrans_rha=data_rha$virusesWithOT_rhabdoviridae-data_rha$htrans_rhabdoviridae
#togaviridae
data_tog$ntrans_tog=data_tog$virusesWithOT_togaviridae-data_tog$htrans_togaviridae
#paramyxoviridae
data_par$ntrans_par=data_par$virusesWithOT_paramyxoviridae-data_par$htrans_paramyxoviridae

## merge: all mammals, all viruses
cdata=comparative.data(phy=tree,data=data,names.col=species,vcv=T,na.omit=F,warn.dropped=T)
#merge: all mammals, 5 viruses
#coronaviridae
cdata_cor=comparative.data(phy=tree,data=data_cor,names.col=species,vcv=T,na.omit=F,warn.dropped=T)
#flaviviridae
cdata_fla=comparative.data(phy=tree,data=data_fla,names.col=species,vcv=T,na.omit=F,warn.dropped=T)
#rhabdoviridae
cdata_rha=comparative.data(phy=tree,data=data_rha,names.col=species,vcv=T,na.omit=F,warn.dropped=T)
#togaviridae
cdata_tog=comparative.data(phy=tree,data=data_tog,names.col=species,vcv=T,na.omit=F,warn.dropped=T)
#paramyxoviridae
cdata_par=comparative.data(phy=tree,data=data_par,names.col=species,vcv=T,na.omit=F,warn.dropped=T)

## taxonomy: all mammals, all viruses
cdata$data$taxonomy=paste(cdata$data$fam,cdata$data$gen,cdata$data$Species,sep='; ')
#taxonomy: all mammals, 5 viruses
#coronaviridae
cdata_cor$data$taxonomy=paste(cdata_cor$data$fam,cdata_cor$data$gen,cdata_cor$data$Species,sep='; ')
#flaviviridae
cdata_fla$data$taxonomy=paste(cdata_fla$data$fam,cdata_fla$data$gen,cdata_fla$data$Species,sep='; ')
#rhabdoviridae
cdata_rha$data$taxonomy=paste(cdata_rha$data$fam,cdata_rha$data$gen,cdata_rha$data$Species,sep='; ')
#togaviridae
cdata_tog$data$taxonomy=paste(cdata_tog$data$fam,cdata_tog$data$gen,cdata_tog$data$Species,sep='; ')
#paramyxoviridae
cdata_par$data$taxonomy=paste(cdata_par$data$fam,cdata_par$data$gen,cdata_par$data$Species,sep='; ')

## separate dataset for onward transmission (remove NA): all mammals, all viruses
##all mammals, all viruses
cdata2=cdata[!is.na(cdata$data$on.frac_all.viruses),]
##all mammals, 5 viruses
#coronaviridae
cdata2_cor=cdata_cor[!is.na(cdata_cor$data$on.frac_coronaviridae),]
#flaviviridae
cdata2_fla=cdata_fla[!is.na(cdata_fla$data$on.frac_flaviviridae),]
#rhabdoviridae
cdata2_rha=cdata_rha[!is.na(cdata_rha$data$on.frac_rhabdoviridae),]
#togaviridae
cdata2_tog=cdata_tog[!is.na(cdata_tog$data$on.frac_togaviridae),]
#paramyxoviridae
cdata2_par=cdata_par[!is.na(cdata_par$data$on.frac_paramyxoviridae),]

#Call mammals, all viruses, using pgls to calculate PS (Pagel's lambda)
mod_me=pgls(meanCFR_all.viruses~1,data=cdata,lambda="ML")
mod_mx=pgls(maxCFR_all.viruses~1,data=cdata,lambda="ML")
mod_ot=pgls(on.frac_all.viruses~1,data=cdata2,lambda="ML")
mod_db=pgls(meanDB_all.viruses~1,data=cdata,lambda="ML")

#all mammals, 5 viruses
#coronaviridae
#mod_me_cor=pgls(meanCFR_coronaviridae~1,data=cdata_cor,lambda="ML") ## see hist(cdata_cor$data$meanCFR_coronaviridae) ... no variation
#mod_mx_cor=pgls(maxCFR_coronaviridae~1,data=cdata_cor,lambda="ML")
#mod_ot_cor=pgls(on.frac_coronaviridae~1,data=cdata2_cor,lambda="ML")
#mod_db_cor=pgls(meanDB_coronaviridae~1,data=cdata_cor,lambda="ML")

#flaviviridae
mod_me_fla=pgls(meanCFR_flaviviridae~1,data=cdata_fla,lambda="ML")
mod_mx_fla=pgls(maxCFR_flaviviridae~1,data=cdata_fla,lambda="ML")
mod_ot_fla=pgls(on.frac_flaviviridae~1,data=cdata2_fla,lambda="ML")
mod_db_fla=pgls(meanDB_flaviviridae~1,data=cdata_fla,lambda="ML")

#rhabdoviridae
mod_me_rha=pgls(meanCFR_rhabdoviridae~1,data=cdata_rha,lambda="ML")
mod_mx_rha=pgls(maxCFR_rhabdoviridae~1,data=cdata_rha,lambda="ML")
#mod_ot_rha=pgls(on.frac_rhabdoviridae~1,data=cdata2_rha,lambda="ML")#see histogram
mod_db_rha=pgls(meanDB_rhabdoviridae~1,data=cdata_rha,lambda="ML")

#togaviridae
mod_me_tog=pgls(meanCFR_togaviridae~1,data=cdata_tog,lambda="ML")
mod_mx_tog=pgls(maxCFR_togaviridae~1,data=cdata_tog,lambda="ML")
mod_ot_tog=pgls(on.frac_togaviridae~1,data=cdata2_tog,lambda="ML")
mod_db_tog=pgls(meanDB_togaviridae~1,data=cdata_tog,lambda="ML")

#paramyxoviridae
mod_me_par=pgls(meanCFR_paramyxoviridae~1,data=cdata_par,lambda="ML")
mod_mx_par=pgls(maxCFR_paramyxoviridae~1,data=cdata_par,lambda="ML")
mod_ot_par=pgls(on.frac_paramyxoviridae~1,data=cdata2_par,lambda="ML")
mod_db_par=pgls(meanDB_paramyxoviridae~1,data=cdata_par,lambda="ML")

## Bloomberg's K
## all mammals, all viruses
psk_me=phylosig(cdata$phy,cdata$data$meanCFR_all.viruses,method="K",test=T)
psk_mx=phylosig(cdata$phy,cdata$data$maxCFR_all.viruses,method="K",test=T)
psk_ot=phylosig(cdata2$phy,cdata2$data$on.frac_all.viruses,method="K",test=T)
psk_db=phylosig(cdata$phy,cdata$data$meanDB_all.viruses,method="K",test=T)

## coronaviridae
psk_me_cor=phylosig(cdata_cor$phy,cdata_cor$data$meanCFR_coronaviridae,method="K",test=T)
psk_mx_cor=phylosig(cdata_cor$phy,cdata_cor$data$maxCFR_coronaviridae,method="K",test=T)
psk_ot_cor=phylosig(cdata2_cor$phy,cdata2_cor$data$on.frac_coronaviridae,method="K",test=T)
psk_db_cor=phylosig(cdata_cor$phy,cdata_cor$data$meanDB_coronaviridae,method="K",test=T)

## flaviviridae
psk_me_fla=phylosig(cdata_fla$phy,cdata_fla$data$meanCFR_flaviviridae,method="K",test=T)
psk_mx_fla=phylosig(cdata_fla$phy,cdata_fla$data$maxCFR_flaviviridae,method="K",test=T)
psk_ot_fla=phylosig(cdata2_fla$phy,cdata2_fla$data$on.frac_flaviviridae,method="K",test=T)
psk_db_fla=phylosig(cdata_fla$phy,cdata_fla$data$meanDB_flaviviridae,method="K",test=T)

## rhabdoviridae
psk_me_rha=phylosig(cdata_rha$phy,cdata_rha$data$meanCFR_rhabdoviridae,method="K",test=T)
psk_mx_rha=phylosig(cdata_rha$phy,cdata_rha$data$maxCFR_rhabdoviridae,method="K",test=T)
#psk_ot_rha=phylosig(cdata2_rha$phy,cdata2_rha$data$on.frac_rhabdoviridae,method="K",test=T) #see histogram
psk_db_rha=phylosig(cdata_rha$phy,cdata_rha$data$meanDB_rhabdoviridae,method="K",test=T)

## togaviridae
psk_me_tog=phylosig(cdata_tog$phy,cdata_tog$data$meanCFR_togaviridae,method="K",test=T)
psk_mx_tog=phylosig(cdata_tog$phy,cdata_tog$data$maxCFR_togaviridae,method="K",test=T)
psk_ot_tog=phylosig(cdata2_tog$phy,cdata2_tog$data$on.frac_togaviridae,method="K",test=T)
psk_db_tog=phylosig(cdata_tog$phy,cdata_tog$data$meanDB_togaviridae,method="K",test=T)

## paramyxoviridae
psk_me_par=phylosig(cdata_par$phy,cdata_par$data$meanCFR_paramyxoviridae,method="K",test=T)
psk_mx_par=phylosig(cdata_par$phy,cdata_par$data$maxCFR_paramyxoviridae,method="K",test=T)
psk_ot_par=phylosig(cdata2_par$phy,cdata2_par$data$on.frac_paramyxoviridae,method="K",test=T)
psk_db_par=phylosig(cdata_par$phy,cdata_par$data$meanDB_paramyxoviridae,method="K",test=T)

## bat analyses: subset to bats, all bats, all viruses
bdata=cdata[cdata$data$ord=="CHIROPTERA",]
bdata2=cdata2[cdata2$data$ord=="CHIROPTERA",]

#all bats, 5 viruses
#coronaviridae
bdata_cor <- bdata
bdata_cor$data <- dplyr::select(bdata_cor$data, meanCFR_coronaviridae, maxCFR_coronaviridae, virusesWithCFR_coronaviridae,
                                meanDB_coronaviridae, htrans_coronaviridae, on.frac_coronaviridae, virusesWithOT_coronaviridae)
#flaviviridae
bdata_fla <- bdata
bdata_fla$data <- dplyr::select(bdata_fla$data, meanCFR_flaviviridae, maxCFR_flaviviridae, virusesWithCFR_flaviviridae,
                                meanDB_flaviviridae, htrans_flaviviridae, on.frac_flaviviridae, virusesWithOT_flaviviridae)
#rhabdoviridae
bdata_rha <- bdata
bdata_rha$data <- dplyr::select(bdata_rha$data, meanCFR_rhabdoviridae, maxCFR_rhabdoviridae, virusesWithCFR_rhabdoviridae,
                                meanDB_rhabdoviridae, htrans_rhabdoviridae, on.frac_rhabdoviridae, virusesWithOT_rhabdoviridae)
#togaviridae
bdata_tog <- bdata
bdata_tog$data <- dplyr::select(bdata_tog$data, meanCFR_togaviridae, maxCFR_togaviridae, virusesWithCFR_togaviridae,
                                meanDB_togaviridae, htrans_togaviridae, on.frac_togaviridae, virusesWithOT_togaviridae)
#paramyxoviridae
bdata_par <- bdata
bdata_par$data <- dplyr::select(bdata_par$data, meanCFR_paramyxoviridae, maxCFR_paramyxoviridae, virusesWithCFR_paramyxoviridae,
                                meanDB_paramyxoviridae, htrans_paramyxoviridae, on.frac_paramyxoviridae, virusesWithOT_paramyxoviridae)

##all bats, all viruses
bdata2=bdata[!is.na(bdata$data$on.frac_all.viruses),]
#coronaviridae
bdata2_cor=bdata_cor[!is.na(bdata_cor$data$on.frac_coronaviridae),]
#flaviviridae
bdata2_fla=bdata_fla[!is.na(bdata_fla$data$on.frac_flaviviridae),]
#rhabdoviridae
bdata2_rha=bdata_rha[!is.na(bdata_rha$data$on.frac_rhabdoviridae),]
#togaviridae
bdata2_tog=bdata_tog[!is.na(bdata_tog$data$on.frac_togaviridae),]
#paramyxoviridae
bdata2_par=bdata_par[!is.na(bdata_par$data$on.frac_paramyxoviridae),]

#all bats, all viruses pagel's lambda
bmod_me=pgls(meanCFR_all.viruses~1,data=bdata,lambda="ML")
bmod_mx=pgls(maxCFR_all.viruses~1,data=bdata,lambda="ML")
bmod_ot=pgls(on.frac_all.viruses~1,data=bdata2,lambda="ML")
bmod_db=pgls(meanDB_all.viruses~1,data=bdata,lambda="ML")

##all bats, 5 virus families
#coronaviridae
# bmod_me_cor=pgls(meanCFR_coronaviridae~1,data=bdata_cor,lambda="ML")
# bmod_mx_cor=pgls(maxCFR_coronaviridae~1,data=bdata_cor,lambda="ML")
# bmod_ot_cor_bat=pgls(on.frac_coronaviridae~1,data=bdata2_cor,lambda="ML") #see histogram
# bmod_db_cor=pgls(meanDB_coronaviridae~1,data=bdata_cor,lambda="ML")

#flaviviridae
bmod_me_fla=pgls(meanCFR_flaviviridae~1,data=bdata_fla,lambda="ML")
bmod_mx_fla=pgls(maxCFR_flaviviridae~1,data=bdata_fla,lambda="ML")
bmod_ot_fla=pgls(on.frac_flaviviridae~1,data=bdata2_fla,lambda="ML")
bmod_db_fla=pgls(meanDB_flaviviridae~1,data=bdata_fla,lambda="ML")

#rhabdoviridae
bmod_me_rha=pgls(meanCFR_rhabdoviridae~1,data=bdata_rha,lambda="ML")
bmod_mx_rha=pgls(maxCFR_rhabdoviridae~1,data=bdata_rha,lambda="ML")
#bmod_ot_rha=pgls(on.frac_rhabdoviridae~1,data=bdata2_rha,lambda="ML") #see histogram
bmod_db_rha=pgls(meanDB_rhabdoviridae~1,data=bdata_rha,lambda="ML")

#togaviridae
bmod_me_tog=pgls(meanCFR_togaviridae~1,data=bdata_tog,lambda="ML")
bmod_mx_tog=pgls(maxCFR_togaviridae~1,data=bdata_tog,lambda="ML")
bmod_ot_tog=pgls(on.frac_togaviridae~1,data=bdata2_tog,lambda="ML")
bmod_db_tog=pgls(meanDB_togaviridae~1,data=bdata_tog,lambda="ML")

#paramyxoviridae
bmod_me_par=pgls(meanCFR_paramyxoviridae~1,data=bdata_par,lambda="ML")
bmod_mx_par=pgls(maxCFR_paramyxoviridae~1,data=bdata_par,lambda="ML")
bmod_ot_par=pgls(on.frac_paramyxoviridae~1,data=bdata2_par,lambda="ML")
bmod_db_par=pgls(meanDB_paramyxoviridae~1,data=bdata_par,lambda="ML")

## Bloomberg's K: bats
## bats only, all viruses
bpsk_me=phylosig(bdata$phy,bdata$data$meanCFR_all.viruses,method="K",test=T)
bpsk_mx=phylosig(bdata$phy,bdata$data$maxCFR_all.viruses,method="K",test=T)
bpsk_ot=phylosig(bdata2$phy,bdata2$data$on.frac_all.viruses,method="K",test=T)
bpsk_db=phylosig(bdata$phy,bdata$data$meanDB_all.viruses,method="K",test=T)

## coronaviridae
bpsk_me_cor=phylosig(bdata_cor$phy,bdata_cor$data$meanCFR_coronaviridae,method="K",test=T)
bpsk_mx_cor=phylosig(bdata_cor$phy,bdata_cor$data$maxCFR_coronaviridae,method="K",test=T)
bpsk_ot_cor=phylosig(bdata2_cor$phy,bdata2_cor$data$on.frac_coronaviridae,method="K",test=T)
bpsk_db_cor=phylosig(bdata_cor$phy,bdata_cor$data$meanDB_coronaviridae,method="K",test=T)

## flaviviridae
bpsk_me_fla=phylosig(bdata_fla$phy,bdata_fla$data$meanCFR_flaviviridae,method="K",test=T)
bpsk_mx_fla=phylosig(bdata_fla$phy,bdata_fla$data$maxCFR_flaviviridae,method="K",test=T)
bpsk_ot_fla=phylosig(bdata2_fla$phy,bdata2_fla$data$on.frac_flaviviridae,method="K",test=T)
bpsk_db_fla=phylosig(bdata_fla$phy,bdata_fla$data$meanDB_flaviviridae,method="K",test=T)

##rhabdoviridae
bpsk_me_rha=phylosig(bdata_rha$phy,bdata_rha$data$meanCFR_rhabdoviridae,method="K",test=T)
bpsk_mx_rha=phylosig(bdata_rha$phy,bdata_rha$data$maxCFR_rhabdoviridae,method="K",test=T)
#bpsk_ot_rha=phylosig(bdata2_rha$phy,bdata2_rha$data$on.frac,method="K",test=T) #see histogram
bpsk_db_rha=phylosig(bdata_rha$phy,bdata_rha$data$meanDB_rhabdoviridae,method="K",test=T)

## togaviridae
bpsk_me_tog=phylosig(bdata_tog$phy,bdata_tog$data$meanCFR_togaviridae,method="K",test=T)
bpsk_mx_tog=phylosig(bdata_tog$phy,bdata_tog$data$maxCFR_togaviridae,method="K",test=T)
bpsk_ot_tog=phylosig(bdata2_tog$phy,bdata2_tog$data$on.frac_togaviridae,method="K",test=T)
bpsk_db_tog=phylosig(bdata_tog$phy,bdata_tog$data$meanDB_togaviridae,method="K",test=T)

## paramyxoviridae
bpsk_me_par=phylosig(bdata_par$phy,bdata_par$data$meanCFR_paramyxoviridae,method="K",test=T)
bpsk_mx_par=phylosig(bdata_par$phy,bdata_par$data$maxCFR_paramyxoviridae,method="K",test=T)
bpsk_ot_par=phylosig(bdata2_par$phy,bdata2_par$data$on.frac_paramyxoviridae,method="K",test=T)
bpsk_db_par=phylosig(bdata_par$phy,bdata_par$data$meanDB_paramyxoviridae,method="K",test=T)

## summarize Lambda estimates: mammals, bats, all viruses
mlist=list(mod_me,mod_mx,mod_ot,mod_db,bmod_me,bmod_mx,bmod_ot,bmod_db)
pdata=data.frame(vfamily=rep("all viruses",8),
                 dataset=c(rep("all mammals",4),rep("bats only",4)),
                 variable=c(rep(c("meanCFR","maxCFR","on.frac","meanDB"),2)),
                 lambda=sapply(mlist,function(x) x$param["lambda"]),
                 lambda_lower=sapply(mlist,function(x) x$param.CI$lambda$ci.val[1]),
                 lambda_lower_p=sapply(mlist,function(x) x$param.CI$lambda$bounds.p[1]),
                 lambda_upper=sapply(mlist,function(x) x$param.CI$lambda$ci.val[2]),
                 lambda_upper_p=sapply(mlist,function(x) x$param.CI$lambda$bounds.val[1]))
pdata$variable=factor(pdata$variable,levels=c("meanCFR","maxCFR","on.frac","meanDB"))

## mammals, bats, coronaviridae
# mlist_cor=list(mod_me_cor,mod_mx_cor,mod_db_cor,bmod_me_cor,bmod_mx_cor,bmod_db_cor)
# pdata_cor=data.frame(vfamily=rep("coronaviridae",6),
#                      dataset=c(rep("all mammals",3),rep("bats only",3)),
#                      variable=c(rep(c("meanCFR","maxCFR","meanDB"),2)),
#                      lambda=sapply(mlist_cor,function(x) x$param["lambda"]),
#                      lambda_lower=sapply(mlist_cor,function(x) x$param.CI$lambda$ci.val[1]),
#                      lambda_lower_p=sapply(mlist_cor,function(x) x$param.CI$lambda$bounds.p[1]),
#                      lambda_upper=sapply(mlist_cor,function(x) x$param.CI$lambda$ci.val[2]),
#                      lambda_upper_p=sapply(mlist_cor,function(x) x$param.CI$lambda$bounds.val[1]))
# pdata_cor$variable=factor(pdata_cor$variable,levels=c("meanCFR","maxCFR","meanDB"))

#flaviviridae
mlist_fla=list(mod_me_fla,mod_mx_fla,mod_ot_fla,mod_db_fla,bmod_me_fla,bmod_mx_fla,bmod_ot_fla,bmod_db_fla)
pdata_fla=data.frame(vfamily=rep("flaviviridae",8),
                     dataset=c(rep("all mammals",4),rep("bats only",4)),
                     variable=c(rep(c("meanCFR","maxCFR","on.frac","meanDB"),2)),
                     lambda=sapply(mlist_fla,function(x) x$param["lambda"]),
                     lambda_lower=sapply(mlist_fla,function(x) x$param.CI$lambda$ci.val[1]),
                     lambda_lower_p=sapply(mlist_fla,function(x) x$param.CI$lambda$bounds.p[1]),
                     lambda_upper=sapply(mlist_fla,function(x) x$param.CI$lambda$ci.val[2]),
                     lambda_upper_p=sapply(mlist_fla,function(x) x$param.CI$lambda$bounds.val[1]))
pdata_fla$variable=factor(pdata_fla$variable,levels=c("meanCFR","maxCFR","on.frac","meanDB"))

#rhabdoviridae
mlist_rha=list(mod_me_rha,mod_mx_rha,mod_db_rha,bmod_me_rha,bmod_mx_rha,bmod_db_rha)
pdata_rha=data.frame(vfamily=rep("rhabdoviridae",6),
                     dataset=c(rep("all mammals",3),rep("bats only",3)),
                     variable=c(rep(c("meanCFR","maxCFR","meanDB"),2)),
                     lambda=sapply(mlist_rha,function(x) x$param["lambda"]),
                     lambda_lower=sapply(mlist_rha,function(x) x$param.CI$lambda$ci.val[1]),
                     lambda_lower_p=sapply(mlist_rha,function(x) x$param.CI$lambda$bounds.p[1]),
                     lambda_upper=sapply(mlist_rha,function(x) x$param.CI$lambda$ci.val[2]),
                     lambda_upper_p=sapply(mlist_rha,function(x) x$param.CI$lambda$bounds.val[1]))
pdata_rha$variable=factor(pdata_rha$variable,levels=c("meanCFR","maxCFR","meanDB"))

#togaviridae
mlist_tog=list(mod_me_tog,mod_mx_tog,mod_ot_tog,mod_db_tog,bmod_me_tog,bmod_mx_tog,bmod_ot_tog,bmod_db_tog)
pdata_tog=data.frame(vfamily=rep("togaviridae",8),
                     dataset=c(rep("all mammals",4),rep("bats only",4)),
                     variable=c(rep(c("meanCFR","maxCFR","on.frac","meanDB"),2)),
                     lambda=sapply(mlist_tog,function(x) x$param["lambda"]),
                     lambda_lower=sapply(mlist_tog,function(x) x$param.CI$lambda$ci.val[1]),
                     lambda_lower_p=sapply(mlist_tog,function(x) x$param.CI$lambda$bounds.p[1]),
                     lambda_upper=sapply(mlist_tog,function(x) x$param.CI$lambda$ci.val[2]),
                     lambda_upper_p=sapply(mlist_tog,function(x) x$param.CI$lambda$bounds.val[1]))
pdata_tog$variable=factor(pdata_tog$variable,levels=c("meanCFR","maxCFR","on.frac","meanDB"))

#paramyxoviridae
mlist_par=list(mod_me_par, mod_mx_par, mod_ot_par, mod_db_par, bmod_me_par, bmod_mx_par, bmod_ot_par, bmod_db_par)
pdata_par=data.frame(vfamily=rep("paramyxoviridae",8),
                     dataset=c(rep("all mammals",4),rep("bats only",4)),
                     variable=c(rep(c("meanCFR","maxCFR","on.frac","meanDB"),2)),
                     lambda=sapply(mlist_par,function(x) x$param["lambda"]),
                     lambda_lower=sapply(mlist_par,function(x) x$param.CI$lambda$ci.val[1]),
                     lambda_lower_p=sapply(mlist_par,function(x) x$param.CI$lambda$bounds.p[1]),
                     lambda_upper=sapply(mlist_par,function(x) x$param.CI$lambda$ci.val[2]),
                     lambda_upper_p=sapply(mlist_par,function(x) x$param.CI$lambda$bounds.val[1]))
pdata_par$variable=factor(pdata_par$variable,levels=c("meanCFR","maxCFR","on.frac","meanDB"))


#save
pagel<- rbind(pdata, pdata_fla, pdata_rha, pdata_tog, pdata_par)
setwd("~/Desktop/GitHub/phylofatality/csv files")
#write.csv(pagel,"03_PS data.csv")
write.csv(pagel,"03_PS data_20250609.csv")

## summarize bloomberg's K
klist=list(psk_me,psk_mx,psk_ot,psk_db,bpsk_me,bpsk_mx,bpsk_ot, bpsk_db)
kdata=data.frame(vfamily=rep("all viruses",8),
                 dataset=c(rep("all mammals",4), rep("bats only", 4)),
                 variable=c(rep(c("meanCFR", "maxCFR", "on.frac","meanDB"), 2)),
                 K=sapply(klist,function(x) x$"K"),
                 P=sapply(klist,function(x) x$"P"))
kdata$variable=factor(kdata$variable,levels=c("meanCFR","maxCFR","on.frac","meanDB"))

# ## mammals, bats, coronaviridae
# klist=list(psk_me_cor,psk_mx_cor,psk_ot_cor,psk_db_cor,bpsk_me_cor,bpsk_mx_cor,bpsk_ot_cor,bpsk_db_cor)
# kdata_cor=data.frame(vfamily=rep("coronaviridae",8), 
#                      dataset=c(rep("all mammals",4), rep("bats only", 4)),
#                  variable=c(rep(c("meanCFR", "maxCFR", "on.frac","meanDB"), 2)),
#                  K=sapply(klist,function(x) x$"K"),
#                  P=sapply(klist,function(x) x$"P"))
# kdata_cor$variable=factor(kdata_cor$variable,levels=c("meanCFR","maxCFR","on.frac","meanDB"))

## flaviviridae
klist=list(psk_me_fla,psk_mx_fla,psk_ot_fla,psk_db_fla,bpsk_me_fla,bpsk_mx_fla,bpsk_ot_fla, bpsk_db_fla)
kdata_fla=data.frame(vfamily=rep("flaviviridae",8),
                 dataset=c(rep("all mammals",4), rep("bats only", 4)),
                 variable=c(rep(c("meanCFR", "maxCFR", "on.frac","meanDB"), 2)),
                 K=sapply(klist,function(x) x$"K"),
                 P=sapply(klist,function(x) x$"P"))
kdata_fla$variable=factor(kdata$variable,levels=c("meanCFR","maxCFR","on.frac","meanDB"))

## rhabdoviridae
klist=list(psk_me_rha,psk_mx_rha,psk_db_rha,bpsk_me_rha,bpsk_mx_rha, bpsk_db_rha)
kdata_rha=data.frame(vfamily=rep("rhabdoviridae",6),
                 dataset=c(rep("all mammals",3), rep("bats only", 3)),
                 variable=c(rep(c("meanCFR", "maxCFR","meanDB"), 2)),
                 K=sapply(klist,function(x) x$"K"),
                 P=sapply(klist,function(x) x$"P"))
kdata_rha$variable=factor(kdata_rha$variable,levels=c("meanCFR","maxCFR","meanDB"))

## togaviridae
klist=list(psk_me_tog,psk_mx_tog,psk_ot_tog,psk_db_tog,bpsk_me_tog,bpsk_mx_tog,bpsk_ot_tog,bpsk_db_tog)
kdata_tog=data.frame(vfamily=rep("togaviridae",8),
                     dataset=c(rep("all mammals",4), rep("bats only", 4)),
                    variable=c(rep(c("meanCFR", "maxCFR", "on.frac","meanDB"), 2)),
                     K=sapply(klist,function(x) x$"K"),
                     P=sapply(klist,function(x) x$"P"))
kdata_tog$variable=factor(kdata_tog$variable,levels=c("meanCFR","maxCFR","on.frac","meanDB"))

## paramyxoviridae
klist=list(psk_me_par,psk_mx_par,psk_ot_par,psk_db_par,bpsk_me_par,bpsk_mx_par,bpsk_ot_par,bpsk_db_par)
kdata_par=data.frame(vfamily=rep("paramyxoviridae",8),
                     dataset=c(rep("all mammals",4), rep("bats only", 4)),
                     variable=c(rep(c("meanCFR", "maxCFR", "on.frac","meanDB"), 2)),
                     K=sapply(klist,function(x) x$"K"),
                     P=sapply(klist,function(x) x$"P"))
kdata_par$variable=factor(kdata_par$variable,levels=c("meanCFR","maxCFR","on.frac","meanDB"))

## save
bloombergk<- rbind(kdata, kdata_fla, kdata_rha, kdata_tog, kdata_par)
setwd("~/Desktop/GitHub/phylofatality/csv files")
#write.csv(bloombergk,"03_K data.csv")
write.csv(bloombergk,"03_K data_20250609.csv")


#plotting (can start here and reload in data)
#don't forget to reload in packages
setwd("~/Desktop/GitHub/phylofatality/csv files")
#ps=read.csv("03_PS data.csv")
ps=read.csv("03_PS data_20250609.csv")
ps$X=NULL

#fix up variable names
ps$vfamily <- str_to_title(ps$vfamily)
ps$variable<- as.character(ps$variable)

ps$variable<- ifelse(ps$variable=="meanCFR", "mean CFR", ps$variable)
ps$variable<- ifelse(ps$variable=="maxCFR", "maximum CFR", ps$variable)
ps$variable<- ifelse(ps$variable=="on.frac", "% with onward transmission", ps$variable)
ps$variable<- ifelse(ps$variable=="meanDB", "mean death burden", ps$variable)
ps$vfamily <- ifelse(ps$vfamily=="All Viruses", "all viruses", ps$vfamily)

## remove maxCFR
raw<- ps

#reorder 
ps$variable <- factor(ps$variable, levels = c("mean CFR",
                                              "% with onward transmission",
                                              "mean death burden",
                                              "maximum CFR"))

ps <- ps %>% filter(variable != "maximum CFR")

#plot
#use ggpubfigs colorblind friendly colors "nickel_five"
plot <- ggplot(ps, aes(vfamily, lambda, color = vfamily)) +
  theme_bw() +
  facet_grid(variable~dataset)+ ##get rid of Mac
  theme(strip.text = element_text(size = 12))+
  theme(legend.position = "none")+
  theme(axis.text.y = element_text(size = 10))+
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1,
                                   size=11, 
                                   color=c("black", "#648FFF", "#FE6100", "#785EF0","#DC267F")))+
  guides(color = guide_legend(title = "virus family"))+
  geom_errorbar(
    aes(ymin = lambda_lower, ymax = lambda_upper),
    position = position_dodge(width = 0.2), 
    width = 0, 
    size = 1)+
  geom_point(position = position_dodge(width = 0.2), size = 3) +
  ylim(0, 1) +
  scale_color_manual(values = c("black", "#648FFF", "#FE6100", "#785EF0","#DC267F")) +
  labs(x = "Virus Family", y = expression(paste("Pagel's ", lambda)))+
  theme(axis.title.x = element_text(size = 16, margin = margin(t = 20)))+
  theme(axis.title.y = element_text(size = 16, margin = margin(r = 18)))+
  scale_x_discrete(labels=c("all viruses",
                            expression(italic(Flaviviridae)),
                            expression(italic(Paramyxoviridae)),
                            expression(italic(Rhabdoviridae)),
                            expression(italic(Togaviridae))))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
plot(plot)

#save
setwd("~/Desktop/GitHub/phylofatality/figs")
ggsave("fig1_20250609.jpg", plot, device = "jpeg", width = 7, height =10, units = "in")
