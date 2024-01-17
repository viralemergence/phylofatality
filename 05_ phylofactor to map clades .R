## phylofatality
##05_phylofactor to map clades
## danbeck@ou.edu carolinecummings@ou.edu
## last update 12/18/2023

## clean environment & plots
rm(list=ls()) 
graphics.off()

## packages
library(ape)
library(caper)
library(plyr)
library(ggtree)
library(ggplot2)
library(data.table)
library(treeio)
library(Hmisc)
library(phylofactor)
library(parallel)
library(emmeans)
#library(ade4)
library(phytools)
library(dplyr)

## load in virulence data
setwd("~/Desktop/PCM Class/phylofatality")
data=read.csv("CFRbySpecies.csv")

## load Upham phylogeny
setwd("~/Desktop/PCM Class/phylofatality/phylo")
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

## save
data$label=data$species
data$Species=data$species

## define non-onward viruses
data$ntrans_all.viruses=data$virusesWithOT_all.viruses-data$htrans_all.viruses
#5 virus families
data$ntrans_coronaviridae=data$virusesWithOT_coronaviridae-data$htrans_coronaviridae
data$ntrans_flaviviridae=data$virusesWithOT_flaviviridae-data$htrans_flaviviridae
data$ntrans_rhabdoviridae=data$virusesWithOT_rhabdoviridae-data$htrans_rhabdoviridae
data$ntrans_togaviridae=data$virusesWithOT_togaviridae-data$htrans_togaviridae
data$ntrans_paramyxoviridae=data$virusesWithOT_paramyxoviridae-data$htrans_paramyxoviridae

## merge
cdata=comparative.data(phy=tree,data=data,names.col=species,vcv=T,na.omit=F,warn.dropped=T)

## taxonomy
cdata$data$taxonomy=paste(cdata$data$fam,cdata$data$gen,cdata$data$Species,sep='; ')
## set taxonomy
taxonomy=data.frame(cdata$data$taxonomy)
names(taxonomy)="taxonomy"
taxonomy$Species=rownames(cdata$data)
taxonomy=taxonomy[c("Species","taxonomy")]
taxonomy$taxonomy=as.character(taxonomy$taxonomy)

## separate dataset for onward transmission (remove NA)
cdata2=cdata[!is.na(cdata$data$on.frac_all.viruses),]
#5 virus families
cdata2_cov<-cdata[!is.na(cdata$data$on.frac_coronaviridae),]
cdata2_fla<-cdata[!is.na(cdata$data$on.frac_flaviviridae),]
cdata2_rha<-cdata[!is.na(cdata$data$on.frac_rhabdoviridae),]
cdata2_tog<-cdata[!is.na(cdata$data$on.frac_togaviridae),]
cdata2_par<-cdata[!is.na(cdata$data$on.frac_paramyxoviridae),]

## subset to bats
bdata=cdata[cdata$data$ord=="CHIROPTERA",]
bdata2=cdata2[cdata2$data$ord=="CHIROPTERA",]

#save dataframes for each virus family for meanCFR and maxCFR
cdata_cov<-cdata[!is.na(cdata$data$meanCFR_coronaviridae),]
cdata_fla<-cdata[!is.na(cdata$data$meanCFR_flaviviridae),]
cdata_rha<-cdata[!is.na(cdata$data$meanCFR_rhabdoviridae),]
cdata_tog<-cdata[!is.na(cdata$data$meanCFR_togaviridae),]
cdata_par<-cdata[!is.na(cdata$data$meanCFR_paramyxoviridae),]

## Holm rejection procedure
HolmProcedure <- function(pf,FWER=0.05){
  
  ## get split variable
  cs=names(coef(pf$models[[1]]))[-1]
  split=ifelse(length(cs)>1,cs[3],cs[1])
  split=ifelse(is.na(split),cs[1],split)
  
  ## obtain p values
  if (pf$models[[1]]$family$family%in%c('gaussian',"Gamma","quasipoisson")){
    pvals <- sapply(pf$models,FUN=function(fit) summary(fit)$coefficients[split,'Pr(>|t|)'])
  } else {
    pvals <- sapply(pf$models,FUN=function(fit) summary(fit)$coefficients[split,'Pr(>|z|)'])
  }
  D <- length(pf$tree$tip.label)
  
  ## this is the line for Holm's sequentially rejective cutoff
  keepers <- pvals<=(FWER/(2*D-3 - 2*(0:(pf$nfactors-1))))
  
  
  if (!all(keepers)){
    nfactors <- min(which(!keepers))-1
  } else {
    nfactors <- pf$nfactors
  }
  return(nfactors)
}

## get species in a clade
cladeget=function(pf,factor){
  spp=pf$tree$tip.label[pf$groups[[factor]][[1]]]
  return(spp)
}

## summarize pf object 
pfsum=function(pf){
  
  ## get formula
  chars=as.character(pf$frmla.phylo)[-1]
  
  ## response
  resp=chars[1]
  
  ## holm
  hp=HolmProcedure(pf)
  
  ## save model
  model=chars[2]
  
  ## set key
  setkey(pf$Data,'Species')
  
  ## make data
  dat=data.frame(pf$Data)
  
  ## make clade columns in data
  for(i in 1:hp){
    
    dat[,paste0(resp,'_pf',i)]=ifelse(dat$Species%in%cladeget(pf,i),'factor','other')
    
  }
  
  ## make data frame to store taxa name, response, mean, and other
  results=data.frame(matrix(ncol=6, nrow = hp))
  colnames(results)=c('factor','taxa','tips','node',"clade",'other')
  
  ## set taxonomy
  taxonomy=dat[c('Species','taxonomy')]
  taxonomy$taxonomy=as.character(taxonomy$taxonomy)
  
  ## loop
  for(i in 1:hp){
    
    ## get taxa
    tx=pf.taxa(pf,taxonomy,factor=i)$group1
    
    ## get tail
    tx=sapply(strsplit(tx,'; '),function(x) tail(x,1))
    
    ## combine
    tx=paste(tx,collapse=', ')
    
    # save
    results[i,'factor']=i
    results[i,'taxa']=tx
    
    ## get node
    tips=cladeget(pf,i)
    node=ggtree::MRCA(pf$tree,tips)
    results[i,'tips']=length(tips)
    results[i,'node']=ifelse(is.null(node) & length(tips)==1,'species',
                             ifelse(is.null(node) & length(tips)!=1,NA,node))
    
    ## get glm 
    mod=pf$models[[i]]
    mod=glm(formula(mod),data=mod$data,family=summary(mod)$family$family)
    em=data.frame(emmeans(mod,"phylo",type="response"))
    
    ## add category
    em$cat=revalue(em$phylo,
                   c("R"="factor",
                     "S"="other"))
    
    ## get means
    #ms=(tapply(dat[,resp],dat[,paste0(resp,'_pf',i)],mean))
    
    ## add in
    #results[i,'clade']=ms['factor']
    #results[i,'other']=ms['other']
    results[i,'clade']=em[em$cat=="factor",2]
    results[i,'other']=em[em$cat=="other",2]
  }
  
  ## return
  return(list(set=dat,results=results))
}

## NOTE THAT gpf() WON'T WORK UNDER R 4.2 OR HIGHER DUE TO THE FOLLOWING CHANGE:
## https://stackoverflow.com/questions/72848442/r-warning-lengthx-2-1-in-coercion-to-logical1/72848495#72848495

## all viruses meanCFR mammals
set.seed(1)
cmean_pf=gpf(Data=cdata$data,tree=cdata$phy,
             frmla.phylo=meanCFR_all.viruses~phylo+virusesWithCFR_all.viruses,
             family=gaussian,algorithm='phylo',nfactors=5,min.group.size=10)
HolmProcedure(cmean_pf) #4

#coronaviridae meanCFR mammals
set.seed(1)
cmean_pf_cov=gpf(Data=cdata_cov$data,tree=cdata_cov$phy,
                 frmla.phylo=meanCFR_coronaviridae~phylo+virusesWithCFR_coronaviridae,
                 family=gaussian,algorithm='phylo',nfactors=5,min.group.size=10)
HolmProcedure(cmean_pf_cov) #2

#flaviviridae meanCFR mammals
set.seed(1)
cmean_pf_fla=gpf(Data=cdata_fla$data,tree=cdata_fla$phy,
                 frmla.phylo=meanCFR_flaviviridae~phylo+virusesWithCFR_flaviviridae,
                 family=gaussian,algorithm='phylo',nfactors=5,min.group.size=10)
HolmProcedure(cmean_pf_fla) #4

#rhabdoviridae meanCFR mammals
set.seed(1)
cmean_pf_rha=gpf(Data=cdata_rha$data,tree=cdata_rha$phy,
                 frmla.phylo=meanCFR_rhabdoviridae~phylo+virusesWithCFR_rhabdoviridae,
                 family=gaussian,algorithm='phylo',nfactors=5,min.group.size=10)
HolmProcedure(cmean_pf_rha) #2

#togaviridae meanCFR mammals
set.seed(1)
cmean_pf_tog=gpf(Data=cdata_tog$data,tree=cdata_tog$phy,
                 frmla.phylo=meanCFR_togaviridae~phylo+virusesWithCFR_togaviridae,
                 family=gaussian,algorithm='phylo',nfactors=5,min.group.size=10)
HolmProcedure(cmean_pf_tog) #1

#paramyxoviridae meanCFR mammals
set.seed(1)
cmean_pf_par=gpf(Data=cdata_par$data,tree=cdata_par$phy,
                 frmla.phylo=meanCFR_paramyxoviridae~phylo+virusesWithCFR_paramyxoviridae,
                 family=gaussian,algorithm='phylo',nfactors=5,min.group.size=10)
HolmProcedure(cmean_pf_par) #0 

## summarize mammals, all viruses, 5 virus families
#save what each of the clades are for each vfamily
cmean_pf_results=pfsum(cmean_pf)$results
cmean_pf_results_cov=pfsum(cmean_pf_cov)$results #4
cmean_pf_results_fla=pfsum(cmean_pf_fla)$results #2
cmean_pf_results_rha=pfsum(cmean_pf_rha)$results #4
cmean_pf_results_tog=pfsum(cmean_pf_tog)$results #2
#cmean_pf_results_par=pfsum(cmean_pf_par)$results #0 factors

## all viruses maxCFR mammals
set.seed(1)
cmax_pf=gpf(Data=cdata$data,tree=cdata$phy,
            frmla.phylo=maxCFR_all.viruses~phylo+virusesWithCFR_all.viruses,
            family=gaussian,algorithm='phylo',nfactors=6,min.group.size=10)
HolmProcedure(cmax_pf) #5

#coronaviridae maxCFR  mammals
set.seed(1)
cmax_pf_cov=gpf(Data=cdata_cov$data,tree=cdata_cov$phy,
                frmla.phylo=maxCFR_coronaviridae~phylo+virusesWithCFR_coronaviridae,
                family=gaussian,algorithm='phylo',nfactors=5,min.group.size=10)
HolmProcedure(cmax_pf_cov) #2

#flaviviridae maxCFR mammals
set.seed(1)
cmax_pf_fla=gpf(Data=cdata_fla$data,tree=cdata_fla$phy,
                frmla.phylo=maxCFR_flaviviridae~phylo+virusesWithCFR_flaviviridae,
                family=gaussian,algorithm='phylo',nfactors=6,min.group.size=10)
HolmProcedure(cmax_pf_fla) #5

#rhabdoviridae maxCFR mammals
set.seed(1)
cmax_pf_rha=gpf(Data=cdata_rha$data,tree=cdata_rha$phy,
                frmla.phylo=maxCFR_rhabdoviridae~phylo+virusesWithCFR_rhabdoviridae,
                family=gaussian,algorithm='phylo',nfactors=5,min.group.size=10)
HolmProcedure(cmax_pf_rha) #2

#togaviridae maxCFR mammals
set.seed(1)
cmax_pf_tog=gpf(Data=cdata_tog$data,tree=cdata_tog$phy,
                frmla.phylo=maxCFR_togaviridae~phylo+virusesWithCFR_togaviridae,
                family=gaussian,algorithm='phylo',nfactors=5,min.group.size=10)
HolmProcedure(cmean_pf_tog) #1

#paramyxoviridae maxCFR mammals
set.seed(1)
cmax_pf_par=gpf(Data=cdata_par$data,tree=cdata_par$phy,
                frmla.phylo=maxCFR_paramyxoviridae~phylo+virusesWithCFR_paramyxoviridae,
                family=gaussian,algorithm='phylo',nfactors=5,min.group.size=10)
HolmProcedure(cmean_pf_par) #0 

## summarize all viruses, 5 virus families
#save what each of the clades are for each vfamily
cmax_pf_results=pfsum(cmax_pf)$results #5
cmax_pf_results_cov=pfsum(cmax_pf_cov)$results #2
cmax_pf_results_fla=pfsum(cmax_pf_fla)$results #5
cmax_pf_results_rha=pfsum(cmax_pf_rha)$results #2
cmax_pf_results_tog=pfsum(cmax_pf_tog)$results #1
#cmax_pf_results_par=pfsum(cmean_pf_par)$results #0 factors

#onward transmission all viruses, mammals
set.seed(1)
cot_pf=gpf(Data=cdata2$data,tree=cdata2$phy,
           frmla.phylo=cbind(htrans_all.viruses,ntrans_all.viruses)~phylo,
           family=binomial,algorithm='phylo',nfactors=5,min.group.size=10)
HolmProcedure(cot_pf) #4

#save dataframes for 5 virus families
cdata2_cov_ot<-cdata2[!is.na(cdata2$data$htrans_coronaviridae),]
cdata2_cov_ot<-cdata2_cov_ot[!is.na(cdata2_cov_ot$data$ntrans_coronaviridae),]
cdata2_fla_ot<-cdata2[!is.na(cdata2$data$htrans_flaviviridae),]
cdata2_fla_ot<-cdata2_fla_ot[!is.na(cdata2_fla_ot$data$ntrans_flaviviridae),]
cdata2_rha_ot<-cdata2[!is.na(cdata2$data$htrans_rhabdoviridae),]
cdata2_rha_ot<-cdata2_rha_ot[!is.na(cdata2_rha_ot$data$ntrans_rhabdoviridae),]
cdata2_tog_ot<-cdata2[!is.na(cdata2$data$htrans_togaviridae),]
cdata2_tog_ot<-cdata2_tog_ot[!is.na(cdata2_tog_ot$data$ntrans_togaviridae),]
cdata2_par_ot<-cdata2[!is.na(cdata2$data$htrans_paramyxoviridae),]
cdata2_par_ot<-cdata2_par_ot[!is.na(cdata2_par_ot$data$ntrans_paramyxoviridae),]

#coronaviridae ot mammals
set.seed(1)
cot_pf_cov=gpf(Data=cdata2_cov_ot$data,tree=cdata2_cov_ot$phy,
               frmla.phylo=cbind(htrans_coronaviridae,ntrans_coronaviridae)~phylo,
               family=binomial,algorithm='phylo',nfactors=5,min.group.size=10)
HolmProcedure(cot_pf_cov) #0
#sanity check for 0
cdata2_cov_ot$data$on.frac_coronaviridae 
# no variance, 0 is correct, all cor have onward transmission

#flaviviridae ot mammals
set.seed(1)
cot_pf_fla=gpf(Data=cdata2_fla_ot$data,tree=cdata2_fla_ot$phy,
               frmla.phylo=cbind(htrans_flaviviridae,ntrans_flaviviridae)~phylo,
               family=binomial,algorithm='phylo',nfactors=5,min.group.size=10)
HolmProcedure(cot_pf_fla) #3

#rhabdoviridae ot mammals
set.seed(1)
cot_pf_rha=gpf(Data=cdata2_rha_ot$data,tree=cdata2_rha_ot$phy,
               frmla.phylo=cbind(htrans_rhabdoviridae,ntrans_rhabdoviridae)~phylo,
               family=binomial,algorithm='phylo',nfactors=5,min.group.size=10)
HolmProcedure(cot_pf_rha) #0
#sanity check for 0
cdata2_rha_ot$data$on.frac_rhabdoviridae
# no variance, 0 is correct, rhabdoviridae don't have onward transmission

#togaviridae ot mammals
set.seed(1)
cot_pf_tog=gpf(Data=cdata2_tog_ot$data,tree=cdata2_tog_ot$phy,
               frmla.phylo=cbind(htrans_togaviridae,ntrans_togaviridae)~phylo,
               family=binomial,algorithm='phylo',nfactors=5,min.group.size=10)
HolmProcedure(cot_pf_tog) #1

#paramyxoviridae ot mammals
set.seed(1)
cot_pf_par=gpf(Data=cdata2_par$data,tree=cdata2_par$phy,
               frmla.phylo=cbind(htrans_paramyxoviridae,ntrans_paramyxoviridae)~phylo,
               family=binomial,algorithm='phylo',nfactors=5,min.group.size=10)
HolmProcedure(cot_pf_par) #0

#summarize all viruses 5 virus families
cot_pf_results=pfsum(cot_pf)$results #4
#cot_pf_results_cov=pfsum(cot_pf_cov)$results #0
cot_pf_results_fla=pfsum(cot_pf_fla)$results #3
#cot_pf_results_rha=pfsum(cot_pf_rha)$results #0
cot_pf_results_tog=pfsum(cot_pf_tog)$results #1
#cot_pf_results_par=pfsum(cot_pf_par)$results #0

## all viruses bat CFR mean
set.seed(1)
bmean_pf=gpf(Data=bdata$data,tree=bdata$phy,
             frmla.phylo=meanCFR_all.viruses~phylo+virusesWithCFR_all.viruses,
             family=gaussian,algorithm='phylo',nfactors=2,min.group.size=10)
HolmProcedure(bmean_pf) #2

#coronaviridae meanCFR bat 
bdata_cov<-bdata[!is.na(bdata$data$meanCFR_coronaviridae),]
set.seed(1)
bmean_pf_cov=gpf(Data=bdata_cov$data,tree=bdata_cov$phy,
                 frmla.phylo=meanCFR_coronaviridae~phylo+virusesWithCFR_coronaviridae,
                 family=gaussian,algorithm='phylo',nfactors=5,min.group.size=10)
HolmProcedure(bmean_pf_cov) #1

#flaviviridae meanCFR bat
bdata_fla<-bdata[!is.na(bdata$data$meanCFR_flaviviridae),]
set.seed(1)
bmean_pf_fla=gpf(Data=bdata_fla$data,tree=bdata_fla$phy,
                 frmla.phylo=meanCFR_flaviviridae~phylo+virusesWithCFR_flaviviridae,
                 family=gaussian,algorithm='phylo',nfactors=5,min.group.size=10)
HolmProcedure(bmean_pf_fla) #2

#rhabdoviridae meanCFR bat
bdata_rha<-bdata[!is.na(bdata$data$meanCFR_rhabdoviridae),]
set.seed(1)
bmean_pf_rha=gpf(Data=bdata_rha$data,tree=bdata_rha$phy,
                 frmla.phylo=meanCFR_rhabdoviridae~phylo+virusesWithCFR_rhabdoviridae,
                 family=gaussian,algorithm='phylo',nfactors=5,min.group.size=10)
HolmProcedure(bmean_pf_rha) #2

#togaviridae meanCFR bat
bdata_tog<-bdata[!is.na(bdata$data$meanCFR_togaviridae),]
set.seed(1)
bmean_pf_tog=gpf(Data=bdata_tog$data,tree=bdata_tog$phy,
                 frmla.phylo=meanCFR_togaviridae~phylo+virusesWithCFR_togaviridae,
                 family=gaussian,algorithm='phylo',nfactors=5,min.group.size=10)
HolmProcedure(bmean_pf_tog) #0

#paramyxoviridae meanCFR bat
bdata_par<-bdata[!is.na(bdata$data$meanCFR_paramyxoviridae),]
set.seed(1)
bmean_pf_par=gpf(Data=bdata_par$data,tree=bdata_par$phy,
                 frmla.phylo=meanCFR_paramyxoviridae~phylo+virusesWithCFR_paramyxoviridae,
                 family=gaussian,algorithm='phylo',nfactors=5,min.group.size=10)
HolmProcedure(bmean_pf_par) #0 

## summarize meanCFR bat
bmean_pf_results=pfsum(bmean_pf)$results #2
bmean_pf_results_cov=pfsum(bmean_pf_cov)$results #1
bmean_pf_results_fla=pfsum(bmean_pf_fla)$results #2
bmean_pf_results_rha=pfsum(bmean_pf_rha)$results #2
#bmean_pf_results_tog=pfsum(bmean_pf_tog)$results #0 factors
#bmean_pf_results_par=pfsum(bmean_pf_par)$results #0 factors

## all viruses bat CFR max
set.seed(1)
bmax_pf=gpf(Data=bdata$data,tree=bdata$phy,
            frmla.phylo=maxCFR_all.viruses~phylo+virusesWithCFR_all.viruses,
            family=gaussian,algorithm='phylo',nfactors=3,min.group.size=10)
HolmProcedure(bmax_pf) #2

#coronaviridae maxCFR bat
set.seed(1)
bmax_pf_cov=gpf(Data=bdata_cov$data,tree=bdata_cov$phy,
                frmla.phylo=maxCFR_coronaviridae~phylo+virusesWithCFR_coronaviridae,
                family=gaussian,algorithm='phylo',nfactors=5,min.group.size=10)
HolmProcedure(bmax_pf_cov) #1

#flaviviridae maxCFR bat
set.seed(1)
bmax_pf_fla=gpf(Data=bdata_fla$data,tree=bdata_fla$phy,
                frmla.phylo=maxCFR_flaviviridae~phylo+virusesWithCFR_flaviviridae,
                family=gaussian,algorithm='phylo',nfactors=5,min.group.size=10)
HolmProcedure(bmax_pf_fla) #1

#rhabdoviridae maxCFR bat
set.seed(1)
bmax_pf_rha=gpf(Data=bdata_rha$data,tree=bdata_rha$phy,
                frmla.phylo=maxCFR_rhabdoviridae~phylo+virusesWithCFR_rhabdoviridae,
                family=gaussian,algorithm='phylo',nfactors=5,min.group.size=10)
HolmProcedure(bmax_pf_rha) #1

#togaviridae maxCFR bat
set.seed(1)
bmax_pf_tog=gpf(Data=bdata_tog$data,tree=bdata_tog$phy,
                frmla.phylo=maxCFR_togaviridae~phylo+virusesWithCFR_togaviridae,
                family=gaussian,algorithm='phylo',nfactors=5,min.group.size=10)
HolmProcedure(bmax_pf_tog) #0

#paramyxoviridae maxCFR bat
set.seed(1)
bmax_pf_par=gpf(Data=bdata_par$data,tree=bdata_par$phy,
                frmla.phylo=maxCFR_paramyxoviridae~phylo+virusesWithCFR_paramyxoviridae,
                family=gaussian,algorithm='phylo',nfactors=5,min.group.size=10)
HolmProcedure(bmax_pf_par) #0 

## summarize maxCFR bat
bmax_pf_results=pfsum(bmax_pf)$results #2
bmax_pf_results_cov=pfsum(bmax_pf_cov)$results #1
bmax_pf_results_fla=pfsum(bmax_pf_fla)$results #1
bmax_pf_results_rha=pfsum(bmax_pf_rha)$results #1
#bmax_pf_results_tog=pfsum(bmean_pf_tog)$results #0 factors
#bmax_pf_results_par=pfsum(bmean_pf_par)$results #0 factors

## fraction of viruses with onward transmission bats
set.seed(1)
bot_pf=gpf(Data=bdata2$data,tree=bdata2$phy,
           frmla.phylo=cbind(htrans_all.viruses,ntrans_all.viruses)~phylo,
           family=binomial,algorithm='phylo',nfactors=5,min.group.size=10)
HolmProcedure(bot_pf) #2

#dataframes for 5 virus families
bdata2_cov<-bdata[!is.na(bdata$data$htrans_coronaviridae),]
bdata2_cov<-bdata2_cov[!is.na(bdata2_cov$data$ntrans_coronaviridae),]
bdata2_fla<-bdata[!is.na(bdata$data$htrans_flaviviridae),]
bdata2_fla<-bdata2_fla[!is.na(bdata2_fla$data$ntrans_flaviviridae),]
bdata2_rha<-bdata[!is.na(bdata$data$htrans_rhabdoviridae),]
bdata2_rha<-bdata2_rha[!is.na(bdata2_rha$data$ntrans_rhabdoviridae),]
bdata2_tog<-bdata[!is.na(bdata$data$htrans_togaviridae),]
bdata2_tog<-bdata2_tog[!is.na(bdata2_tog$data$ntrans_togaviridae),]
bdata2_par<-bdata[!is.na(bdata$data$htrans_paramyxoviridae),]
bdata2_par<-bdata2_par[!is.na(bdata2_par$data$ntrans_paramyxoviridae),]

#coronaviridae ot bats
set.seed(1)
bot_pf_cov=gpf(Data=bdata2_cov$data,tree=bdata2_cov$phy,
               frmla.phylo=cbind(htrans_coronaviridae,ntrans_coronaviridae)~phylo,
               family=binomial,algorithm='phylo',nfactors=5,min.group.size=10)
HolmProcedure(bot_pf_cov) #0

#flaviviridae ot bats
set.seed(1)
bot_pf_fla=gpf(Data=bdata2_fla$data,tree=bdata2_fla$phy,
               frmla.phylo=cbind(htrans_flaviviridae,ntrans_flaviviridae)~phylo,
               family=binomial,algorithm='phylo',nfactors=5,min.group.size=10)
HolmProcedure(bot_pf_fla) #0

#rhabdoviridae ot bats
set.seed(1)
bot_pf_rha=gpf(Data=bdata2_rha$data,tree=bdata2_rha$phy,
               frmla.phylo=cbind(htrans_rhabdoviridae,ntrans_rhabdoviridae)~phylo,
               family=binomial,algorithm='phylo',nfactors=5,min.group.size=10)
HolmProcedure(bot_pf_rha) #0

#togaviridae ot bats
set.seed(1)
bot_pf_tog=gpf(Data=bdata2_tog$data,tree=bdata2_tog$phy,
               frmla.phylo=cbind(htrans_togaviridae,ntrans_togaviridae)~phylo,
               family=binomial,algorithm='phylo',nfactors=5,min.group.size=10)
HolmProcedure(bot_pf_tog) #0

#paramyxoviridae ot bats
set.seed(1)
bot_pf_par=gpf(Data=bdata2_par$data,tree=bdata2_par$phy,
               frmla.phylo=cbind(htrans_paramyxoviridae,ntrans_paramyxoviridae)~phylo,
               family=binomial,algorithm='phylo',nfactors=5,min.group.size=10)
HolmProcedure(bot_pf_par) #0

bot_pf_results=pfsum(bot_pf)$results #2
#bot_pf_results_cov=pfsum(bot_pf_cov)$results #0
#bot_pf_results_fla=pfsum(bot_pf_fla)$results #0
#bot_pf_results_rha=pfsum(bot_pf_rha)$results #0
#bot_pf_results_tog=pfsum(bot_pf_tog)$results #0
#bot_pf_results_par=pfsum(bot_pf_par)$results #0

#make a dataframe
#add an ID variable
cmean_pf_results$ID<- "cmean"
cmean_pf_results_cov$ID<- "cmean_cov"
cmean_pf_results_fla$ID<- "cmeans_fla"
cmean_pf_results_rha$ID<-"cmean_rha"
cmean_pf_results_tog$ID<- "cmean_tog"

cmax_pf_results$ID<-"cmax"
cmax_pf_results_cov$ID<-"cmax_cov"
cmax_pf_results_fla$ID<-"cmax_fla"
cmax_pf_results_rha$ID<- "cmax_rha"
cmax_pf_results_tog$ID<- "cmax_tog"

cot_pf_results$ID<-"cot"
cot_pf_results_fla$ID<-"cot_fla"
cot_pf_results_tog$ID<-"cot_tog"

bmean_pf_results$ID<-"bmean"
bmean_pf_results_cov$ID<-"bmean_cov"
bmean_pf_results_fla$ID<-"bmean_fla"
bmean_pf_results_rha$ID<-"bmean_rha"

bmax_pf_results$ID<-"bmax"
bmax_pf_results_cov$ID<-"bmax_cov"
bmax_pf_results_fla$ID<-"bmax_fla"
bmax_pf_results_rha$ID<-"bmax_rha"

bot_pf_results$ID<-"bot"

#bind everything together
results<- do.call("rbind", list(cmean_pf_results,cmean_pf_results_cov,cmean_pf_results_fla,
                                cmean_pf_results_rha,cmean_pf_results_tog,cmax_pf_results, cmax_pf_results_cov,
                                cmax_pf_results_fla,cmax_pf_results_rha,cmax_pf_results_tog,cot_pf_results,cot_pf_results_fla,
                                cot_pf_results_tog,bmean_pf_results,bmean_pf_results_cov,bmean_pf_results_fla,
                                bmean_pf_results_rha,bmax_pf_results,bmax_pf_results_cov, bmax_pf_results_fla,bmax_pf_results_rha,
                                bot_pf_results))
setwd("~/Desktop/PCM Class/phylofatality/")
#write.csv(results,"pf_allclades.csv")

## visualize risky clades
## load in risky species data from 03_data mining
#setwd("~/Desktop/PCM Class/phylofatality")
#species=read.csv("pf_riskyspecies.csv")
#clades=read.csv("pf_allclades.csv")
#rclades=read.csv("pf_riskyclades.csv")

#you get the legend in your console
#mammal meanCFR
pf.tree(cmean_pf,factors=1:4,size=0.1,alphas=rep(0.75,5))
pf.tree(cmean_pf_cov,factors=1:2,size=0.1,alphas=rep(0.75,5))
pf.tree(cmean_pf_fla,factors=c(1:2,4),size=0.1,alphas=rep(0.75,5)) 
pf.tree(cmean_pf_rha,factors=1:2,size=0.1,alphas=rep(0.75,5))
pf.tree(cmean_pf_tog,factors=1,size=0.1,alphas=rep(0.75,5))
#mammals maxCFR
pf.tree(cmax_pf,factors=c(1,3),size=0.1,alphas=rep(0.75,5)) 
pf.tree(cmax_pf_cov,factors=1:2,size=0.1,alphas=rep(0.75,5))
pf.tree(cmax_pf_fla,factors=c(1:4),size=0.1,alphas=rep(0.75,5))
pf.tree(cmax_pf_rha,factors=1:2,size=0.1,alphas=rep(0.75,5))
pf.tree(cmax_pf_tog,factors=1,size=0.1,alphas=rep(0.75,5))
#mammal ot
pf.tree(cot_pf,factors=1:4,size=0.1,alphas=rep(0.75,5)) 
pf.tree(cot_pf_fla,factors=1:2,size=0.1,alphas=rep(0.75,5))
pf.tree(cot_pf_tog,factors=1,size=0.1,alphas=rep(0.75,5))
#bat meanCFR
pf.tree(bmean_pf,factors=1,size=0.1,alphas=rep(0.75,5))
pf.tree(bmean_pf_cov,factors=1,size=0.1,alphas=rep(0.75,5))
pf.tree(bmean_pf_fla,factors=1:2,size=0.1,alphas=rep(0.75,5)) 
#bat maxCFR
pf.tree(bmax_pf_cov,factors=1,size=0.1,alphas=rep(0.75,5))
pf.tree(bmax_pf_fla,factors=1,size=0.1,alphas=rep(0.75,5)) 
#bat ot
#none are risky


## save trees
dtree=treeio::full_join(as.treedata(cdata$phy),cdata$data,by="label")
btree=treeio::full_join(as.treedata(bdata$phy),bdata$data,by="label")

#save  trees for each virus family
dtree_cov=treeio::full_join(as.treedata(cdata_cov$phy),cdata_cov$data,by="label")
dtree_fla=treeio::full_join(as.treedata(cdata_fla$phy),cdata_fla$data,by="label")
dtree_rha=treeio::full_join(as.treedata(cdata_rha$phy),cdata_rha$data,by="label")
dtree_tog=treeio::full_join(as.treedata(cdata_tog$phy),cdata_tog$data,by="label")
dtree_par=treeio::full_join(as.treedata(cdata_par$phy),cdata_par$data,by="label")


## fix palette
AlberPalettes <- c("YlGnBu","Reds","BuPu", "PiYG")
AlberColours <- sapply(AlberPalettes, function(a) RColorBrewer::brewer.pal(5, a)[4])
afun=function(x){
  a=AlberColours[1:x]
  return(a)
}

## make low and high
pcols=afun(2)

## set x max
plus=1
pplus=plus+1


## cfr mean: mammal all viruses
gg=ggtree(dtree,size=0.2,layout="circular",
          aes(colour=meanCFR_all.viruses, group=node))+ #change virus
  #scale_colour_manual(values=c("grey80","black"))+
  scale_color_gradient(low="grey90",high="black")+
  guides(colour=F)

#fix labels for the plot below
#cmean_pf_results$factor[1]="1: subclade~of~italic(Natalidae, Molossidae, Vespertilionidae, Nycteridae, Emballonuridae)" 
#cmean_pf_results$factor[3]="3: subclade~of~italic(Rhinopomatidae, Megadeermatidae, Rhinolophidae, Hipposideridae, Pteropodidae, Noctilionidae, Morpoopidae, Phyllostomidae)"

## add clades
for(i in 1:nrow(cmean_pf_results)){ ##cmean_pf changes
  
  gg=gg+
    geom_hilight(node=cmean_pf_results$node[i], ##cmean_pf changes
                 alpha=0.15,
                 fill=ifelse(cmean_pf_results$clade[i]> ##cmean_pf changes
                               cmean_pf_results$other[i],pcols[2],pcols[1]))+ ##cmean_pf changes
    geom_cladelabel(node=cmean_pf_results$node[i], ##cmean_pf changes
                    label=cmean_pf_results$factor[i], ##cmean_pf changes
                    offset=pplus+1,
                    hjust=0.75,
                    offset.text=pplus*5,
                    parse=T,
                    angle=20)
}
#gg+geom_tippoint(aes(colour=meanCFR),shape=15)

## save
gg_cmean=gg
plot(gg_cmean)


## cfr max: mammals, all viruses
gg=ggtree(dtree,size=0.2,layout="circular",
          aes(colour=maxCFR_all.viruses, group=node))+
  #scale_colour_manual(values=c("grey80","black"))+
  scale_color_gradient(low="grey90",high="black")+
  guides(colour=F)

## add clades
for(i in 1:nrow(cmax_pf_results)){
  
  gg=gg+
    geom_hilight(node=cmax_pf_results$node[i],
                 alpha=0.15,
                 fill=ifelse(cmax_pf_results$clade[i]>
                               cmax_pf_results$other[i],pcols[2],pcols[1]))+
    geom_cladelabel(node=cmax_pf_results$node[i],
                    label=cmax_pf_results$factor[i],
                    offset=pplus,
                    hjust=0.75,
                    offset.text=pplus*5,
                    parse=T,
                    angle=20)
}
#gg+geom_tippoint(aes(colour=meanCFR),shape=15)

## save
gg_cmax=gg
plot(gg_cmax)


## ot, mammals, all viruses
gg=ggtree(dtree,size=0.2,layout="circular",
          aes(colour=on.frac_all.viruses, group=node))+
  #scale_colour_manual(values=c("grey80","black"))+
  scale_color_gradient(low="grey90",high="black")+
  guides(colour=F)

## add clades
for(i in 1:nrow(cot_pf_results)){
  
  gg=gg+
    geom_hilight(node=cot_pf_results$node[i],
                 alpha=0.15,
                 fill=ifelse(cot_pf_results$clade[i]>
                               cot_pf_results$other[i],pcols[2],pcols[1]))+
    geom_cladelabel(node=cot_pf_results$node[i],
                    label=cot_pf_results$factor[i],
                    offset=pplus,
                    hjust=0.75,
                    offset.text=pplus*5,
                    parse=T,
                    angle=20)
}
#gg+geom_tippoint(aes(colour=meanCFR),shape=15)

## save
gg_cot=gg
plot(gg_cot)

####coronaviridae
## cfr mean: mammals, coronaviridae
gg=ggtree(dtree_cov,size=0.2,layout="circular",
          aes(colour=meanCFR_coronaviridae, group=node))+
  #scale_colour_manual(values=c("grey80","black"))+
  scale_color_gradient(low="grey90",high="black")+
  guides(colour=F)

## add clades
for(i in 1:nrow(cmean_pf_results_cov)){ 
  
  gg=gg+
    geom_hilight(node=cmean_pf_results_cov$node[i],
                 alpha=0.15,
                 fill=ifelse(cmean_pf_results_cov$clade[i]> 
                               cmean_pf_results_cov$other[i],pcols[2],pcols[1]))+ 
    geom_cladelabel(node=cmean_pf_results_cov$node[i],
                    label=cmean_pf_results_cov$factor[i], 
                    offset=pplus+1,
                    hjust=0.75,
                    offset.text=pplus*5,
                    parse=T,
                    angle=20)
}
#gg+geom_tippoint(aes(colour=meanCFR),shape=15)

## save
gg_cmean=gg
plot(gg_cmean)


## cfr max: mammals, coronaviridaee
gg=ggtree(dtree_cov,size=0.2,layout="circular",
          aes(colour=maxCFR_coronaviridae, group=node))+
  #scale_colour_manual(values=c("grey80","black"))+
  scale_color_gradient(low="grey90",high="black")+
  guides(colour=F)

## add clades
for(i in 1:nrow(cmax_pf_results_cov)){
  
  gg=gg+
    geom_hilight(node=cmax_pf_results_cov$node[i],
                 alpha=0.15,
                 fill=ifelse(cmax_pf_results_cov$clade[i]>
                               cmax_pf_results_cov$other[i],pcols[2],pcols[1]))+
    geom_cladelabel(node=cmax_pf_results_cov$node[i],
                    label=cmax_pf_results_cov$factor[i],
                    offset=pplus,
                    hjust=0.75,
                    offset.text=pplus*5,
                    parse=T,
                    angle=20)
}
#gg+geom_tippoint(aes(colour=meanCFR),shape=15)

## save
gg_cmax=gg
plot(gg_cmax)

####flaviviridae
## cfr mean: mammals, flaviviridae
gg=ggtree(dtree_fla,size=0.2,layout="circular",
          aes(colour=meanCFR_flaviviridae, group=node))+ 
  #scale_colour_manual(values=c("grey80","black"))+
  scale_color_gradient(low="grey90",high="black")+
  guides(colour=F)

## add clades
for(i in 1:nrow(cmean_pf_results_fla)){ 

  gg=gg+
    geom_hilight(node=cmean_pf_results_fla$node[i], 
                 alpha=0.15,
                 fill=ifelse(cmean_pf_results_fla$clade[i]> 
                               cmean_pf_results_fla$other[i],pcols[2],pcols[1]))+ 
    geom_cladelabel(node=cmean_pf_results_fla$node[i], 
                    label=cmean_pf_results_fla$factor[i], 
                    offset=pplus+1,
                    hjust=0.75,
                    offset.text=pplus*5,
                    parse=T,
                    angle=20)
}
#gg+geom_tippoint(aes(colour=meanCFR),shape=15)

## save
gg_cmean=gg
plot(gg_cmean)


## cfr max: mammals, all viruses
gg=ggtree(dtree_fla,size=0.2,layout="circular",
          aes(colour=maxCFR_flaviviridae, group=node))+
  #scale_colour_manual(values=c("grey80","black"))+
  scale_color_gradient(low="grey90",high="black")+
  guides(colour=F)

## add clades
for(i in 1:nrow(cmax_pf_results_fla)){
  
  gg=gg+
    geom_hilight(node=cmax_pf_results_fla$node[i],
                 alpha=0.15,
                 fill=ifelse(cmax_pf_results_fla$clade[i]>
                               cmax_pf_results_fla$other[i],pcols[2],pcols[1]))+
    geom_cladelabel(node=cmax_pf_results_fla$node[i],
                    label=cmax_pf_results_fla$factor[i],
                    offset=pplus,
                    hjust=0.75,
                    offset.text=pplus*5,
                    parse=T,
                    angle=20)
}
#gg+geom_tippoint(aes(colour=meanCFR),shape=15)

## save
gg_cmax=gg
plot(gg_cmax)

## ot, mammals, flaviviridae
gg=ggtree(dtree_fla,size=0.2,layout="circular",
          aes(colour=on.frac_flaviviridae, group=node))+
  #scale_colour_manual(values=c("grey80","black"))+
  scale_color_gradient(low="grey90",high="black")+
  guides(colour=F)

## add clades
for(i in 1:nrow(cot_pf_results_fla)){
  
  gg=gg+
    geom_hilight(node=cot_pf_results_fla$node[i],
                 alpha=0.15,
                 fill=ifelse(cot_pf_results_fla$clade[i]>
                               cot_pf_results_fla$other[i],pcols[2],pcols[1]))+
    geom_cladelabel(node=cot_pf_results_fla$node[i],
                    label=cot_pf_results_fla$factor[i],
                    offset=pplus,
                    hjust=0.75,
                    offset.text=pplus*5,
                    parse=T,
                    angle=20)
}
#gg+geom_tippoint(aes(colour=meanCFR),shape=15)

## save
gg_cot=gg
plot(gg_cot)


