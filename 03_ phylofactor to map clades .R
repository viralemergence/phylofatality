## phylofatality
##03_phylofactor to map clades
## danbeck@ou.edu carolinecummings@ou.edu
## last update 4/3/2024

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
setwd("~/Desktop/GitHub/phylofatality/csv files")
data=read.csv("CFRbySpecies.csv")

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

## save
data$label=data$species
data$Species=data$species

## define non-onward viruses (adding 6 columns)
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


###Run Phylofactor
##1 all viruses meanCFR mammals
set.seed(1)
cmean_pf=gpf(Data=cdata$data,tree=cdata$phy,
             frmla.phylo=meanCFR_all.viruses~phylo+virusesWithCFR_all.viruses,
             family=gaussian,algorithm='phylo',nfactors=5,min.group.size=10)
HolmProcedure(cmean_pf) #4

#2 coronaviridae meanCFR mammals
set.seed(1)
cmean_pf_cov=gpf(Data=cdata_cov$data,tree=cdata_cov$phy,
                 frmla.phylo=meanCFR_coronaviridae~phylo+virusesWithCFR_coronaviridae,
                 family=gaussian,algorithm='phylo',nfactors=5,min.group.size=10)
HolmProcedure(cmean_pf_cov) #2

#3 flaviviridae meanCFR mammals
set.seed(1)
cmean_pf_fla=gpf(Data=cdata_fla$data,tree=cdata_fla$phy,
                 frmla.phylo=meanCFR_flaviviridae~phylo+virusesWithCFR_flaviviridae,
                 family=gaussian,algorithm='phylo',nfactors=5,min.group.size=10)
HolmProcedure(cmean_pf_fla) #4

#4 rhabdoviridae meanCFR mammals
set.seed(1)
cmean_pf_rha=gpf(Data=cdata_rha$data,tree=cdata_rha$phy,
                 frmla.phylo=meanCFR_rhabdoviridae~phylo+virusesWithCFR_rhabdoviridae,
                 family=gaussian,algorithm='phylo',nfactors=5,min.group.size=10)
HolmProcedure(cmean_pf_rha) #2

#5 togaviridae meanCFR mammals
set.seed(1)
cmean_pf_tog=gpf(Data=cdata_tog$data,tree=cdata_tog$phy,
                 frmla.phylo=meanCFR_togaviridae~phylo+virusesWithCFR_togaviridae,
                 family=gaussian,algorithm='phylo',nfactors=5,min.group.size=10)
HolmProcedure(cmean_pf_tog) #1

#6 paramyxoviridae meanCFR mammals
set.seed(1)
cmean_pf_par=gpf(Data=cdata_par$data,tree=cdata_par$phy,
                 frmla.phylo=meanCFR_paramyxoviridae~phylo+virusesWithCFR_paramyxoviridae,
                 family=gaussian,algorithm='phylo',nfactors=5,min.group.size=10)
HolmProcedure(cmean_pf_par) #0 

## summarize Mean CFR mammals, all viruses, 5 virus families
#save what each of the clades are for each vfamily
cmean_pf_results=pfsum(cmean_pf)$results
cmean_pf_results_cov=pfsum(cmean_pf_cov)$results #4
cmean_pf_results_fla=pfsum(cmean_pf_fla)$results #2
cmean_pf_results_rha=pfsum(cmean_pf_rha)$results #4
cmean_pf_results_tog=pfsum(cmean_pf_tog)$results #2
#cmean_pf_results_par=pfsum(cmean_pf_par)$results #0 factors

##1 all viruses MaxCFR mammals
set.seed(1)
cmax_pf=gpf(Data=cdata$data,tree=cdata$phy,
            frmla.phylo=maxCFR_all.viruses~phylo+virusesWithCFR_all.viruses,
            family=gaussian,algorithm='phylo',nfactors=6,min.group.size=10)
HolmProcedure(cmax_pf) #5

#2 coronaviridae MaxCFR mammals
set.seed(1)
cmax_pf_cov=gpf(Data=cdata_cov$data,tree=cdata_cov$phy,
                frmla.phylo=maxCFR_coronaviridae~phylo+virusesWithCFR_coronaviridae,
                family=gaussian,algorithm='phylo',nfactors=5,min.group.size=10)
HolmProcedure(cmax_pf_cov) #2

#3 flaviviridae maxCFR mammals
set.seed(1)
cmax_pf_fla=gpf(Data=cdata_fla$data,tree=cdata_fla$phy,
                frmla.phylo=maxCFR_flaviviridae~phylo+virusesWithCFR_flaviviridae,
                family=gaussian,algorithm='phylo',nfactors=6,min.group.size=10)
HolmProcedure(cmax_pf_fla) #5

#4 rhabdoviridae maxCFR mammals
set.seed(1)
cmax_pf_rha=gpf(Data=cdata_rha$data,tree=cdata_rha$phy,
                frmla.phylo=maxCFR_rhabdoviridae~phylo+virusesWithCFR_rhabdoviridae,
                family=gaussian,algorithm='phylo',nfactors=5,min.group.size=10)
HolmProcedure(cmax_pf_rha) #2

#5 togaviridae maxCFR mammals
set.seed(1)
cmax_pf_tog=gpf(Data=cdata_tog$data,tree=cdata_tog$phy,
                frmla.phylo=maxCFR_togaviridae~phylo+virusesWithCFR_togaviridae,
                family=gaussian,algorithm='phylo',nfactors=5,min.group.size=10)
HolmProcedure(cmean_pf_tog) #1

#6 paramyxoviridae maxCFR mammals
set.seed(1)
cmax_pf_par=gpf(Data=cdata_par$data,tree=cdata_par$phy,
                frmla.phylo=maxCFR_paramyxoviridae~phylo+virusesWithCFR_paramyxoviridae,
                family=gaussian,algorithm='phylo',nfactors=5,min.group.size=10)
HolmProcedure(cmean_pf_par) #0 

## summarize Max CFR all viruses, 5 virus families
#save what each of the clades are for each vfamily
cmax_pf_results=pfsum(cmax_pf)$results #5
cmax_pf_results_cov=pfsum(cmax_pf_cov)$results #2
cmax_pf_results_fla=pfsum(cmax_pf_fla)$results #5
cmax_pf_results_rha=pfsum(cmax_pf_rha)$results #2
cmax_pf_results_tog=pfsum(cmax_pf_tog)$results #1
#cmax_pf_results_par=pfsum(cmean_pf_par)$results #0 factors

#1 onward transmission all viruses, mammals
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

#2 coronaviridae ot mammals
set.seed(1)
cot_pf_cov=gpf(Data=cdata2_cov_ot$data,tree=cdata2_cov_ot$phy,
               frmla.phylo=cbind(htrans_coronaviridae,ntrans_coronaviridae)~phylo,
               family=binomial,algorithm='phylo',nfactors=5,min.group.size=10)
HolmProcedure(cot_pf_cov) #0
#sanity check for 0
cdata2_cov_ot$data$on.frac_coronaviridae 
# no variance, 0 is correct, all cor have onward transmission

#3 flaviviridae ot mammals
set.seed(1)
cot_pf_fla=gpf(Data=cdata2_fla_ot$data,tree=cdata2_fla_ot$phy,
               frmla.phylo=cbind(htrans_flaviviridae,ntrans_flaviviridae)~phylo,
               family=binomial,algorithm='phylo',nfactors=5,min.group.size=10)
HolmProcedure(cot_pf_fla) #3

#4 rhabdoviridae ot mammals
set.seed(1)
cot_pf_rha=gpf(Data=cdata2_rha_ot$data,tree=cdata2_rha_ot$phy,
               frmla.phylo=cbind(htrans_rhabdoviridae,ntrans_rhabdoviridae)~phylo,
               family=binomial,algorithm='phylo',nfactors=5,min.group.size=10)
HolmProcedure(cot_pf_rha) #0
#sanity check for 0
cdata2_rha_ot$data$on.frac_rhabdoviridae
# no variance, 0 is correct, rhabdoviridae don't have onward transmission

#5 togaviridae ot mammals
set.seed(1)
cot_pf_tog=gpf(Data=cdata2_tog_ot$data,tree=cdata2_tog_ot$phy,
               frmla.phylo=cbind(htrans_togaviridae,ntrans_togaviridae)~phylo,
               family=binomial,algorithm='phylo',nfactors=5,min.group.size=10)
HolmProcedure(cot_pf_tog) #1

#6 paramyxoviridae ot mammals
set.seed(1)
cot_pf_par=gpf(Data=cdata2_par$data,tree=cdata2_par$phy,
               frmla.phylo=cbind(htrans_paramyxoviridae,ntrans_paramyxoviridae)~phylo,
               family=binomial,algorithm='phylo',nfactors=5,min.group.size=10)
HolmProcedure(cot_pf_par) #0

#summarize OT all viruses 5 virus families
cot_pf_results=pfsum(cot_pf)$results #4
#cot_pf_results_cov=pfsum(cot_pf_cov)$results #0
cot_pf_results_fla=pfsum(cot_pf_fla)$results #3
#cot_pf_results_rha=pfsum(cot_pf_rha)$results #0
cot_pf_results_tog=pfsum(cot_pf_tog)$results #1
#cot_pf_results_par=pfsum(cot_pf_par)$results #0

##1 all viruses bat CFR mean
set.seed(1)
bmean_pf=gpf(Data=bdata$data,tree=bdata$phy,
             frmla.phylo=meanCFR_all.viruses~phylo+virusesWithCFR_all.viruses,
             family=gaussian,algorithm='phylo',nfactors=2,min.group.size=10)
HolmProcedure(bmean_pf) #2

#2 coronaviridae meanCFR bat 
bdata_cov<-bdata[!is.na(bdata$data$meanCFR_coronaviridae),]
set.seed(1)
bmean_pf_cov=gpf(Data=bdata_cov$data,tree=bdata_cov$phy,
                 frmla.phylo=meanCFR_coronaviridae~phylo+virusesWithCFR_coronaviridae,
                 family=gaussian,algorithm='phylo',nfactors=5,min.group.size=10)
HolmProcedure(bmean_pf_cov) #1

#3 flaviviridae meanCFR bat
bdata_fla<-bdata[!is.na(bdata$data$meanCFR_flaviviridae),]
set.seed(1)
bmean_pf_fla=gpf(Data=bdata_fla$data,tree=bdata_fla$phy,
                 frmla.phylo=meanCFR_flaviviridae~phylo+virusesWithCFR_flaviviridae,
                 family=gaussian,algorithm='phylo',nfactors=5,min.group.size=10)
HolmProcedure(bmean_pf_fla) #2

#4 rhabdoviridae meanCFR bat
bdata_rha<-bdata[!is.na(bdata$data$meanCFR_rhabdoviridae),]
set.seed(1)
bmean_pf_rha=gpf(Data=bdata_rha$data,tree=bdata_rha$phy,
                 frmla.phylo=meanCFR_rhabdoviridae~phylo+virusesWithCFR_rhabdoviridae,
                 family=gaussian,algorithm='phylo',nfactors=5,min.group.size=10)
HolmProcedure(bmean_pf_rha) #2

#5 togaviridae meanCFR bat
bdata_tog<-bdata[!is.na(bdata$data$meanCFR_togaviridae),]
set.seed(1)
bmean_pf_tog=gpf(Data=bdata_tog$data,tree=bdata_tog$phy,
                 frmla.phylo=meanCFR_togaviridae~phylo+virusesWithCFR_togaviridae,
                 family=gaussian,algorithm='phylo',nfactors=5,min.group.size=10)
HolmProcedure(bmean_pf_tog) #0

#6 paramyxoviridae meanCFR bat
bdata_par<-bdata[!is.na(bdata$data$meanCFR_paramyxoviridae),]
set.seed(1)
bmean_pf_par=gpf(Data=bdata_par$data,tree=bdata_par$phy,
                 frmla.phylo=meanCFR_paramyxoviridae~phylo+virusesWithCFR_paramyxoviridae,
                 family=gaussian,algorithm='phylo',nfactors=5,min.group.size=10)
HolmProcedure(bmean_pf_par) #0 

## summarize MeanCFR bat
bmean_pf_results=pfsum(bmean_pf)$results #2
bmean_pf_results_cov=pfsum(bmean_pf_cov)$results #1
bmean_pf_results_fla=pfsum(bmean_pf_fla)$results #2
bmean_pf_results_rha=pfsum(bmean_pf_rha)$results #2
#bmean_pf_results_tog=pfsum(bmean_pf_tog)$results #0 factors
#bmean_pf_results_par=pfsum(bmean_pf_par)$results #0 factors

##1 all viruses bat CFR max
set.seed(1)
bmax_pf=gpf(Data=bdata$data,tree=bdata$phy,
            frmla.phylo=maxCFR_all.viruses~phylo+virusesWithCFR_all.viruses,
            family=gaussian,algorithm='phylo',nfactors=3,min.group.size=10)
HolmProcedure(bmax_pf) #2

#2 coronaviridae maxCFR bat
set.seed(1)
bmax_pf_cov=gpf(Data=bdata_cov$data,tree=bdata_cov$phy,
                frmla.phylo=maxCFR_coronaviridae~phylo+virusesWithCFR_coronaviridae,
                family=gaussian,algorithm='phylo',nfactors=5,min.group.size=10)
HolmProcedure(bmax_pf_cov) #1

#3 flaviviridae maxCFR bat
set.seed(1)
bmax_pf_fla=gpf(Data=bdata_fla$data,tree=bdata_fla$phy,
                frmla.phylo=maxCFR_flaviviridae~phylo+virusesWithCFR_flaviviridae,
                family=gaussian,algorithm='phylo',nfactors=5,min.group.size=10)
HolmProcedure(bmax_pf_fla) #1

#4 rhabdoviridae maxCFR bat
set.seed(1)
bmax_pf_rha=gpf(Data=bdata_rha$data,tree=bdata_rha$phy,
                frmla.phylo=maxCFR_rhabdoviridae~phylo+virusesWithCFR_rhabdoviridae,
                family=gaussian,algorithm='phylo',nfactors=5,min.group.size=10)
HolmProcedure(bmax_pf_rha) #1

#5 togaviridae maxCFR bat
set.seed(1)
bmax_pf_tog=gpf(Data=bdata_tog$data,tree=bdata_tog$phy,
                frmla.phylo=maxCFR_togaviridae~phylo+virusesWithCFR_togaviridae,
                family=gaussian,algorithm='phylo',nfactors=5,min.group.size=10)
HolmProcedure(bmax_pf_tog) #0

#6 paramyxoviridae maxCFR bat
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

##1 fraction of viruses with onward transmission bats
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

#2 coronaviridae ot bats
set.seed(1)
bot_pf_cov=gpf(Data=bdata2_cov$data,tree=bdata2_cov$phy,
               frmla.phylo=cbind(htrans_coronaviridae,ntrans_coronaviridae)~phylo,
               family=binomial,algorithm='phylo',nfactors=5,min.group.size=10)
HolmProcedure(bot_pf_cov) #0

#3 flaviviridae ot bats
set.seed(1)
bot_pf_fla=gpf(Data=bdata2_fla$data,tree=bdata2_fla$phy,
               frmla.phylo=cbind(htrans_flaviviridae,ntrans_flaviviridae)~phylo,
               family=binomial,algorithm='phylo',nfactors=5,min.group.size=10)
HolmProcedure(bot_pf_fla) #0

#4 rhabdoviridae ot bats
set.seed(1)
bot_pf_rha=gpf(Data=bdata2_rha$data,tree=bdata2_rha$phy,
               frmla.phylo=cbind(htrans_rhabdoviridae,ntrans_rhabdoviridae)~phylo,
               family=binomial,algorithm='phylo',nfactors=5,min.group.size=10)
HolmProcedure(bot_pf_rha) #0

#5 togaviridae ot bats
set.seed(1)
bot_pf_tog=gpf(Data=bdata2_tog$data,tree=bdata2_tog$phy,
               frmla.phylo=cbind(htrans_togaviridae,ntrans_togaviridae)~phylo,
               family=binomial,algorithm='phylo',nfactors=5,min.group.size=10)
HolmProcedure(bot_pf_tog) #0

#6 paramyxoviridae ot bats
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
setwd("~/Desktop/GitHub/phylofatality/csv files")
write.csv(results,"pf_allclades.csv")

#================================================================================
##quick intermission --> go to 03a_data mining and run this script --> come back

## save trees
dtree=treeio::full_join(as.treedata(cdata$phy),cdata$data,by="label")
btree=treeio::full_join(as.treedata(bdata$phy),bdata$data,by="label")

#save  trees for each virus family
dtree_cov=treeio::full_join(as.treedata(cdata_cov$phy),cdata_cov$data,by="label")
dtree_fla=treeio::full_join(as.treedata(cdata_fla$phy),cdata_fla$data,by="label")
dtree_rha=treeio::full_join(as.treedata(cdata_rha$phy),cdata_rha$data,by="label")
dtree_tog=treeio::full_join(as.treedata(cdata_tog$phy),cdata_tog$data,by="label")
dtree_par=treeio::full_join(as.treedata(cdata_par$phy),cdata_par$data,by="label")

#===============================================================================
####Plotting
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


##1 CFR mean: mammal_all viruses
gg=ggtree(dtree,size=0.2, layout="circular",
          aes(colour=meanCFR_all.viruses, group=node))+ 
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
                    offset.text=pplus*10,
                    parse=T,
                    angle=20)
}
#plot
gg_cmean<- gg+ 
  ggtitle("MeanCFR-All Viruses")+ 
  theme(plot.title = element_text(hjust = 0.5, size=8))
plot(gg_cmean)

#save
setwd("~/Desktop/GitHub/phylofatality/figs")
#ggsave("MeanCFR_allviruses.jpg", gg_cmean, device = "jpeg", width = 10, height = 6, units = "in")

###2 CFR max: mammals_all viruses
gg=ggtree(dtree,size=0.2,layout="circular",
          aes(colour=maxCFR_all.viruses, group=node))+
  #scale_colour_manual(values=c("grey80","black"))+
  scale_color_gradient(low="grey90",high="black")+
  guides(colour=F)

## add clades
for(i in 1:nrow(cmax_pf_results)){ #cmax_pf_results
  
  gg=gg+
    geom_hilight(node=cmax_pf_results$node[i],
                 alpha=0.15,
                 fill=ifelse(cmax_pf_results$clade[i]>
                               cmax_pf_results$other[i],pcols[2],pcols[1]))+
    geom_cladelabel(node=cmax_pf_results$node[i],
                    label=cmax_pf_results$factor[i],
                    offset=pplus,
                    hjust=0.75,
                    offset.text=pplus*10,
                    parse=T,
                    angle=20)
}

#plot
gg_cmax<- gg+
  geom_tippoint(aes(colour=maxCFR_all.viruses),shape=15)+
  ggtitle("MaxCFR-All Viruses")+ 
  theme(plot.title = element_text(hjust = 0.5, size=8))

## save
plot(gg_cmax)
#ggsave("MaxCFR_allviruses.jpg", gg_cmax, device = "jpeg", width = 10, height = 6, units = "in")


##3 ot, mammals_all viruses
gg=ggtree(dtree,size=0.2,layout="circular",
          aes(colour=on.frac_all.viruses, group=node))+
  #scale_colour_manual(values=c("grey80","black"))+
  scale_color_gradient(low="grey90",high="black")+
  guides(colour=F)

## add clades
for(i in 1:nrow(cot_pf_results)){ #cot_pf_results
  
  gg=gg+
    geom_hilight(node=cot_pf_results$node[i],
                 alpha=0.15,
                 fill=ifelse(cot_pf_results$clade[i]>
                               cot_pf_results$other[i],pcols[2],pcols[1]))+
    geom_cladelabel(node=cot_pf_results$node[i],
                    label=cot_pf_results$factor[i],
                    offset=pplus,
                    hjust=0.75,
                    offset.text=pplus*10,
                    parse=T,
                    angle=20)
}

#plot
gg_cot=gg+
  ggtitle("Fraction with Onward Transmission-All Viruses")+
  geom_tippoint(aes(colour=on.frac_all.viruses),shape=15)+
  theme(plot.title = element_text(hjust = 0.5, size=8))

## save
plot(gg_cot)
#ggsave("OT_allviruses.jpg", gg_cot, device = "jpeg", width = 10, height = 6, units = "in")


####4 coronaviridae
##4 CFR mean/max: mammals, coronaviridae
gg=ggtree(dtree_cov,size=0.2,layout="circular",
          aes(colour=meanCFR_coronaviridae, group=node))+
  #scale_colour_manual(values=c("grey80","black"))+
  scale_color_gradient(low="grey90",high="black")+
  guides(colour=F)

## add clades
for(i in 1:nrow(cmean_pf_results_cov)){ #cmean_pf_results_cov
  
  gg=gg+
    geom_hilight(node=cmean_pf_results_cov$node[i],
                 alpha=0.15,
                 fill=ifelse(cmean_pf_results_cov$clade[i]> 
                               cmean_pf_results_cov$other[i],pcols[2],pcols[1]))+ 
    geom_cladelabel(node=cmean_pf_results_cov$node[i],
                    label=cmean_pf_results_cov$factor[i], 
                    offset=pplus+1,
                    hjust=0.75,
                    offset.text=pplus*10,
                    parse=T,
                    angle=20,
                    label.size=12)
}
#plot
gg_cmean_cov <- gg+
  ggtitle(expression("MeanCFR/MaxCFR-" ~ italic("Coronaviridae")))+
  geom_tippoint(aes(colour=meanCFR_coronaviridae),shape=15)+
  theme(plot.title = element_text(hjust = 0.5, size=8))

## save
plot(gg_cmean_cov)
#ggsave("MeanCFR/MaxCFR_coronaviridae.jpg", gg_cmean, device = "jpeg", width = 10, height = 6, units = "in")


####5 flaviviridae
##5 CFR mean: mammals, flaviviridae
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
                    offset.text=pplus*10,
                    parse=T,
                    angle=20)
}
#plot
gg_cmean_fla <- gg+ 
  ggtitle(expression("MeanCFR-" ~ italic("Flaviviridae")))+
  geom_tippoint(aes(colour=meanCFR_flaviviridae),shape=15)+
  theme(plot.title = element_text(hjust = 0.5, size=8))

## save
plot(gg_cmean_fla)
#ggsave("MeanCFR_flaviviridae.jpg", gg_cmean, device = "jpeg", width = 10, height = 6, units = "in")


##6 CFR max: mammals, flaviviridae
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
                    offset.text=pplus*10,
                    parse=T,
                    angle=20)
}
#gg
gg_cmax_fla <- gg+
  ggtitle(expression("MaxCFR-" ~ italic("Flaviviridae")))+
  geom_tippoint(aes(colour=maxCFR_flaviviridae),shape=15)+
  theme(plot.title = element_text(hjust = 0.5, size=8))

## save
plot(gg_cmax_fla)
#ggsave("MaxCFR_flaviviridae.jpg", gg_cmax, device = "jpeg", width = 10, height = 6, units = "in")

##7 OT, mammals, flaviviridae
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
                    offset.text=pplus*10,
                    parse=T,
                    angle=20)
}
#plot
gg_ot_fla <- gg+
  ggtitle(expression("% with Onward Transmission-"~italic("Flaviviridae")))+
  geom_tippoint(aes(colour=on.frac_flaviviridae),shape=15)+
  theme(plot.title = element_text(hjust = 0.5, size=8))

## save
plot(gg_ot_fla)
#ggsave("OT_flaviviridae.jpg", gg_ot_fla, device = "jpeg", width = 10, height = 6, units = "in")

#try to  wrap my plots
giant_phylo<- gg_cmean+gg_cmax+gg_cot+gg_cmean_cov+gg_cmean_fla+gg_cmax_fla+gg_ot_fla

print(giant_phylo)
#ggsave("03_giant_phylofactor.jpg", giant_phylo, device = "jpeg", width = 8, height = 6, units = "in")

