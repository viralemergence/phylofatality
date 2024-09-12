## phylofatality
## 03_pgls and phylofactors
## danbeck@ou.edu, carolinecummings2018@gmail.com
## last update 09/12/204

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
library(ade4)
library(phytools)

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

## define non-onward viruses
data$ntrans_all.viruses=data$virusesWithOT_all.viruses-data$htrans_all.viruses

## merge
cdata=comparative.data(phy=tree,data=data,names.col=species,vcv=T,na.omit=F,warn.dropped=T)

## taxonomy
cdata$data$taxonomy=paste(cdata$data$fam,cdata$data$gen,cdata$data$Species,sep='; ')

## separate dataset for onward transmission (remove NA)
cdata2=cdata[!is.na(cdata$data$on.frac_all.viruses),]

## pagel's lambda on mammals
mod_me=pgls(meanCFR_all.viruses~1,data=cdata,lambda="ML")
mod_mx=pgls(maxCFR_all.viruses~1,data=cdata,lambda="ML")
mod_ot=pgls(on.frac_all.viruses~1,data=cdata2,lambda="ML")

## we can also implement these tests using phylosig
#psl_me=phylosig(cdata$phy,cdata$data$meanCFR_all.viruses,method="lambda",test=T)
#psl_mx=phylosig(cdata$phy,cdata$data$maxCFR_all.viruses,method="lambda",test=T)
#psl_ot=phylosig(cdata2$phy,cdata2$data$on.frac_all.viruses,method="lambda",test=T)

## compare that values are equivalent
#round(mod_me$param["lambda"],3)==round(psl_me$lambda,3)
#round(mod_mx$param["lambda"],3)==round(psl_mx$lambda,3)
#round(mod_ot$param["lambda"],3)==round(psl_ot$lambda,3)

## we can also use phylosig to estimate Bloomberg's K
#psk_me=phylosig(cdata$phy,cdata$data$meanCFR_all.viruses,method="K",test=T)
#psk_mx=phylosig(cdata$phy,cdata$data$maxCFR_all.viruses,method="K",test=T)
#psk_ot=phylosig(cdata2$phy,cdata2$data$on.frac_all.viruses,method="K",test=T)

## Moran's I will use inverse distances
#d=1/cophenetic(cdata$phy)
#diag(d)=0

## reduced dataset
#d2=1/cophenetic(cdata2$phy)
#diag(d2)=0

## Moran's I
#Moran.I(cdata$data$meanCFR_all.viruses,d)
#Moran.I(cdata$data$maxCFR_all.viruses,d)
#Moran.I(cdata2$data$on.frac_all.viruses,d2)

## correlograms
#form=meanCFR_all.viruses+maxCFR_all.viruses~ord/fam/gen
#form2=on.frac_all.viruses~ord/fam/gen

## run
#plot(correlogram.formula(form,data=cdata$data))
#plot(correlogram.formula(form2,data=cdata2$data))

## similar in ade4
#gearymoran(d,data.frame(cdata$data$meanCFR_all.viruses,
#                       cdata$data$maxCFR_all.viruses),alter="two-sided")
#gearymoran(d2,cdata2$data$on.frac_all.viruses,alter="two-sided")

## subset to bats
bdata=cdata[cdata$data$ord=="CHIROPTERA",]
bdata2=cdata2[cdata2$data$ord=="CHIROPTERA",]

## pagel's lambda
bmod_me=pgls(meanCFR_all.viruses~1,data=bdata,lambda="ML")
bmod_mx=pgls(maxCFR_all.viruses~1,data=bdata,lambda="ML")
bmod_ot=pgls(on.frac_all.viruses~1,data=bdata2,lambda="ML")

## summarize estimates
mlist=list(mod_me,mod_mx,mod_ot,bmod_me,bmod_mx,bmod_ot)
pdata=data.frame(dataset=c(rep("all mammals",3),rep("bats only",3)),
                 variable=c(rep(c("meanCFR","maxCFR","on.frac"),2)),
                 lambda=sapply(mlist,function(x) x$param["lambda"]),
                 lambda_lower=sapply(mlist,function(x) x$param.CI$lambda$ci.val[1]),
                 lambda_upper=sapply(mlist,function(x) x$param.CI$lambda$ci.val[2]),
                 lambda_lower_p=sapply(mlist,function(x) x$param.CI$lambda$bounds.p[1]),
                 lambda_upper_p=sapply(mlist,function(x) x$param.CI$lambda$bounds.val[1]))
pdata$variable=factor(pdata$variable,levels=c("meanCFR","maxCFR","on.frac"))

#save to open later 
setwd("~/Desktop/GitHub/phylofatality/csv files")
write.csv(pdata, "PS_allviruses.csv")

## reopen csv
setwd("~/Desktop/GitHub/phylofatality/csv files")
pdata=read.csv("PS_allviruses.csv")

##plot
##PS of viral virulence measures across all mammals and within bats (all viruses)
ggplot(pdata,aes(variable,lambda))+
  theme_bw()+
  geom_segment(aes(x=variable,xend=variable,y=lambda_lower,yend=lambda_upper))+
  geom_point(size=2)+
  ylim(0,1)+
  facet_wrap(~dataset)+
  labs(y=expression(paste("Pagel's ",lambda)),
       x="viral virulence measures in human hosts")+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(axis.title.x=element_text(margin=margin(t=10,r=0,b=0,l=0)))+
  theme(axis.title.y=element_text(margin=margin(t=0,r=10,b=0,l=0)))

## set taxonomy
taxonomy=data.frame(cdata$data$taxonomy)
names(taxonomy)="taxonomy"
taxonomy$Species=rownames(cdata$data)
taxonomy=taxonomy[c("Species","taxonomy")]
taxonomy$taxonomy=as.character(taxonomy$taxonomy)

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
    em$cat=plyr::revalue(em$phylo,
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

## CFR mean
set.seed(1)
cmean_pf=gpf(Data=cdata$data,tree=cdata$phy,
           frmla.phylo=meanCFR_all.viruses~phylo+virusesWithCFR_all.viruses,
           family=gaussian,algorithm='phylo',nfactors=5,min.group.size=10)
HolmProcedure(cmean_pf) #4

## summarize
cmean_pf_results=pfsum(cmean_pf)$results

## CFR max
set.seed(1)
cmax_pf=gpf(Data=cdata$data,tree=cdata$phy,
             frmla.phylo=maxCFR_all.viruses~phylo+virusesWithCFR_all.viruses,
             family=gaussian,algorithm='phylo',nfactors=6,min.group.size=10)
HolmProcedure(cmax_pf) #5

## summarize
cmax_pf_results=pfsum(cmax_pf)$results

## fraction of viruses with onward transmission
set.seed(1)
cot_pf=gpf(Data=cdata2$data,tree=cdata2$phy,
            frmla.phylo=cbind(htrans_all.viruses,ntrans_all.viruses)~phylo,
            family=binomial,algorithm='phylo',nfactors=5,min.group.size=10)
HolmProcedure(cot_pf) #4

## summarize
cot_pf_results=pfsum(cot_pf)$results

## bat CFR mean
set.seed(1)
bmean_pf=gpf(Data=bdata$data,tree=bdata$phy,
             frmla.phylo=meanCFR_all.viruses~phylo+virusesWithCFR_all.viruses,
             family=gaussian,algorithm='phylo',nfactors=5,min.group.size=10)
HolmProcedure(bmean_pf) #2

## summarize
bmean_pf_results=pfsum(bmean_pf)$results

## bat CFR max
set.seed(1)
bmax_pf=gpf(Data=bdata$data,tree=bdata$phy,
            frmla.phylo=maxCFR_all.viruses~phylo+virusesWithCFR_all.viruses,
            family=gaussian,algorithm='phylo',nfactors=3,min.group.size=10)
HolmProcedure(bmax_pf) #2

## summarize
bmax_pf_results=pfsum(bmax_pf)$results

## fraction of viruses with onward transmission
set.seed(1)
bot_pf=gpf(Data=bdata$data,tree=bdata$phy,
           frmla.phylo=cbind(htrans_all.viruses,ntrans_all.viruses)~phylo,
           family=binomial,algorithm='phylo',nfactors=5,min.group.size=10)
HolmProcedure(bot_pf) #2

## summarize
bot_pf_results=pfsum(bot_pf)$results

## save trees
dtree=treeio::full_join(as.treedata(cdata$phy),cdata$data,by="label")
btree=treeio::full_join(as.treedata(bdata$phy),bdata$data,by="label")

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

## mammal
gg=ggtree(dtree,size=0.2,aes(colour=meanCFR_all.viruses))+
  #scale_colour_manual(values=c("grey80","black"))+
  scale_color_gradient(low="grey90",high="black")+
  guides(colour=F)

## add clades
for(i in 1:nrow(cmean_pf_results)){
  
  gg=gg+
    geom_hilight(node=cmean_pf_results$node[i],
                 alpha=0.15,
                 fill=ifelse(cmean_pf_results$clade>
                               cmean_pf_results$other,pcols[2],pcols[1])[i])+
    geom_cladelabel(node=cmean_pf_results$node[i],
                    label=cmean_pf_results$factor[i],
                    offset=pplus,
                    hjust=0.75,
                    offset.text=pplus*2,
                    parse=T,
                    angle=90)
}
#gg+geom_tippoint(aes(colour=meanCFR),shape=15)

## save
gg_cmean=gg

## cfr max
gg=ggtree(dtree,size=0.2,aes(colour=maxCFR_all.viruses))+
  #scale_colour_manual(values=c("grey80","black"))+
  scale_color_gradient(low="grey90",high="black")+
  guides(colour=F)

## add clades
for(i in 1:nrow(cmax_pf_results)){
  
  gg=gg+
    geom_hilight(node=cmax_pf_results$node[i],
                 alpha=0.15,
                 fill=ifelse(cmax_pf_results$clade>
                               cmax_pf_results$other,pcols[2],pcols[1])[i])+
    geom_cladelabel(node=cmax_pf_results$node[i],
                    label=cmax_pf_results$factor[i],
                    offset=pplus,
                    hjust=0.75,
                    offset.text=pplus*2,
                    parse=T,
                    angle=90)
}
#gg+geom_tippoint(aes(colour=meanCFR),shape=15)

## save
gg_cmax=gg
