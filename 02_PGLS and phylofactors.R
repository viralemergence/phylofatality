## phylofatality
## 02_PGLS and phylofactor
## danbeck@ou.edu
## last update 8/25/23

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

## load data
setwd("~/Desktop/phylofatality")
data=read.csv("CFRbySpecies.csv")

## load Upham phylogeny
setwd("~/Desktop/phylofatality/phylo")
tree=read.nexus('MamPhy_fullPosterior_BDvr_Completed_5911sp_topoCons_NDexp_MCC_v2_target.tre')

## load in taxonomy
taxa=read.csv('taxonomy_mamPhy_5911species.csv',header=T)
taxa$tip=taxa$Species_Name

## fix tip
tree$tip.label=sapply(strsplit(tree$tip.label,'_'),function(x) paste(x[1],x[2],sep=' '))
taxa$species=sapply(strsplit(taxa$tip,'_'),function(x) paste(x[1],x[2],sep=' '))

## species in data
data$species=capitalize(data$Host)

## match
miss=setdiff(data$species,taxa$species)

## flag
data$flag=ifelse(data$species%in%miss,1,0)

## fix data names from CLOVER
setwd("~/Desktop/clover/clover/clover_0.1_mammalviruses/phylogenies")
tdata=read.csv("mammal_phylo_translations.csv")
tdata=tdata[!duplicated(tdata$Host),]

## merge
tdata$X=NULL
tdata$species=capitalize(tdata$Host)
data=merge(data,tdata[c("species","Host_Upham")],by="species",all.x=T)

## fix
data$species2=ifelse(data$flag==1,data$Host_Upham,data$species)

## clean
data$species=data$species2
data$species2=NULL
data$Host_Upham=NULL
data$flag=NULL

## reflag
data$flag=ifelse(is.na(data$species),1,0)
data$species=ifelse(data$flag==1,capitalize(data$Host),data$species)

# ## manual fix
# data$species=revalue(data$species,
#                      c("Allochrocebus preussi"="Cercopithecus preussi",
#                        "Apodemus chejuensis"="Apodemus agrarius",
#                        "Bos taurus x bison bison"="Bos taurus",
#                        "Cavia cutleri"="Cavia tschudii",
#                        "Cercopithecus doggetti"="Cercopithecus mitis",
#                        "Cercopithecus kandti"="Cercopithecus mitis",
#                        "Cercopithecus roloway"="Cercopithecus diana",
#                        ))

## rematch
setdiff(data$species,taxa$species)

## remove missing
data=data[!data$species%in%setdiff(data$species,taxa$species),]

## merge data and taxa
data=merge(data,taxa[c("species","tiplabel","gen","fam","ord","clade")],by="species",all.x=T)
rm(tdata,taxa)

## remove duplicates
data=data[!duplicated(data$species),]

## trim tree
tree=keep.tip(tree,data$species)

## save
data$label=data$species
data$Species=data$species

## merge
cdata=comparative.data(phy=tree,data=data,names.col=species,vcv=T,na.omit=F,warn.dropped=T)

## taxonomy
cdata$data$taxonomy=paste(cdata$data$fam,cdata$data$gen,cdata$data$Species,sep='; ')

## pagel's lambda on mammals
mod_me=pgls(meanCFR~1,data=cdata,lambda="ML")
mod_mx=pgls(maxCFR~1,data=cdata,lambda="ML")

## subset to bats
bdata=cdata[cdata$data$ord=="CHIROPTERA",]

## pagel's lambda
bmod_me=pgls(meanCFR~1,data=bdata,lambda="ML")
bmod_mx=pgls(maxCFR~1,data=bdata,lambda="ML")

## summarize estimates
mlist=list(mod_me,mod_mx,bmod_me,bmod_mx)
pdata=data.frame(dataset=c("all mammals","all mammals","only bats","only bats"),
                 variable=c(rep(c("mean","max"),2)),
                 lambda=sapply(mlist,function(x) x$param["lambda"]))
pdata$variable=factor(pdata$variable,levels=c("mean","max"))

## plot
ggplot(pdata,aes(variable,lambda))+
  theme_bw()+
  geom_segment(aes(x=variable,xend=variable,y=0,yend=lambda))+
  geom_point(size=2)+
  ylim(0,1)+
  facet_wrap(~dataset)+
  labs(y=expression(paste("Pagel's ",lambda)),
       x="case fatality rate of host viruses in humans")+
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
    
    ## get means
    ms=(tapply(dat[,resp],dat[,paste0(resp,'_pf',i)],mean))
    
    ## add in
    results[i,'clade']=ms['factor']
    results[i,'other']=ms['other']
    
  }
  
  ## return
  return(list(set=dat,results=results))
}

## CFR mean
set.seed(1)
cmean_pf=gpf(Data=cdata$data,tree=cdata$phy,
           frmla.phylo=meanCFR~phylo,
           family=gaussian,algorithm='phylo',nfactors=5,min.group.size=5)

## summarize
cmean_pf_results=pfsum(cmean_pf)$results

## CFR max
set.seed(1)
cmax_pf=gpf(Data=cdata$data,tree=cdata$phy,
             frmla.phylo=maxCFR~phylo,
             family=gaussian,algorithm='phylo',nfactors=6,min.group.size=5)

## summarize
cmax_pf_results=pfsum(cmax_pf)$results

## bat CFR mean
set.seed(1)
bmean_pf=gpf(Data=bdata$data,tree=bdata$phy,
             frmla.phylo=meanCFR~phylo,
             family=gaussian,algorithm='phylo',nfactors=5,min.group.size=5)

## summarize
bmean_pf_results=pfsum(bmean_pf)$results

## bat CFR max
set.seed(1)
bmax_pf=gpf(Data=bdata$data,tree=bdata$phy,
            frmla.phylo=maxCFR~phylo,
            family=gaussian,algorithm='phylo',nfactors=6,min.group.size=5)

## summarize
bmax_pf_results=pfsum(bmax_pf)$results

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
gg=ggtree(dtree,size=0.2,aes(colour=meanCFR))+
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
gg=ggtree(dtree,size=0.2,aes(colour=maxCFR))+
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