library(dplyr)
library(tidyr)
library(lubridate)
library(ggplot2)
library(lsmeans)
library(mgcv)
library(itsadug)
library(lme4)
library(lsmeans)
library(stats)
library(psych)
library(LNCDR)
library(FactoMineR)
library(corrplot)
library(lavaan)
########load all data#######
load(file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/Factorfigs/allfactordataLuna.Rdata")
Lunascreebaseline<-alleigssaves$screedata
Lunascreebaseline$type<-"baseline"
Lunascreelongitudinal<-alleigssaves$screedatawithin
Lunascreelongitudinal$type<-"within"
Lunaallscree<-plyr::rbind.fill(Lunascreebaseline,Lunascreelongitudinal)
Lunaallscree$dataset<-"Luna"

###Variance
Lunavarbaseline<-data.frame(alleigssaves$varaccount)
Lunavarbaseline<-data.frame(t(Lunavarbaseline[row.names(Lunavarbaseline)=="Proportion Var",]))
Lunavarbaseline$index<-seq(1:nrow(Lunavarbaseline))
Lunavarbaseline$type<-"baseline"

Lunavarlongitudinal<-data.frame(alleigssaves$varaccountwithin)
Lunavarlongitudinal<-data.frame(t(Lunavarlongitudinal[row.names(Lunavarlongitudinal)=="Proportion Var",]))
Lunavarlongitudinal$index<-seq(1:nrow(Lunavarlongitudinal))
Lunavarlongitudinal$type<-"longitudinal"

allunavar<-plyr::rbind.fill(Lunavarbaseline,Lunavarlongitudinal) 
allunavar$dataset<-"Luna"

###loading
Lunaloading<-alleigssaves$loadingbetween
Lunaloading$dataset<-"Luna\nbaseline"

Lunaloadingwithin<-alleigssaves$loadingwithin
Lunaloadingwithin$dataset<-"Luna\nlongitudinal"
######NCANDA#####
load(file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/Factorfigs/allfactordataNCANDA.Rdata")
NCANDAscreebaseline<-alleigssaves$screedata
NCANDAscreebaseline$type<-"baseline"
NCANDAscreelongitudinal<-alleigssaves$screedatawithin
NCANDAscreelongitudinal$type<-"within"
NCANDAallscree<-plyr::rbind.fill(NCANDAscreebaseline,NCANDAscreelongitudinal)
NCANDAallscree$dataset<-"NCANDA"

NCANDAvarbaseline<-data.frame(alleigssaves$varaccount)
NCANDAvarbaseline<-data.frame(t(NCANDAvarbaseline[row.names(NCANDAvarbaseline)=="Proportion Var",]))
NCANDAvarbaseline$index<-seq(1:nrow(NCANDAvarbaseline))
NCANDAvarbaseline$type<-"baseline"

NCANDAvarlongitudinal<-data.frame(alleigssaves$varaccountwithin)
NCANDAvarlongitudinal<-data.frame(t(NCANDAvarlongitudinal[row.names(NCANDAvarlongitudinal)=="Proportion Var",]))
NCANDAvarlongitudinal$index<-seq(1:nrow(NCANDAvarlongitudinal))
NCANDAvarlongitudinal$type<-"longitudinal"

allNCANDAvar<-plyr::rbind.fill(NCANDAvarbaseline,NCANDAvarlongitudinal) 
allNCANDAvar$dataset<-"NCANDA"

###loading
NCANDAloading<-alleigssaves$loadingbetween
NCANDAloading$dataset<-"NCANDA\nbaseline"

NCANDAloadinglongitudinal<-alleigssaves$loadingwithin
NCANDAloadinglongitudinal$dataset<-"NCANDA\nlongitudinal"

####
load(file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/Factorfigs/allfactordataNKI.Rdata")
NKIscree<-alleigssaves$screedata
NKIscree$dataset<-"NKI"

NKIvar<-data.frame(alleigssaves$varaccount)
NKIvar<-data.frame(t(NKIvar[row.names(NKIvar)=="Proportion Var",]))
NKIvar$index<-seq(1:nrow(NKIvar))
NKIvar$dataset<-"NKI"

NKIloading<-alleigssaves$loadingbetween
NKIloading$dataset<-"NKI"

####
load(file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/Factorfigs/allfactordataPNC.Rdata")
PNCscree<-alleigssaves$screedata
PNCscree$dataset<-"PNC"

PNCvar<-data.frame(alleigssaves$varaccount)
PNCvar<-data.frame(t(PNCvar[row.names(PNCvar)=="Proportion Var",]))
PNCvar$index<-seq(1:nrow(PNCvar))
PNCvar$dataset<-"PNC"

PNCloading<-alleigssaves$loadingbetween
PNCloading$dataset<-"PNC"

#################
###all scree####
allscree<-plyr::rbind.fill(Lunaallscree,NCANDAallscree) %>% plyr::rbind.fill(.,NKIscree) %>% plyr::rbind.fill(.,PNCscree)
allscreebar<-allscree %>% group_by(vals) %>% dplyr::summarize(n=n())
allscreebar<-allscreebar[rev(1:nrow(allscreebar)),] %>% dplyr::mutate(cumusumn=cumsum(n)) ##reverse cumsum so can plot % supported by criteria
allscreebar$percentn<-allscreebar$n/sum(allscreebar$n)
allscreebar$percentcumusumn<-allscreebar$cumusumn/sum(allscreebar$n)
ggallscree<-ggplot()+geom_bar(data=allscreebar,aes(x=vals,y=percentcumusumn*100),stat="identity",colour="black",fill="black",width = .5,position=position_nudge(x = 0.25))+scale_x_continuous(limits=c(.9,12),breaks=seq(1:12))+scale_y_continuous(expand=c(0,0))
ggallscree<-LNCDR::lunaize(ggallscree)+xlab("Factor #")+ylab("%")

###all EF variance#####
allvar<-plyr::rbind.fill(allunavar,allNCANDAvar) %>% plyr::rbind.fill(.,NKIvar) %>% plyr::rbind.fill(.,PNCvar)
allvar$label<-allvar$dataset
allvar$label[allvar$dataset %in% c("Luna","NCANDA")]<-paste(allvar$label[allvar$dataset %in% c("Luna","NCANDA")],allvar$type[allvar$dataset %in% c("Luna","NCANDA")])

allvaraverage<-allvar %>% group_by(index) %>% dplyr::summarize(medianvar=median(Proportion.Var))
ggallvar<-ggplot(allvaraverage,aes(x=index,y=medianvar))+geom_line(size=1.15)+geom_jitter(data=allvar,aes(x=index,y=Proportion.Var,shape=label),width=.005)+
  scale_x_continuous(breaks=seq(1:12),limits=c(.9,12.1))+scale_y_continuous(breaks=seq(0,.35,by=.05))
ggallvar<-LNCDR::lunaize(ggallvar)+ylab("Total EF Variance Explained")+xlab("Factor #")+theme(legend.title = element_blank(),legend.position=c(.75,.5),legend.text=element_text(size=15))

##combine plots
ggallscreeforpanel<-ggallscree

require(patchwork)
panelvarandscree<-ggallvar/ggallscreeforpanel+plot_layout(heights=c(1,.25))

####loadingsf or first facto########
####individual plot calls (i.e., not faceted) because facet wrapping will squish expand bars to fill same area...

allloadings<-plyr::rbind.fill(Lunaloading,Lunaloadingwithin) %>% plyr::rbind.fill(.,NCANDAloading) %>% plyr::rbind.fill(.,NCANDAloadinglongitudinal) %>% plyr::rbind.fill(.,NKIloading) %>% plyr::rbind.fill(.,PNCloading)
allloadingsfirstfac<-allloadings[allloadings$factor=="Factor1",]

allloadingplots<-lapply(unique(allloadingsfirstfac$dataset),function(di){
  diloading<-allloadingsfirstfac[allloadingsfirstfac$dataset==di,]
  diloading$typebygroup<-paste(diloading$type,diloading$group,sep="_") ###type first for ordering in plot
  diloading$groupbytypef<-factor(diloading$groupbytype,levels=rev(unique(diloading$groupbytype[order(diloading$typebygroup)])))
  
  if(length(unique(diloading$groupbytypef))<12){
  scalexdiscretelimits=c(levels(factor(diloading$groupbytypef)),sprintf("Null%s",seq(12-length(levels(factor(diloading$groupbytypef))))))###adding nulls to ensure equal bar width across datasets
  }else{
    scalexdiscretelimits<-levels(diloading$groupbytypef)
  }
  
  ggloadingdi<-ggplot(diloading, aes(y=groupbytypef, x=dataset, fill=loadingflip)) + 
    geom_tile()+
    scale_fill_gradient2(name = "loading", 
                         high = "#c9270e", mid = "white", low = "#2f42bde6", 
                         midpoint=0, guide="none",
                         limits=c(-1,1)) +
    xlab("") + #improve y-axis label
    ylab("")+
    #facet_wrap(~dataset,scales="free")+
    #facet_grid(.~dataset,scales="free_y")+
    scale_y_discrete(limits = scalexdiscretelimits)+
    scale_x_discrete(expand = c(0,0))+
    geom_text(aes(label=loadinglabel),size=5)
  ggloadingdi<-LNCDR::lunaize(ggloadingdi)+ylab("")+theme(axis.line.y = element_blank())+theme(axis.text.x =element_text(size=15))+
    theme(axis.text.y =element_text(size=12))+theme(plot.margin = unit(c(0, 0, 0, 0), "null"))
    
  
  return(ggloadingdi)
  
})

allloaings<-allloadingplots[[1]]+allloadingplots[[2]]+allloadingplots[[3]]+allloadingplots[[4]]+allloadingplots[[5]]+allloadingplots[[6]]+plot_layout(ncol = 6)

allfactorfigs<-ggallvar/ggallscreeforpanel/(allloadingplots[[1]]+allloadingplots[[2]]+allloadingplots[[3]]+allloadingplots[[4]]+allloadingplots[[5]]+allloadingplots[[6]]+plot_layout(ncol = 6))+plot_layout(heights=c(1.5,.25,1))
allsaveplot<-list(ggallvar=ggallvar,ggallscreeforpanel=ggallscreeforpanel,allloaings=allloaings)
save(allsaveplot,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/Factorfigs/allfactormainfigs.Rdata")
####commented out, switched to above loop function for plotting simplicity 
# ggloadingluna<-ggplot(allloadingsfirstfac, aes(y=groupbytypef, x=1, fill=loadingflip)) + 
#   geom_tile()+
#   scale_fill_gradient2(name = "loading", 
#                        high = "#c9270e", mid = "white", low = "#2f42bde6", 
#                        midpoint=0, guide="none",
#                        limits=c(-1,1)) +
#   xlab("") + #improve y-axis label
#   ylab("")+
#   #facet_wrap(~dataset,scales="free")+
#   facet_grid(.~dataset,scales="free_y")+
#   #scale_y_discrete(limits = scalexdiscretelimits)+
#   scale_x_discrete(expand = c(0,0))+
#   geom_text(aes(label=loadinglabel),size=5)
# ggloadingluna<-LNCDR::lunaize(ggloadingluna)+ylab("")+theme(axis.line.y = element_blank())+theme(axis.text.x=element_blank())
# 
# ###
# NCANDAloading<-allloadingsfirstfac[allloadingsfirstfac$dataset=="NCANDA",]
# NCANDAloading$groupbytypeff<-factor(NCANDAloading$groupbytypef)
# scalexdiscretelimits=c(levels(factor(NCANDAloading$groupbytypef)),sprintf("Null%s",seq(12-length(levels(factor(NCANDAloading$groupbytypef))))))###adding nulls to ensure equal bar width across datasets
# 
# ggloadingNCANDA<-ggplot(NCANDAloading, aes(y=groupbytypeff, x=1, fill=loadingflip)) + 
#   geom_tile()+
#   scale_fill_gradient2(name = "loading", 
#                        high = "#c9270e", mid = "white", low = "#2f42bde6", 
#                        midpoint=0, guide="none",
#                        limits=c(-1,1)) +
#   xlab("") + #improve y-axis label
#   ylab("")+
#   #facet_wrap(~dataset,scales="free")+
#   scale_y_discrete(limits = scalexdiscretelimits)+
#   scale_x_discrete(expand = c(0,0))+
#   geom_text(aes(label=loadinglabel),size=5)
# ggloadingNCANDA<-LNCDR::lunaize(ggloadingNCANDA)+ylab("")+theme(axis.line.y = element_blank())+theme(axis.text.x=element_blank())
# 
# ####
# 
# 
# ggloadingNKI<-ggplot(NKIloading, aes(y=groupbytypef, x=1, fill=loadingflip)) + 
#   geom_tile()+
#   scale_fill_gradient2(name = "loading", 
#                        high = "#c9270e", mid = "white", low = "#2f42bde6", 
#                        midpoint=0, guide="none",
#                        limits=c(-1,1)) +
#   xlab("") + #improve y-axis label
#   ylab("")+
#   #facet_wrap(~dataset,scales="free")+
#   scale_y_discrete(limits = scalexdiscretelimits)+
#   scale_x_discrete(expand = c(0,0))+
#   geom_text(aes(label=loadinglabel),size=5)
# ggloadingNKI<-LNCDR::lunaize(ggloadingNKI)+ylab("")+theme(axis.line.y = element_blank())+theme(axis.text.x=element_blank())
# 
# ####
# 
# PNCloading<-allloadingsfirstfac[allloadingsfirstfac$dataset=="PNC",]
# PNCloading$typebygroup<-paste(PNCloading$type,PNCloading$group,sep="_") ###type first for ordering in plot
# PNCloading$groupbytypef<-factor(PNCloading$groupbytype,levels=rev(unique(PNCloading$groupbytype[order(PNCloading$typebygroup)])))
# scalexdiscretelimits=c(levels(factor(PNCloading$groupbytypef)),sprintf("Null%s",seq(12-length(levels(factor(PNCloading$groupbytypef))))))###adding nulls to ensure equal bar width across datasets
# 
# ggloadingPNC<-ggplot(PNCloading, aes(y=groupbytypef, x=1, fill=loadingflip)) + 
#   geom_tile()+
#   scale_fill_gradient2(name = "loading", 
#                        high = "#c9270e", mid = "white", low = "#2f42bde6", 
#                        midpoint=0, guide="none",
#                        limits=c(-1,1)) +
#   xlab("") + #improve y-axis label
#   ylab("")+
#   #facet_wrap(~dataset,scales="free")+
#   scale_y_discrete(limits = scalexdiscretelimits)+
#   scale_x_discrete(expand = c(0,0))+
#   geom_text(aes(label=loadinglabel),size=5)
# ggloadingPNC<-LNCDR::lunaize(ggloadingPNC)+ylab("")+theme(axis.line.y = element_blank())+theme(axis.text.x=element_blank())
# 
