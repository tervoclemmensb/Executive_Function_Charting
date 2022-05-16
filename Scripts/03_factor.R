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
########
############
coglongdata<-read.csv("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/data/btc_R03cleaneddata_20220306.outlierremoved.compositeacclat.csv")
coglongdata$id<-as.factor(coglongdata$id)
####vars############
####################
allefvars<-c("Anti_CRLat","Anti_CRR","Mix_CRLat","Mix_CRR","SOC.Problems.solved.in.minimum.moves",
             "SOC.Overallmeaninitialthinkingtime","DMS.PC","DMS.Median.correct.latency","SSP.Span.length",
             "nfixbreak","best_acc_m_exclude","first_lat_m_exclude")

accvars<-c("Anti_CRR","Mix_CRR","DMS.PC","SSP.Span.length","nfixbreak_fl","SOC.Problems.solved.in.minimum.moves","best_acc_m_exclude_fl")

latvars<-c("Anti_CRLat","Mix_CRLat","SOC.Overallmeaninitialthinkingtime","DMS.Median.correct.latency","first_lat_m_exclude")

allefvars<-c(accvars,latvars)

coglongdata$visitnum<-as.numeric(coglongdata$visit)

groups<-data.frame(acc=c("Accuracycomposite","Anti_CRR","Mix_CRR","DMS.PC","SSP.Span.length","nfixbreak_fl","SOC.Problems.solved.in.minimum.moves","best_acc_m_exclude_fl"),
                   lat=c("Latencycomposite","Anti_CRLat","Mix_CRLat","DMS.Median.correct.latency",NA,NA,"SOC.Overallmeaninitialthinkingtime","first_lat_m_exclude"),
                   group=c("Composite","ANTI","MIX","DMS","SSP","FIX","SOC","MGS"))
groupslong<-gather(groups, type, outcome, acc:lat, factor_key=TRUE)
groupslong$groupbytype<-paste(groupslong$group,groupslong$type,sep="_")
######################
####EFA###############
###baselin#####
coglongdata_baseline<-coglongdata[coglongdata$visitnum==1,c("id",allefvars)]
coglongdata_baseline[,allefvars]<-lapply(coglongdata_baseline[,allefvars],function(x){as.numeric(as.character(x))})###all continuous
#parallel<-psych::fa.parallel(coglongdata_baseline[,allefvars],n.iter=1000,fm="ml",quant=.95)
#save(parallel,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Data/parallelsavefactor.Rdata") 
load(file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Data/parallelsavefactor.Rdata") ##already run load in to save computation time
nscree<-nFactors::nScree(x=parallel$fa.values,aparallel=parallel$fa.sim,model="factors")
screedata<-data.frame(vars=unlist(lapply(strsplit(names(nscree$Components),"n"),"[[",2)),vals=unlist(nscree$Components[1,]))
screedata$label<-paste(screedata$vars,screedata$vals,sep=" : ")
screedata$finlabel<-paste(c("N Factors","Baseline",screedata$label),collapse="\n")
####3 factors by parallel analysis ###1 factors by eigen value of 1 rule#
favals<-fa(coglongdata_baseline[,allefvars],nfactors=3,rotate ="bifactor",fm="ml",scores="tenBerge")
favalsforeigvis<-fa(coglongdata_baseline[,allefvars],nfactors=length(allefvars),rotate ="bifactor",fm="ml",scores="tenBerge")
####
eigvals<-data.frame(eigvals=parallel$fa.values,resamp=parallel$fa.sim)
eigvals$index<-seq(1:nrow(eigvals))
eigvals$type<-"full"
#####################
# #####eig plots######
# eigvals$var<-paste0(as.character(round(favalsforeigvis$Vaccounted[row.names(favalsforeigvis$Vaccounted)=="Proportion Var",],2)*100),"%")
# eigvals$var[eigvals$var %in% "0%"]<-"<1%"
# eigvals$finallabel<-paste(eigvals$index,eigvals$var,sep="\n")
# eigvals$finallabelf<-factor(eigvals$finallabel,levels=eigvals$finallabel)
# ggeigs<-ggplot(eigvals,aes(x=finallabelf,y=eigvals,group=1))+geom_line(size=1)+geom_point(size=3)+geom_line(data=eigvals,aes(x=finallabelf,y=resamp,group=1),linetype="dashed",colour="black",size=1.15)
# ggeigs<-lunaize(ggeigs)+theme(legend.position = "top",legend.title = element_blank())+xlab("Factor #")+ylab("Eigenvalue")+geom_hline(yintercept = mean(eigvals$eigvals),linetype="solid")+
#   #scale_x_continuous(breaks = seq(from = 1, to = 12, by = 1))+#,limits=c(1,12))+
#   geom_text(data=screedata,aes(x=7,y=2.5,label=finlabel),lineheight = .9,size = 5)
# ggsave(ggeigs,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/Factorfigs/Lunaeigsfig.pdf",width=6,height=5)

####loading plots####
loadingplot<-data.frame(unclass(favals$loadings),outcome=row.names(unclass(favals$loadings)))
names(loadingplot)[grep("ML",names(loadingplot))]<-c("Factor1","Factor2","Factor3")
loadingplot_long<-gather(loadingplot, factor, loading, Factor1:Factor3,factor_key=TRUE)               
loadingplot_long$factorlong<-as.character(loadingplot_long$factor)
loadingplot_long$factorlong[loadingplot_long$factor=="Factor1"]<-"Fact. 1\n(DG EF)"
loadingplot_long$factorlong[loadingplot_long$factor=="Factor2"]<-"Fact. 2\n"
loadingplot_long$factorlong[loadingplot_long$factor=="Factor3"]<-"Fact. 3\n"

loadingplot_long<-merge(loadingplot_long,groupslong,by=c("outcome"))
loadingplot_long$typebygroup<-paste(loadingplot_long$type,loadingplot_long$group,sep="_") ###type first for ordering in plot
loadingplot_long$groupbytypef<-factor(loadingplot_long$groupbytype,levels=rev(unique(loadingplot_long$groupbytype[order(loadingplot_long$typebygroup)])))
loadingplot_long$loadingflip<-loadingplot_long$loading

ggloading<-ggplot(loadingplot_long, aes(groupbytypef, loadingflip, fill=loadingflip)) + 
  facet_grid(cols=vars(factorlong),scales="free") + #place the factors in separate facets
  geom_bar(stat="identity") + #make the bars
  coord_flip(expand = FALSE) + #flip the axes so the test names can be horizontal  
  #define the fill color gradient: blue=positive, red=negative
  scale_fill_gradient2(name = "loading", 
                       high = "#c9270e", mid = "white", low = "#2f42bde6", 
                       midpoint=0, guide="none") +
  ylab("Loading") + #improve y-axis label
  xlab("")#+
  #scale_x_discrete(limits = scalexdiscretelimits) 
ggloading<-LNCDR::lunaize(ggloading)+ylab("Factor Loading")+theme(strip.background = element_blank())+
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

####geom_tile version####
no_leadzero <- function(x, integer_zero = TRUE, digits = 2) {
  x <- round(x, digits = digits)
  y <- sprintf(paste0('%.', digits,'f'),x)
  y[x > 0 & x < 1] <- sprintf('.%s',lapply(strsplit(as.character(x[x > 0 & x < 1]),"[.]"),"[[",2))
  y[x > -1 & x < 0] <- sprintf('-.%s',lapply(strsplit(as.character(x[x > -1 & x < 0]),"[.]"),"[[",2))
  if(integer_zero){
    y[(x%%1) == 0] <- sprintf("%s", x[(x%%1) == 0])  
  } else {
    y[x == 0] <- '0'  
  }
  return(y)
}
loadingplot_long$loadinglabel<-no_leadzero(round(loadingplot_long$loadingflip,2))
ggloadingtile<-ggplot(loadingplot_long, aes(y=groupbytypef, x=factorlong, fill=loadingflip)) + 
  geom_tile()+
  scale_fill_gradient2(name = "loading", 
                       high = "#c9270e", mid = "white", low = "#2f42bde6", 
                       midpoint=0,
                       limits = c(-1, 1)) +
  xlab("") + #improve y-axis label
  ylab("")+
  #scale_y_discrete(limits = scalexdiscretelimits)+
  scale_x_discrete(expand = c(0,0))#+
  #geom_text(aes(label=loadinglabel),size=5)
tilelegend<-cowplot::get_legend(ggloadingtile)
ggloadingtile<-LNCDR::lunaize(ggloadingtile)+ylab("")+theme(axis.line.y = element_blank())+theme(legend.position = "none")+theme(axis.text.x=element_blank())

# fullpanel<-ggeigs+ggloadingtile
# ggsave(fullpanel,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/Factorfigs/Lunaeigwithloading.panel.pdf",width=12,height=5)
# ggsave(tilelegend,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/Factorfigs/tileloading.legend.pdf",width=12,height=5)
loadingplot_longbetween<-loadingplot_long
########################
######within############
coglongdata$idnum<-coglongdata$id ###make id numeric for multilevel group in psych statsby
cogdataforwithinfactor<-coglongdata[c("idnum",allefvars)] ##reduce to only efvars and idnumber (grouping var)

stabsycog<-statsBy(cogdataforwithinfactor,group="idnum",cors=TRUE)
withincorrelations<-stabsycog$rwg
ndf<-data.frame(stabsycog$n)
ndf$max<-apply(ndf[, allefvars], 1, max)
ndf$haslongitudinal<-dplyr::if_else(ndf$max>1,1,0)###if n has more than one observation -- longitudinal data
parallelwithin<-psych::fa.parallel(withincorrelations,n.iter=1000,fm="ml",n.obs=length(which(ndf$haslongitudinal==1)),quant=.95)
save(parallelwithin,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Data/parallelwithinsavefactor.Rdata")
load(file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Data/parallelwithinsavefactor.Rdata") ##already run load in to save computation time
nscreewithin<-nFactors::nScree(x=parallelwithin$fa.values,aparallel=parallelwithin$fa.sim,model="factors")
screedatawithin<-data.frame(vars=unlist(lapply(strsplit(names(nscreewithin$Components),"n"),"[[",2)),vals=unlist(nscreewithin$Components[1,]))
screedatawithin$label<-paste(screedatawithin$vars,screedatawithin$vals,sep=" : ")
screedatawithin$finlabel<-paste(c("N Factors","Longitudinal",screedatawithin$label),collapse="\n")

favalswithin<-fa(withincorrelations,nfactors=3,n.obs=length(which(ndf$haslongitudinal==1)),rotate ="bifactor",fm="ml",scores="tenBerge")
favalswithinforeigvis<-fa(withincorrelations,nfactors=length(allefvars),n.obs=length(which(ndf$haslongitudinal==1)),rotate ="bifactor",fm="ml",scores="tenBerge")
####
eigvalswithin<-data.frame(eigvals=parallelwithin$fa.values,resamp=parallelwithin$fa.sim)
eigvalswithin$index<-seq(1:nrow(eigvalswithin))
eigvalswithin$type<-"full"
#####################
#####eig plots######
# eigvalswithin$var<-paste0(as.character(round(favalswithinforeigvis$Vaccounted[row.names(favalswithinforeigvis$Vaccounted)=="Proportion Var",],2)*100),"%")
# eigvalswithin$var[eigvalswithin$var %in% "0%"]<-"<1%"
# eigvalswithin$finallabel<-paste(eigvalswithin$index,eigvalswithin$var,sep="\n")
# eigvalswithin$finallabelf<-factor(eigvalswithin$finallabel,levels=eigvalswithin$finallabel)
# ggeigswithin<-ggplot(eigvalswithin,aes(x=finallabelf,y=eigvals,group=1))+geom_line(size=1)+geom_point(size=3)+geom_line(data=eigvalswithin,aes(x=finallabelf,y=resamp,group=1),linetype="dashed",colour="black",size=1.15)
# ggeigswithin<-lunaize(ggeigswithin)+theme(legend.position = "top",legend.title = element_blank())+xlab("Factor #")+ylab("Eigenvalue")+geom_hline(yintercept = mean(eigvalswithin$eigvals),linetype="solid")+
#   #scale_x_continuous(breaks = seq(from = 1, to = 12, by = 1))+#,limits=c(1,12))+
#   geom_text(data=screedatawithin,aes(x=7,y=1.5,label=finlabel),lineheight = .9,size = 5)
# ggsave(ggeigswithin,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/Factorfigs/Lunaeigsfig.within.pdf",width=6,height=5)

###eig plot together#####
names(eigvalswithin)[!(names(eigvalswithin) %in% c("index"))]<-paste0("within.",names(eigvalswithin)[!(names(eigvalswithin) %in% c("index"))])
alleigs<-merge(eigvals,eigvalswithin,by=c("index"))
alleigs$fulllabel<-paste(alleigs$finallabelf,alleigs$within.var,sep="\n")
alleigs$fulllabelf<-factor(alleigs$fulllabel,levels=alleigs$fulllabel)

ggbaselinelongitudinal<-ggplot(alleigs,aes(x=index,y=eigvals,group=1))+geom_line(size=1)+geom_point(size=3)+geom_line(data=alleigs,aes(x=index,y=resamp,group=1),linetype="dashed",colour="black",size=1.15)+
  geom_hline(yintercept = mean(alleigs$eigvals),linetype="solid")+
  geom_line(aes(x=index,y=within.eigvals),size=1,colour="grey55")+
  geom_point(aes(x=index,y=within.eigvals),size=3,shape=21,fill="white")+
  geom_hline(yintercept = mean(alleigs$within.eigvals ),linetype="solid",colour="grey55")+
  geom_line(data=alleigs,aes(x=index,y=within.resamp,group=1),linetype="dashed",colour="grey55",size=1.15)+
  geom_text(data=screedata,aes(x=5.5,y=1.5,label=finlabel),lineheight = .9,size = 6)+
  geom_text(data=screedatawithin,aes(x=10,y=1.5,label=finlabel),lineheight = .9,size = 6)+
  scale_x_continuous(breaks=seq(1,max(alleigs$index)))

ggbaselinelongitudinal<-lunaize(ggbaselinelongitudinal)+theme(legend.position = "top",legend.title = element_blank())+xlab("Factor #")+ylab("Eigenvalue")

fullpanelwithbotheigs<-ggbaselinelongitudinal+ggloadingtile+plot_layout(widths = c(1,.25))
ggsave(fullpanelwithbotheigs,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/Factorfigs/Lunaeigwithloading.botheigs.panel.pdf",width=8,height=5)

####cors for save###
#base#
cors<-cor(coglongdata_baseline[,allefvars],use="pairwise.complete.obs")
source("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/GeneralScripts/corsave.R")
corssaveclean<-corsaveclean(cors,groupslong)
write.table(lowermerge,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/Factorfigs/Lunabaselinecorrs.csv",row.names=FALSE,col.names=FALSE,sep=",",quote=FALSE)

corslong<-cors %>% rstatix::cor_gather()
corslong<-corslong[corslong$var1!=corslong$var2,] ### remove identity 
corslong$varorder<-unlist(lapply(1:nrow(corslong),function(ri){
  vars<-c(corslong$var1[ri],corslong$var2[ri])
  varorder<-paste(vars[order(vars)],collapse="_")
}))
corslong<-corslong[!duplicated(corslong$varorder),] ###remove symmetric cors, awithb, bwitha (upper/lower tri but preserves labels)
##within##
wgcorsave<-data.frame(stabsycog$rwg)
row.names(wgcorsave)<-gsub(".wg","",row.names(wgcorsave))
names(wgcorsave)<-gsub(".wg","",names(wgcorsave))
corssavecleanwithin<-corsaveclean(wgcorsave,groupslong)
write.table(corssavecleanwithin,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/Factorfigs/Lunalongitudinalcorrs.csv",row.names=FALSE,col.names=FALSE,sep=",",quote=FALSE)

corslongwithin<-stabsycog$rwg %>% rstatix::cor_gather()
corslongwithin<-corslongwithin[corslongwithin$var1!=corslongwithin$var2,] ### remove identity 
corslongwithin$varorder<-unlist(lapply(1:nrow(corslongwithin),function(ri){
  vars<-c(corslongwithin$var1[ri],corslongwithin$var2[ri])
  varorder<-paste(vars[order(vars)],collapse="_")
}))
corslongwithin<-corslongwithin[!duplicated(corslongwithin$varorder),] ###remove symmetric cors, awithb, bwitha (upper/lower tri but preserves labels)

####loading plots####
###within#####
loadingplot<-data.frame(unclass(favalswithin$loadings),outcome=row.names(unclass(favalswithin$loadings)))
names(loadingplot)[grep("ML",names(loadingplot))]<-c("Factor1","Factor2","Factor3")
loadingplot_long<-gather(loadingplot, factor, loading, Factor1:Factor3,factor_key=TRUE)               
loadingplot_long$factorlong<-as.character(loadingplot_long$factor)
loadingplot_long$factorlong[loadingplot_long$factor=="Factor1"]<-"Factor 1\nDomain General EF"
loadingplot_long$factorlong[loadingplot_long$factor=="Factor2"]<-"Factor 2\n"
loadingplot_long$factorlong[loadingplot_long$factor=="Factor3"]<-"Factor 3\n"
loadingplot_long$ogoutcome<-loadingplot_long$outcome
loadingplot_long$outcome<-gsub(".wg","",loadingplot_long$outcome)

loadingplot_long<-merge(loadingplot_long,groupslong,by=c("outcome"))
loadingplot_long$typebygroup<-paste(loadingplot_long$type,loadingplot_long$group,sep="_") ###type first for ordering in plot
loadingplot_long$groupbytypef<-factor(loadingplot_long$groupbytype,levels=rev(unique(loadingplot_long$groupbytype[order(loadingplot_long$typebygroup)])))
loadingplot_long$loadingflip<-loadingplot_long$loading

####geom_tile version####

loadingplot_long$loadinglabel<-no_leadzero(round(loadingplot_long$loadingflip,2))
ggloadingtile<-ggplot(loadingplot_long, aes(y=groupbytypef, x=factorlong, fill=loadingflip)) + 
  geom_tile()+
  scale_fill_gradient2(name = "loading", 
                       high = "#c9270e", mid = "white", low = "#2f42bde6", 
                       midpoint=0, guide="none") +
  xlab("") + #improve y-axis label
  ylab("")+
  #scale_y_discrete(limits = scalexdiscretelimits)+
  scale_x_discrete(expand = c(0,0))+
  geom_text(aes(label=loadinglabel),size=5)
ggloadingtile<-LNCDR::lunaize(ggloadingtile)+ylab("")+theme(axis.line.y = element_blank())

ggsave(ggloadingtile,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/Factorfigs/Lunawithinloadings.pdf",width=6,height=5)


alleigssaves<-list(screedata=screedata,screedatawithin=screedatawithin,
                   varaccount=favalsforeigvis$Vaccounted,
                   varaccountwithin=favalswithinforeigvis$Vaccounted,
                   alleigs=alleigs,
                   corsbaseline=corslong,
                   corswithin=corslongwithin,
                   loadingbetween=loadingplot_longbetween,
                   loadingwithin=loadingplot_long,
                   nbetween=length(unique(coglongdata$id)),
                   nwithin=length(which(ndf$haslongitudinal==1)))
save(alleigssaves,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/Factorfigs/allfactordataLuna.Rdata")
