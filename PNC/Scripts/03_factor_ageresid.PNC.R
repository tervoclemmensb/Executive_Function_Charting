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
source("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/GeneralScripts/growthrate_mgcvgam.R")
############
coglongdata<-read.csv("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/PNC/Data/btc_PNCscoredmeasures_20220306.outlierremoved.compositeacclat.csv")
####vars############
####################
allcogtestsacc<-c("PADT_A","PFMT_IFAC_TOT","PEIT_CR","PWMT_KIWRD_TOT","PVRT_CR","PEDT_A","PMAT_CR","VOLT_SVT","LNB_MCR","PCET_ACC2","PCPT_T_TP","PLOT_TC","WRAT_CR_RAW")
allcogtestslatency<-c("PADT_T","PFMT_IFAC_RTC","PEIT_CRT","PWMT_KIWRD_RTC","PVRT_RTCR","PEDT_T","MP_MP2RTCR","PMAT_RTCR","VOLT_SVTCRT","LNB_MRTC","PCET_RTCR","PCPT_T_TPRT","PLOT_TCRT")

allcogvars<-c(allcogtestsacc,allcogtestslatency)
allfactorvars<-allcogvars

testgroups<-c("PADT","PFMT","PEIT","PWMT","PVRT","PEDT","PMAT","VOLT","LNB","PCET","PCPT","PLOT","WRAT","MP")
efficiencytestgroups<-c("PADT","PFMT","PEIT","PWMT","PVRT","PEDT","PMAT","VOLT","LNB","PCPT","PLOT")
eftestgroups<-c("PCET","PCPT","LNB") ###EXECUTIVE-CONTROL#HTTPS://WWW.NCBI.NLM.NIH.GOV/PMC/ARTICLES/PMC3295891/PDF/NIHMS348849.PDF

accvars<-unlist(lapply(eftestgroups,function(vg){grep(vg,allcogtestsacc,value=TRUE)}))
latvars<-unlist(lapply(eftestgroups,function(vg){grep(vg,allcogtestslatency,value=TRUE)}))

allefvars<-c(accvars,latvars) 
groups<-data.frame(acc=c("Accuracycomposite","PCET_ACC2","PCPT_T_TP","LNB_MCR"),
                   lat=c("Latencycomposite","PCET_RTCR","PCPT_T_TPRT","LNB_MRTC"),
                   group=c("Composite","PCET","PCTP","LNB"))
groupslong<-gather(groups, type, outcome, acc:lat, factor_key=TRUE)
groupslong$groupbytype<-paste(groupslong$group,groupslong$type,sep="_")
######################
####EFA###############
coglongdata[,allefvars]<-lapply(coglongdata[,allefvars],function(x){as.numeric(as.character(x))})###all continuous
coglongdata_baselineageresid<-coglongdata[c("dbGaP_Subject_ID","age",allefvars)]
coglongdata_baselineageresid[,allefvars]<-lapply(coglongdata_baselineageresid[,allefvars],function(x){mgcvresidualizecovar(x,coglongdata_baselineageresid$age,coglongdata_baselineageresid)})
parallel<-psych::fa.parallel(coglongdata_baselineageresid[,allefvars],n.iter=1000,fm="ml",quant=.95)
save(parallel,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/PNC/Data/parallelsavefactor.ageresid.Rdata") 
load(file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/PNC/Data/parallelsavefactor.ageresid.Rdata") ##already run load in to save computation time
nscree<-nFactors::nScree(x=parallel$fa.values,aparallel=parallel$fa.sim,model="factors")
screedata<-data.frame(vars=unlist(lapply(strsplit(names(nscree$Components),"n"),"[[",2)),vals=unlist(nscree$Components[1,]))
screedata$label<-paste(screedata$vars,screedata$vals,sep=" : ")
screedata$finlabel<-paste(c("N Factors",screedata$label),collapse="\n")
####3-4 factors by parallel analysis ###1 factors by eigen value of 1 rule#
favals<-fa(coglongdata_baselineageresid[,allefvars],nfactors=3,rotate ="bifactor",fm="ml",scores="tenBerge")
favalsforeigvis<-fa(coglongdata_baselineageresid[,allefvars],nfactors=max(c(max(screedata$vals),3)),rotate ="bifactor",fm="ml",scores="tenBerge")
####
eigvals<-data.frame(eigvals=parallel$fa.values,resamp=parallel$fa.sim)
eigvals$index<-seq(1:nrow(eigvals))
eigvals$type<-"full"
#####################
#####eig plots######
eigvals$var<-paste0(as.character(round(favalsforeigvis$Vaccounted[row.names(favalsforeigvis$Vaccounted)=="Proportion Var",],2)*100),"%")
eigvals$var[eigvals$var %in% "0%"]<-"<1%"
eigvals$finallabel<-paste(eigvals$index,eigvals$var,sep="\n")
eigvals$finallabelf<-factor(eigvals$finallabel,levels=eigvals$finallabel)
ggeigs<-ggplot(eigvals,aes(x=index,y=eigvals,group=1))+geom_line(size=1)+geom_point(size=3)+geom_line(data=eigvals,aes(x=index,y=resamp,group=1),linetype="dashed",colour="black",size=1.15)
ggeigs<-lunaize(ggeigs)+theme(legend.position = "top",legend.title = element_blank())+xlab("Factor #")+ylab("Eigenvalue")+geom_hline(yintercept = mean(eigvals$eigvals),linetype="solid")+
  #scale_x_continuous(breaks = seq(from = 1, to = 12, by = 1))+#,limits=c(1,12))+
  geom_text(data=screedata,aes(x=4,y=1.15,label=finlabel),lineheight = .9,size = 6)+
  scale_x_continuous(breaks=seq(1,max(eigvals$index)))
ggsave(ggeigs,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/Factorfigs/PNCeigsfig.ageresid.pdf",width=6,height=5)

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
#loadingplot_long$loadingflip[loadingplot_long$factor=="Factor1"]<-loadingplot_long$loading[loadingplot_long$factor=="Factor1"]*-1 ###flip sign for consistency across datasets (postive for accuracy, negative for latency)

scalexdiscretelimits=c(levels(loadingplot_long$groupbytypef),sprintf("Null%s",seq(12-length(levels(loadingplot_long$groupbytypef)))))###adding nulls to ensure equal bar width across datasets

ggloading<-ggplot(loadingplot_long, aes(groupbytypef, loadingflip, fill=loadingflip)) + 
  facet_grid(cols=vars(factorlong),scales="free") + #place the factors in separate facets
  geom_bar(stat="identity") + #make the bars
  coord_flip(expand = FALSE) + #flip the axes so the test names can be horizontal  
  #define the fill color gradient: blue=positive, red=negative
  scale_fill_gradient2(name = "loading", 
                       high = "#c9270e", mid = "white", low = "#2f42bde6", 
                       midpoint=0, guide="none") +
  ylab("Loading") + #improve y-axis label
  xlab("")+
  scale_x_discrete(limits = scalexdiscretelimits) 
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
                       midpoint=0, guide="none",
                       limits=c(-1,1)) +
  xlab("") + #improve y-axis label
  ylab("")+
  scale_y_discrete(limits = scalexdiscretelimits)+
  scale_x_discrete(expand = c(0,0))#+
  #geom_text(aes(label=loadinglabel),size=5)
ggloadingtile<-LNCDR::lunaize(ggloadingtile)+ylab("")+theme(axis.line.y = element_blank())+theme(axis.text.x=element_blank())

fullpanel<-ggeigs+ggloadingtile+plot_layout(widths = c(1,.25))
ggsave(fullpanel,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/Factorfigs/PNCeigwithloading.panel.ageresid.pdf",width=8,height=5)

cors<-cor(coglongdata[,allefvars],use="pairwise.complete.obs")
corslong<-cors %>% rstatix::cor_gather()
corslong<-corslong[corslong$var1!=corslong$var2,] ### remove identity 
corslong$varorder<-unlist(lapply(1:nrow(corslong),function(ri){
  vars<-c(corslong$var1[ri],corslong$var2[ri])
  varorder<-paste(vars[order(vars)],collapse="_")
}))
corslong<-corslong[!duplicated(corslong$varorder),] ###remove symetric cors, awithb, bwitha (upper/lower tri but preserves labels)

alleigssaves<-list(screedata=screedata,
                   varaccount=favalsforeigvis$Vaccounted,
                   eigs=eigvals,
                   corslong=corslong,
                   loadingbetween=loadingplot_long)
save(alleigssaves,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/Factorfigs/allfactordataPNC.ageresid.Rdata")

