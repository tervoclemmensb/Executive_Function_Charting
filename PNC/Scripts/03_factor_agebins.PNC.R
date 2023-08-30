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

#########cut based on quartiles
coglongdata$metaage<-coglongdata$age
quartilecuts<-quantile(coglongdata$metaage,probs=c(0,.25,.5,.75,1))
coglongdata$agebinsfaclabel<-cut(coglongdata$metaage, breaks =  quartilecuts)
coglongdata<-coglongdata %>% group_by(agebinsfaclabel) %>% mutate(agebinmedian=median(unique(round(metaage))))
require(boot)
FAbyagebin<-lapply(unique(coglongdata$agebinsfaclabel[!is.na(coglongdata$agebinsfaclabel)]),function(labeli){
  print(labeli)
  bindata<-coglongdata[coglongdata$agebinsfaclabel==labeli,]
  print(nrow(bindata))
  
  bindata[,allefvars]<-lapply(bindata[,allefvars],function(x){as.numeric(as.character(x))})###all continuous
  # parallel<-psych::fa.parallel(bindata[,allefvars],n.iter=1000,fm="ml",quant=.95)
  # nscree<-nFactors::nScree(x=parallel$fa.values,aparallel=parallel$fa.sim,model="factors")
  # screedata<-data.frame(vars=unlist(lapply(strsplit(names(nscree$Components),"n"),"[[",2)),vals=unlist(nscree$Components[1,]))
  # screedata$label<-paste(screedata$vars,screedata$vals,sep=" : ")
  # screedata$finlabel<-paste(c("N Factors",screedata$label),collapse="\n")
  # favals<-fa(bindata[,allefvars],nfactors=3,rotate ="bifactor",fm="ml",scores="tenBerge")
  bootfunc<-function(data,i){
  d2 <- data[i,]
  favalsforeigvis<-fa(d2,nfactors=length(allefvars),rotate ="bifactor",fm="ml",scores="tenBerge")
  favalsforeigvisdf<-data.frame(favalsforeigvis$Vaccounted)
  propvarbin<-favalsforeigvisdf[row.names(favalsforeigvisdf)=="Proportion Var",]
  return(propvarbin$ML1)
  }
  
  bootpropvarbinML1<-boot(data=bindata[,allefvars],bootfunc,R=1000,parallel="multicore",ncpus=4)
  #propvarbin<-data.frame(summary(bootpropvarbinML1))
  propvarbin<-data.frame(ML1propvarboot=bootpropvarbinML1$t)
  propvarbin$agebinsfaclabel<-labeli
  propvarbin$agebin<-unique(bindata$agebinmedian[!is.na(bindata$agebinmedian)])
  ####
  # screedata<-data.frame(vars=unlist(lapply(strsplit(names(nscree$Components),"n"),"[[",2)),vals=unlist(nscree$Components[1,]))
  # screedata$label<-paste(screedata$vars,screedata$vals,sep=" : ")
  # screedata$finlabel<-paste(c("N Factors",screedata$label),collapse="\n")
  # 
  return(propvarbin)
})

FAbyagebindf<-do.call(rbind,FAbyagebin)
ggplot(FAbyagebindf,aes(x=agebin,y=ML1propvarboot,group=agebin))+geom_boxplot()+geom_smooth(method="gam",formula = y~s(x,k=4),se=FALSE,aes(group=1))
save(FAbyagebindf,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/Factorfigs/FAbyagebin.PNC.Rdata")

