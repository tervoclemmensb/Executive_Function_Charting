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
coglongdata_baseline<-coglongdata[coglongdata$visitnum==1,c("id","Ageatvisit",allefvars)]
coglongdata_baseline[,allefvars]<-lapply(coglongdata_baseline[,allefvars],function(x){as.numeric(as.character(x))})###all continuous

#########cut based on quartiles
coglongdata_baseline$metaage<-coglongdata_baseline$Ageatvisit
quartilecuts<-quantile(coglongdata_baseline$metaage,probs=c(0,.25,.5,.75,1))
coglongdata_baseline$agebinsfaclabel<-cut(coglongdata_baseline$metaage, breaks =  quartilecuts)
coglongdata_baseline<-coglongdata_baseline %>% group_by(agebinsfaclabel) %>% mutate(agebinmedian=median(unique(round(metaage))))
require(boot)
FAbyagebin<-lapply(unique(coglongdata_baseline$agebinsfaclabel[!is.na(coglongdata_baseline$agebinsfaclabel)]),function(labeli){
  print(labeli)
  bindata<-coglongdata_baseline[coglongdata_baseline$agebinsfaclabel==labeli,]
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
save(FAbyagebindf,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/Factorfigs/FAbyagebin.Luna.Rdata")

