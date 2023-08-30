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
coglongdata<-read.csv("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/NCANDA/Data/btc_NCANDAscoredmeasures_20220306.outlierremoved.compositeacclat.csv")
coglongdata$id<-as.factor(coglongdata$subject)
coglongdata$visitnum<-as.numeric(as.factor(coglongdata$visit))
####vars############
####################
allaccvars<-c("cnp_cpf_ifac_tot","cnp_cpw_iwrd_tot","cnp_spcptnl_scpt_tp","cnp_sfnb2_sfnb_mcr","cnp_pmat24a_pmat24_a_cr","cnp_cpfd_dfac_tot","cnp_cpwd_dwrd_tot",
              "cnp_shortvolt_svt","cnp_er40d_er40_cr","cnp_pcet_pcet_acc2","cnp_medf36_medf36_a","cnp_pvoc_pvoccr","cnp_pvrt_pvrt_pc","cnp_svdelay_svt_ld","np_wais4_rawscore_computed","latentdd","Accuracycomposite")   

alllatvars<-c("cnp_cpf_ifac_rtc","cnp_cpw_iwrd_rtc","cnp_spcptnl_scpt_tprt","cnp_sfnb2_sfnb_mrtc","cnp_pmat24a_pmat24_a_rtcr","cnp_cpfd_dfac_rtc",
              "cnp_cpwd_dwrd_rtc","cnp_shortvolt_svtcrt","cnp_er40d_er40_crt","cnp_pcet_pcetrtcr","cnp_medf36_medf36_t","cnp_pvoc_pvocrtcr","cnp_pvrt_pvrtrtcr","cnp_svdelay_svtldrtc","stroop_total_mean","latentgroove","Latencycomposite")

accvars<-c("cnp_pcet_pcet_acc2","cnp_sfnb2_sfnb_mcr","cnp_spcptnl_scpt_tp")
latvars<-c("cnp_pcet_pcetrtcr","cnp_sfnb2_sfnb_mrtc","cnp_spcptnl_scpt_tprt","stroop_total_mean")

allefvars<-c(accvars,latvars)

groups<-data.frame(acc=c("Accuracycomposite","cnp_pcet_pcet_acc2","cnp_spcptnl_scpt_tp","cnp_sfnb2_sfnb_mcr",NA),
                   lat=c("Latencycomposite","cnp_pcet_pcetrtcr","cnp_spcptnl_scpt_tprt","cnp_sfnb2_sfnb_mrtc","stroop_total_mean"),
                   group=c("Composite","PCET","PCTP","PNBK","STRP"))
groupslong<-gather(groups, type, outcome, acc:lat, factor_key=TRUE)
groupslong$groupbytype<-paste(groupslong$group,groupslong$type,sep="_")
######################
####EFA###############
###baselin#####
coglongdata_baseline<-coglongdata[coglongdata$visitnum==1,c("id","cnp_age",allefvars)]
coglongdata_baseline[,allefvars]<-lapply(coglongdata_baseline[,allefvars],function(x){as.numeric(as.character(x))})###all continuous
#########cut based on quartiles
coglongdata_baseline$metaage<-coglongdata_baseline$cnp_age
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
save(FAbyagebindf,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/Factorfigs/FAbyagebin.NCANDA.Rdata")

