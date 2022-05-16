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
library(mgcv)
############
source("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/GeneralScripts/metabyage.R")
########
coglongdata<-read.csv("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/NCANDA/Data/btc_NCANDAscoredmeasures_20220214.outlierremoved.csv")
####vars############
####################
###Exclude crystalized intelligence measures from composites 

allaccvars<-c("cnp_cpf_ifac_tot","cnp_cpw_iwrd_tot","cnp_spcptnl_scpt_tp","cnp_sfnb2_sfnb_mcr","cnp_pmat24a_pmat24_a_cr","cnp_cpfd_dfac_tot","cnp_cpwd_dwrd_tot",
              "cnp_shortvolt_svt","cnp_er40d_er40_cr","cnp_pcet_pcet_acc2","cnp_medf36_medf36_a","cnp_pvoc_pvoccr","cnp_pvrt_pvrt_pc","cnp_svdelay_svt_ld","np_wais4_rawscore_computed","latentdd","Accuracycomposite")   

alllatvars<-c("cnp_cpf_ifac_rtc","cnp_cpw_iwrd_rtc","cnp_spcptnl_scpt_tprt","cnp_sfnb2_sfnb_mrtc","cnp_pmat24a_pmat24_a_rtcr","cnp_cpfd_dfac_rtc",
              "cnp_cpwd_dwrd_rtc","cnp_shortvolt_svtcrt","cnp_er40d_er40_crt","cnp_pcet_pcetrtcr","cnp_medf36_medf36_t","cnp_pvoc_pvocrtcr","cnp_pvrt_pvrtrtcr","cnp_svdelay_svtldrtc","stroop_total_mean","latentgroove","Latencycomposite")

allefvars<-c("cnp_pcet_pcet_acc2","cnp_pcet_pcetrtcr","cnp_sfnb2_sfnb_mcr","cnp_sfnb2_sfnb_mrtc","cnp_spcptnl_scpt_tp","cnp_spcptnl_scpt_tprt","stroop_total_mean")

accvars<-c("cnp_pcet_pcet_acc2","cnp_sfnb2_sfnb_mcr","cnp_spcptnl_scpt_tp")
latvars<-c("cnp_pcet_pcetrtcr","cnp_sfnb2_sfnb_mrtc","cnp_spcptnl_scpt_tprt","stroop_total_mean")

eftestgroups<-c("cnp_pcet","cnp_sfnb2","cnp_spcptnl") 
####################
###composites#######
coglongdata$Latencycomposite<-compositecols(latvars,coglongdata)
coglongdata$Accuracycomposite<-compositecols(accvars,coglongdata)

write.csv(coglongdata,"~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/NCANDA/Data/btc_NCANDAscoredmeasures_20220306.outlierremoved.compositeacclat.csv")

adultidx<-which(coglongdata$cnp_age > 20 & coglongdata$cnp_age < 30 & coglongdata$cnp_age < max(coglongdata$cnp_age))
unique(coglongdata$subject[adultidx])
coglongdata<-scalebyidxrenamereturndf(df=coglongdata,vars=c(allefvars,"Accuracycomposite","Latencycomposite"),adultidx,appendname="20_30")
write.csv(coglongdata,"~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/NCANDA/data/btc_R03cleaneddata_20220306.outlierremoved.compositeacclat.2030scale.csv")


########
#########
