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
library(gamm4)
######functions#####
source("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/GeneralScripts/outliers_mgcvgam.R")
############
NCANDAdata<-read.csv("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/NCANDA/Data/btc_NCANDAscoredmeasures_20220214.csv")
####vars############
####################
allaccvars<-c("cnp_cpf_ifac_tot","cnp_cpw_iwrd_tot","cnp_spcptnl_scpt_tp","cnp_sfnb2_sfnb_mcr","cnp_pmat24a_pmat24_a_cr","cnp_cpfd_dfac_tot","cnp_cpwd_dwrd_tot",
              "cnp_shortvolt_svt","cnp_er40d_er40_cr","cnp_pcet_pcet_acc2","cnp_medf36_medf36_a","cnp_pvoc_pvoccr","cnp_pvrt_pvrt_pc","cnp_svdelay_svt_ld","np_wais4_rawscore_computed","latentdd","Accuracycomposite")   

alllatvars<-c("cnp_cpf_ifac_rtc","cnp_cpw_iwrd_rtc","cnp_spcptnl_scpt_tprt","cnp_sfnb2_sfnb_mrtc","cnp_pmat24a_pmat24_a_rtcr","cnp_cpfd_dfac_rtc",
              "cnp_cpwd_dwrd_rtc","cnp_shortvolt_svtcrt","cnp_er40d_er40_crt","cnp_pcet_pcetrtcr","cnp_medf36_medf36_t","cnp_pvoc_pvocrtcr","cnp_pvrt_pvrtrtcr","cnp_svdelay_svtldrtc","stroop_total_mean","latentgroove","Latencycomposite")

allefvars<-c("cnp_pcet_pcet_acc2","cnp_pcet_pcetrtcr","cnp_sfnb2_sfnb_mcr","cnp_sfnb2_sfnb_mrtc","cnp_spcptnl_scpt_tp","cnp_spcptnl_scpt_tprt","stroop_total_mean")

accvars<-c("cnp_pcet_pcet_acc2","cnp_sfnb2_sfnb_mcr","cnp_spcptnl_scpt_tp")
latvars<-c("cnp_pcet_pcetrtcr","cnp_sfnb2_sfnb_mrtc","cnp_spcptnl_scpt_tprt","stroop_total_mean")

eftestgroups<-c("cnp_pcet","cnp_sfnb2","cnp_spcptnl") 
####################
coglongwinds<-NCANDAdata
coglongwinds<-coglongwinds[!is.na(coglongwinds$cnp_age),]
###visits#######
coglongwinds[,allefvars]<-lapply(coglongwinds[,allefvars],mgcvgamleveragedetect,covar=coglongwinds$cnp_age,idvar="subject",df=coglongwinds)
mvouts<-scale(psych::outlier(coglongwinds[,allefvars]))
coglongmdcensor<-coglongwinds
coglongmdcensor$mvout<-0
coglongmdcensor$mvout[mvouts>4]<-1
coglongmdcensor<-coglongmdcensor[coglongmdcensor$mvout==0,]
###27 visits removed entierly#########
coglongwinds_means<-coglongmdcensor %>% dplyr::select(c(allefvars,"subject")) %>% group_by(subject) %>% dplyr::summarise_all(~mean(.,na.rm=TRUE))
#coglongwinds_means[,allefvars]<-abs(scale(coglongwinds_means[,allefvars]))
coglongwinds_means_vars<-coglongwinds_means[,allefvars]
#coglongwinds_means_vars[coglongwinds_means_vars>4]<-NA #####integrated into distribute function
coglongwinds_means_varswitid<-coglongwinds_means_vars
coglongwinds_means_varswitid$subject<-coglongwinds_means$subject

###
coglongmdcensorwithsubcensor<-coglongmdcensor
coglongmdcensorwithsubcensor<-distributeNAsforoutliermeans(coglongwinds_means_varswitid,coglongmdcensorwithsubcensor,allefvars,idvar="subject")
coglongmdcensorwithsubcensor<-coglongmdcensorwithsubcensor[rowSums(is.na(coglongmdcensorwithsubcensor[,allefvars])) != ncol(coglongmdcensorwithsubcensor[,allefvars]), ]
###complete visits 3412 total visits#####
##########write data########
write.csv(coglongmdcensorwithsubcensor,"~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/NCANDA//Data/btc_NCANDAscoredmeasures_20220214.outlierremoved.csv")

#for(x in 1:1000){print(sprintf("each day you get %s times better",x))}







