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
#####0X NCDANDA cognitive data#######
###get cognitive data, remove bad records, SAVE#######
###Stroop #######
stroop<-read.csv("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/NCANDA_Behavioral/Data/NCANDA_RELEASE_4Y_REDCAP_MEASUREMENTS_V02/summaries/redcap/stroop.csv")
stroop_merge<-stroop[,c("subject","visit","stroop_total_mean")]
#########CNP##########
cnpdatadict<-read.csv("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/NCANDA_Behavioral/Data/NCANDA_RELEASE_4Y_REDCAP_MEASUREMENTS_V02///datadict/redcap/cnp_datadict.csv")
cnpdata<-read.csv("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/NCANDA_Behavioral/Data/NCANDA_RELEASE_4Y_REDCAP_MEASUREMENTS_V02/summaries/redcap/cnp.csv")
########
allaccvars<-c("cnp_cpf_ifac_tot","cnp_cpw_iwrd_tot","cnp_spcptnl_scpt_tp","cnp_sfnb2_sfnb_mcr","cnp_pmat24a_pmat24_a_cr","cnp_cpfd_dfac_tot","cnp_cpwd_dwrd_tot",
              "cnp_shortvolt_svt","cnp_er40d_er40_cr","cnp_pcet_pcet_acc2","cnp_medf36_medf36_a","cnp_pvoc_pvoccr","cnp_pvrt_pvrt_pc","cnp_svdelay_svt_ld","np_wais4_rawscore_computed","latentdd","Accuracycomposite")   

alllatvars<-c("cnp_cpf_ifac_rtc","cnp_cpw_iwrd_rtc","cnp_spcptnl_scpt_tprt","cnp_sfnb2_sfnb_mrtc","cnp_pmat24a_pmat24_a_rtcr","cnp_cpfd_dfac_rtc",
              "cnp_cpwd_dwrd_rtc","cnp_shortvolt_svtcrt","cnp_er40d_er40_crt","cnp_pcet_pcetrtcr","cnp_medf36_medf36_t","cnp_pvoc_pvocrtcr","cnp_pvrt_pvrtrtcr","cnp_svdelay_svtldrtc","stroop_total_mean","latentgroove","Latencycomposite")

allefvars<-c("cnp_pcet_pcet_acc2","cnp_pcet_pcetrtcr","cnp_sfnb2_sfnb_mcr","cnp_sfnb2_sfnb_mrtc","cnp_spcptnl_scpt_tp","cnp_spcptnl_scpt_tprt","stroop_total_mean")

accvars<-c("cnp_pcet_pcet_acc2","cnp_sfnb2_sfnb_mcr","cnp_spcptnl_scpt_tp")
latvars<-c("cnp_pcet_pcetrtcr","cnp_sfnb2_sfnb_mrtc","cnp_spcptnl_scpt_tprt","stroop_total_mean")

eftestgroups<-c("cnp_pcet","cnp_sfnb2","cnp_spcptnl") 

efficiency_function<-function(df,accvar,latvar){
  print("accuracy and latency correlated at:")
  acclatcor<-cor(as.numeric(as.character(df[,accvar])),as.numeric(as.character(df[,latvar])),use="complete")
  print(acclatcor)
  corrsign<-sign(acclatcor)
  print(sprintf("efficiency based on %s sign",corrsign))
  eff<-scale(as.numeric(as.character(df[,accvar])),center=TRUE,scale=TRUE)+scale(corrsign*as.numeric(as.character(df[,latvar])),center=TRUE,scale=TRUE) ####lat sign flipped so higher means "better" on both scales
  return(eff)
}

allefficiency<-lapply(eftestgroups,function(vg){
  print(vg)
  accvar<-grep(vg,allaccvars,value=TRUE)
  nacc<-length(which(!is.na(cnpdata[,accvar])))
  latvar=grep(vg,alllatvars,value=TRUE)
  nlat<-length(which(!is.na(cnpdata[,latvar])))
  print(sprintf("Latency has %s useable datapoints; Accuracy has %s usable datapoints",nlat,nacc))
  eff<-efficiency_function(cnpdata,accvar=accvar,latvar=latvar)
  effdf<-as.data.frame(eff)
  names(effdf)<-paste0(vg,"_eff")
  return(effdf)
})

alleffs<-as.data.frame(do.call(cbind,allefficiency))
cnpdatawitheffs<-cbind(cnpdata,alleffs)

allNCANDA<-merge(cnpdatawitheffs,stroop_merge,by=c("subject","visit"))

write.csv(allNCANDA,"~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/NCANDA/Data/btc_NCANDAscoredmeasures_20220214.csv")

