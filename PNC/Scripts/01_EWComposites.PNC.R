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
########
source("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/GeneralScripts/metabyage.R")
coglongdata<-read.csv("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/PNC/Data/btc_PNCscoredmeasures_20220214.outlierremoved.csv")
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

####################
###composites#######
coglongdata$Latencycomposite<-compositecols(latvars,coglongdata)
coglongdata$Accuracycomposite<-compositecols(accvars,coglongdata)
write.csv(coglongdata,"~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/PNC/Data/btc_PNCscoredmeasures_20220306.outlierremoved.compositeacclat.csv")

adultidx<-which(coglongdata$age > 20 & coglongdata$age < 30 & coglongdata$age < max(coglongdata$age))
coglongdata<-scalebyidxrenamereturndf(df=coglongdata,vars=c(allefvars,"Accuracycomposite","Latencycomposite"),adultidx,appendname="20_30")
write.csv(coglongdata,"~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/PNC/data/btc_R03cleaneddata_20220306.outlierremoved.compositeacclat.2030scale.csv")

