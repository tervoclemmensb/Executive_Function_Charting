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
coglongdata<-read.csv("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/NKI/Data/btc_NKIscoredmeasures_20220208.outlierremoved.csv")
####vars############
####################
accvars<-c("PCET_PCET_ACC2","SLNB2_SLNB_MCR","SPCPTNL_SCPT_TP","TOWacc","DFLacc")
latvars<-c("PCET_PCETRTCR","SLNB2_SLNB_MRTC","SPCPTNL_SCPT_TPRT","CWIlat","TMTlat")

allefvars<-c(accvars,latvars)

####################
###composites#######
coglongdata$Latencycomposite<-compositecols(latvars,coglongdata)
coglongdata$Accuracycomposite<-compositecols(accvars,coglongdata)
write.csv(coglongdata,
          file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/NKI/Data/btc_NKIscoredmeasures_20220306.outlierremoved.compositeacclat.csv")

########
adultidx<-which(coglongdata$age > 20 & coglongdata$age < 30 & coglongdata$age < max(coglongdata$age))
coglongdata<-scalebyidxrenamereturndf(df=coglongdata,vars=c(allefvars,"Accuracycomposite","Latencycomposite"),adultidx,appendname="20_30")
write.csv(coglongdata,"~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/NKI/data/btc_NKIscoredmeasures_20220306.outlierremoved.compositeacclat.2030scale.csv")





