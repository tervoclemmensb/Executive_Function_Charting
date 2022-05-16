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
source("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/GeneralScripts/outliers_mgcvgam.R")
############
NKIdata<-read.csv("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/NKI/Data/btc_NKIscoredmeasures_20220208.csv")
####vars############
####################
accvars<-c("PCET_PCET_ACC2","SLNB2_SLNB_MCR","SPCPTNL_SCPT_TP","TOWacc","DFLacc")
latvars<-c("PCET_PCETRTCR","SLNB2_SLNB_MRTC","SPCPTNL_SCPT_TPRT","CWIlat","TMTlat")

allefvars<-c(accvars,latvars)

eftestgroups<-c("PCET_PCET","SLNB2_SLNB","SPCPTNL_SCPT") 

####################
coglongwinds<-NKIdata
coglongwinds<-coglongwinds[!is.na(coglongwinds$age),]
###visits#######
coglongwinds[,allefvars]<-lapply(coglongwinds[,allefvars],mgcvgamleveragedetect,covar=coglongwinds$age,idvar="subject",df=coglongwinds)
########
mvouts<-scale(psych::outlier(coglongwinds[,allefvars]))
coglongmdcensor<-coglongwinds
coglongmdcensor$mvout<-0
coglongmdcensor$mvout[mvouts>4]<-1
coglongmdcensor<-coglongmdcensor[coglongmdcensor$mvout==0,]
coglongmdcensornoallna<-coglongmdcensor[rowSums(is.na(coglongmdcensor[,allefvars])) != ncol(coglongmdcensor[,allefvars]), ]
###588 participants 
##########write data########
write.csv(coglongmdcensornoallna,"~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/NKI/Data/btc_NKIscoredmeasures_20220208.outlierremoved.csv")





