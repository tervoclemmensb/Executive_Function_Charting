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
coglongdata<-read.csv("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/data/btc_R03cleaneddata_20220214.outlierremoved.csv")
####vars############
####################
allefvars<-c("Anti_CRLat","Anti_CRR","Mix_CRLat","Mix_CRR","SOC.Problems.solved.in.minimum.moves",
             "SOC.Overallmeaninitialthinkingtime","DMS.PC","DMS.Median.correct.latency","SSP.Span.length",
             "nfixbreak","best_acc_m_exclude","first_lat_m_exclude")

accvars<-c("Anti_CRR","Mix_CRR","DMS.PC","SSP.Span.length","nfixbreak_fl","SOC.Problems.solved.in.minimum.moves","best_acc_m_exclude_fl")

latvars<-c("Anti_CRLat","Mix_CRLat","SOC.Overallmeaninitialthinkingtime","DMS.Median.correct.latency","first_lat_m_exclude")

allefvars<-c(accvars,latvars)

coglongdata$visitnum<-as.numeric(coglongdata$visit)
###covert acc measures to be positive=better####
coglongdata$best_acc_m_exclude_fl<-coglongdata$best_acc_m_exclude*-1
coglongdata$nfixbreak_fl<-coglongdata$nfixbreak*-1
####################
###composites#######
coglongdata$Latencycomposite<-compositecols(latvars,coglongdata)
coglongdata$Accuracycomposite<-compositecols(accvars,coglongdata)
####scale vars to adults for cross-data harmonization with different age ranges####
write.csv(coglongdata,"~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/data/btc_R03cleaneddata_20220306.outlierremoved.compositeacclat.csv")
adultidx<-which(coglongdata$Ageatvisit > 20 & coglongdata$Ageatvisit < 30 & coglongdata$Ageatvisit < max(coglongdata$Ageatvisit))
unique(coglongdata$subj[adultidx])
coglongdata<-scalebyidxrenamereturndf(df=coglongdata,vars=c(allefvars,"Accuracycomposite","Latencycomposite"),adultidx,appendname="20_30")
write.csv(coglongdata,"~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/data/btc_R03cleaneddata_20220306.outlierremoved.compositeacclat.2030scale.csv")




