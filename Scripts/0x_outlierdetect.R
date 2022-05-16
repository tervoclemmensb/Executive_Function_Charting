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
######functions#####
source("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/GeneralScripts/outliers_mgcvgam.R")
############
coglongdata<-read.csv("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Data/btc_R03cleaneddata_20220214.csv")
coglongdata<-as.data.frame(coglongdata)
####vars############
####################
allefvars<-c("Anti_CRLat","Anti_CRR","Mix_CRLat","Mix_CRR","SOC.Problems.solved.in.minimum.moves",
                 "SOC.Overallmeaninitialthinkingtime","DMS.PC","DMS.Median.correct.latency","SSP.Span.length",
                 "nfixbreak","best_acc_m_exclude","first_lat_m_exclude")

accvars<-c("Anti_CRR","Mix_CRR","DMS.Percent.correct","SSP.Span.length","nfixbreak")
accvarsperc<-c("Anti_CRR","Mix_CRR","DMS.PC")
accvarscount<-c("SSP.Span.length","nfixbreak")
##convert DMS to perc
coglongdata$DMS.PC<-coglongdata$DMS.Percent.correct/100
###add visit column###
coglongdata<-coglongdata %>% group_by(id) %>% mutate(visit=rank(d8))
#coglongdata<-coglongdata[coglongdata$id!=10408,]###
coglongdata<-as.data.frame(coglongdata)
####################
coglongwinds<-coglongdata
###visits#######
coglongwinds$id<-as.factor(coglongwinds$id)
coglongwinds[,allefvars]<-lapply(coglongwinds[,allefvars],mgcvgamleveragedetect,covar=coglongwinds$Ageatvisit,idvar="id",df=coglongwinds,
                                 unioutfirst=FALSE)
mvouts<-scale(psych::outlier(coglongwinds[,allefvars]))
coglongmdcensor<-coglongwinds
coglongmdcensor$mvout<-0
coglongmdcensor$mvout[mvouts>4]<-1
coglongmdcensor<-coglongmdcensor[coglongmdcensor$mvout==0,]
###visits removed entierly#########
####invidiaul subjects########
coglongwinds_means<-coglongmdcensor %>% dplyr::select(c(allefvars,"id")) %>% group_by(id) %>% dplyr::summarise_all(~mean(.,na.rm=TRUE))
#coglongwinds_means[,allefvars]<-abs(scale(coglongwinds_means[,allefvars]))
coglongwinds_means_vars<-coglongwinds_means[,allefvars]
#coglongwinds_means_vars[coglongwinds_means_vars>4]<-NA #####integrated into distribute function
coglongwinds_means_varswitid<-coglongwinds_means_vars
coglongwinds_means_varswitid$id<-coglongwinds_means$id

###
coglongmdcensorwithsubcensor<-coglongmdcensor
coglongmdcensorwithsubcensor<-distributeNAsforoutliermeans(coglongwinds_means_varswitid,coglongmdcensorwithsubcensor,allefvars,idvar="id")
coglongmdcensorwithsubcensor<-coglongmdcensorwithsubcensor[rowSums(is.na(coglongmdcensorwithsubcensor[,allefvars])) != ncol(coglongmdcensorwithsubcensor[,allefvars]), ]
###196 participants 668 total visits#####
##########write data########
write.csv(coglongmdcensorwithsubcensor,"~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/data/btc_R03cleaneddata_20220214.outlierremoved.csv")








