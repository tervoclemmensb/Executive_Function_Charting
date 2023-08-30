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
source("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/GeneralScripts/metabyage.R")
#####CANTAB############################
CANTAB<-read.csv("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Data/CogCANTAB_20190128.csv")
CANTAB<-CANTAB[CANTAB$Age>7,] ##remove datapoitns with invalid ages
CANTAB<-CANTAB[CANTAB$Age<=35,] ##remove datapoitns with invalid ages
CANTAB_nowarn<-CANTAB[,!(names(CANTAB) %in% grep("Warning",names(CANTAB),value=TRUE))]
CANTAB_nowarn$SOC.Overallmeaninitialthinkingtime<- apply(CANTAB_nowarn[c("SOC.Mean.initial.thinking.time..2.moves.","SOC.Mean.initial.thinking.time..3.moves.","SOC.Mean.initial.thinking.time..4.moves.","SOC.Mean.initial.thinking.time..5.moves.")], 1, median)
CANTAB_nowarn$SOC.Overallmeansubsequentthinkingtime<-apply(CANTAB_nowarn[c("SOC.Mean.subsequent.thinking.time..2.moves.","SOC.Mean.subsequent.thinking.time..3.moves.","SOC.Mean.subsequent.thinking.time..4.moves.","SOC.Mean.subsequent.thinking.time..5.moves.")],1, median)
CANTAB_nowarn$d8<-as.numeric(format(mdy_hms(CANTAB_nowarn$Session.start.time),"%Y%m%d"))
CANTAB_nowarn$id<-gsub("_.","",gsub("[xX]([0-9])","",CANTAB_nowarn$Subject.ID))
###gender write from CANTAB
GENDER_R<-CANTAB_nowarn[,c("id","Gender")]
write.csv(GENDER_R,"~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Data/gender_20220214.csv")
#####
#MOT Median Latency MOT Mean Error###
###SOC.Problems.solved.in.minimum.moves" #
###MEan SOC.Mean.initial.thinking.time..3.moves###
#"DMS.Median.correct.latency"    
#DMS.Percent.correct
#SSP.Span.length

###Measures
CANTAB_Tmeasures<-CANTAB_nowarn[,c("id","Age","d8","MOT.Median.latency","MOT.Mean.error","SOC.Problems.solved.in.minimum.moves","SOC.Overallmeaninitialthinkingtime","SOC.Overallmeansubsequentthinkingtime","DMS.Percent.correct","DMS.Median.correct.latency","SSP.Span.length","SSP.Mean.time.to.first.response","SSP.Mean.time.to.last.response")]
CANTAB_Tmeasures<-CANTAB_Tmeasures[!duplicated(CANTAB_Tmeasures[,c("id","d8")]),]
CANTABdemo<-CANTAB_nowarn[,c("Age","SOC.Problems.solved.in.minimum.moves","DMS.Percent.correct","SSP.Span.length")]
CANTABdemo_d<-CANTABdemo[,c("SOC.Problems.solved.in.minimum.moves","DMS.Percent.correct","SSP.Span.length")]

###########################
####MGS scored#####
###data pull from server --> /Volumes/Hera/Projects/autoeyescore/mgs/
MGS<-read.csv("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Data//btc_summaries_20200817.csv")
MGS$nfixbreakrate<-MGS$nfixbreak/MGS$ntrial
MGS$id<-MGS$subj
MGS$d8<-MGS$date
MGS$MGS_ntrials<-MGS$ntrial
###Eye Tracking ############
library(dplyr)
library(tidyr)
Eyedata<-read.csv("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Data/all_et_score.csv")
Eyedata<-Eyedata[Eyedata$lat>60,] ###limit trials where there is impossibly short saccade (60ms)
Eyedata<-Eyedata[!is.na(Eyedata$id),]
Eyedata$scoringtype<-"auto"
Eyedata$scoringtype[grep("manual",Eyedata$fn)]<-"manual"
###all taks but mix
Eyedata_wide <- 
  Eyedata %>%
  filter(Count!=-1,task!="Mix") %>% # remove drop,no mix
  # summarise
  group_by(id,d8,vtype,task,scoringtype) %>%
  summarise(CRR=length(which(Count==1))/n(),
            CRLat=median(lat[Count==1]),
            Ageatvisit=mean(age),
            nvisits=length(unique(age)),
            ntrials=n()) %>% 
  gather(metric,value,  CRR, CRLat,ntrials) %>%
  unite(task_metric,task,metric) %>%
  spread(task_metric, value )
###mix####
##as from mix####
Eyedata_widemixAS <- 
  Eyedata %>%
  filter(Count!=-1,task=="Mix",AS=="AS") %>% # remove drop,only mix, only antisaccade trials
  # summarise
  group_by(id,d8,vtype,task,scoringtype,AS) %>%
  summarise(CRR=length(which(Count==1))/n(),
            CRLat=median(lat[Count==1]),
            ntrials=n()) %>% 
  gather(metric,value,  CRR, CRLat,ntrials) %>%
  unite(task_metric,task,metric) %>%
  spread(task_metric,value )

########
Eyedata_wide<-Eyedata_wide[Eyedata_wide$Ageatvisit<45,] ###remove incorrectly coded age visits (> oldest possible age in study)
Eyedata_wide$antiminusvgslat<-Eyedata_wide$Anti_CRLat-Eyedata_wide$VGS_CRLat
Eyedata_wide$antiminusvgslat_normvgslat<-(Eyedata_wide$Anti_CRLat-Eyedata_wide$VGS_CRLat)/Eyedata_wide$VGS_CRLat
Eyedata_wide$antiminusvgslat_normantiplusvgs<-(Eyedata_wide$Anti_CRLat-Eyedata_wide$VGS_CRLat)/(Eyedata_wide$Anti_CRLat+Eyedata_wide$VGS_CRLat)

Eyedata_wide$iddate<-paste(Eyedata_wide$id,Eyedata_wide$d8,sep="_")
#######check manual and autoscoring#######
a<-reshape2::melt(Eyedata_wide, id.vars=c("iddate","vtype","scoringtype"))
m <-
  a %>% filter(scoringtype=='auto') %>%
  merge(a %>% filter(scoringtype=='manual'),
        by=c("iddate","vtype","variable"),
        suffixes=c('.auto','.man') )
m %>% group_by(variable) %>% summarise(cor(value.auto,value.man, use='pairwise.complete.obs',method="kendall"))

autoeyescorevalidationplot<-ggplot(m %>% filter(variable %in% c("Anti_CRR","Anti_CRLat")),aes(x=value.man,y=value.auto,colour=variable))+geom_point()+facet_wrap(~variable,scales="free")
autoeyescorevalidationplot<-LNCDR::lunaize(autoeyescorevalidationplot)+xlab("\nManual Scoring")+ylab("Automatic Scoring\n")+theme(legend.position = "none") 
###########
Eyedata_wide_auto<-m[m$scoringtype=="auto",]
Eyedata_wide_manual<-m[m$scoringtype=="manual",]
eyedata_wide_automanualtest<-merge(m,m,by=c("iddate","vtype"),suffixes=c(".auto",".man"))
#########
#########################merge all data####
Eyedata_wide_autoonly<-Eyedata_wide[Eyedata_wide$scoringtype=="auto",]
Eyedata_wide_autoonlywithmix<-merge(Eyedata_wide_autoonly,Eyedata_widemixAS,by=c("id","d8","vtype","scoringtype"),all=TRUE)

Eyedata_wide_autoonlywithcantab<-merge(Eyedata_wide_autoonlywithmix,CANTAB_Tmeasures,by=c("id","d8"),all=TRUE)

allcoglongmeasures<-merge(Eyedata_wide_autoonlywithcantab,MGS,by=c("id","d8"),all = TRUE)
allcoglongmeasures<-allcoglongmeasures[as.numeric(allcoglongmeasures$id)>10000 & as.numeric(allcoglongmeasures$id)<20000,] ###limit ids to expected values (remove RA's coded test subjects)
allcoglongmeasures<-allcoglongmeasures[!is.na(allcoglongmeasures$id) & !is.na(allcoglongmeasures$d8),] ###limit to records with a real id and d8

allcoglongmeasuresminimumtrials<-allcoglongmeasures[allcoglongmeasures$Anti_ntrials>.7*max(allcoglongmeasures$Anti_ntrials,na.rm=TRUE) & allcoglongmeasures$Mix_ntrials>.7*max(allcoglongmeasures$Mix_ntrials,na.rm=TRUE) & allcoglongmeasures$MGS_ntrials>.7*max(allcoglongmeasures$MGS_ntrials,na.rm=TRUE),]
allcoglongmeasuresminimumtrials<-allcoglongmeasuresminimumtrials[!is.na(allcoglongmeasuresminimumtrials$id) & !is.na(allcoglongmeasuresminimumtrials$d8),] ###limit to records with a real id and d8

write.csv(allcoglongmeasuresminimumtrials,"~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Data/btc_R03cleaneddata_20220214.csv")





