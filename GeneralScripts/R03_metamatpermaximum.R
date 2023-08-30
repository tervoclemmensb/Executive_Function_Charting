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
library(cowplot)
library(gamm4)
library(parallel)
library(scales)
library(RColorBrewer)
###
source("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/GeneralScripts/growthrate_mgcvgam.R")
source("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/GeneralScripts/metabyage.R")
###Read fit data######adult scaled (ultimately this analysis scales to maximum so orig scale doesn't matter but using to minimize different versions of same data)
####
###fixed age grid###
agegridforinterp<-data.frame(ages=as.numeric(seq(8,35,by=.1)))
####Read common adult scale fits######## 
load("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/AdultScaledFits/LUNAallfitsadultscaled.Rdata")
LUNAscaledfits<-scalefitsall2030sig
lunanonsig<-unique(LUNAscaledfits$outcome[LUNAscaledfits$sig!="1"])
acclatmeasures<-LUNAscaledfits %>% dplyr::group_by(outcome) %>% dplyr::summarize(type=unique(type))
acclatmeasures$var<-acclatmeasures$outcome
LUNAinterpolatedfits<-interpolatebyage(LUNAscaledfits,agegrid=agegridforinterp,agegridvar="ages",agevar="pred",vars=unique(LUNAscaledfits$outcome),longformat = TRUE,
                                       datanamefromlongformat="fit",varnameforlongformat="outcome",namemodifier="fit",returnlongformat=TRUE)
LUNAinterpolatedfitsse<-interpolatebyage(LUNAscaledfits,agegrid=agegridforinterp,agegridvar="ages",agevar="pred",vars=unique(LUNAscaledfits$outcome),longformat = TRUE,
                                       datanamefromlongformat="se",varnameforlongformat="outcome",namemodifier="se",returnlongformat=TRUE)
LUNAallinterpfits<-merge(LUNAinterpolatedfits,LUNAinterpolatedfitsse,by=c("ages","var"))
LUNAallinterpfits$dataset<-"LUNA"
LUNAallinterpfits<-merge(LUNAallinterpfits,acclatmeasures,by=c("var"))
#####
load("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/AdultScaledFits/NCANDAallfitsadultscaled.Rdata")
NCANDAscaledfits<-scalefitsall2030sig
ncandanonsig<-unique(NCANDAscaledfits$outcome[NCANDAscaledfits$sig!="1"])
acclatmeasures<-NCANDAscaledfits %>% dplyr::group_by(outcome) %>% dplyr::summarize(type=unique(type))
acclatmeasures$var<-acclatmeasures$outcome
NCANDAinterpolatedfits<-interpolatebyage(NCANDAscaledfits,agegrid=agegridforinterp,agegridvar="ages",agevar="pred",vars=unique(NCANDAscaledfits$outcome),longformat = TRUE,
                                       datanamefromlongformat="fit",varnameforlongformat="outcome",namemodifier="fit",returnlongformat=TRUE)
NCANDAinterpolatedfitsse<-interpolatebyage(NCANDAscaledfits,agegrid=agegridforinterp,agegridvar="ages",agevar="pred",vars=unique(NCANDAscaledfits$outcome),longformat = TRUE,
                                         datanamefromlongformat="se",varnameforlongformat="outcome",namemodifier="se",returnlongformat=TRUE)
NCANDAallinterpfits<-merge(NCANDAinterpolatedfits,NCANDAinterpolatedfitsse,by=c("ages","var"))
NCANDAallinterpfits$dataset<-"NCANDA"
NCANDAallinterpfits<-merge(NCANDAallinterpfits,acclatmeasures,by=c("var"))
####
load("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/AdultScaledFits/NKIallfitsadultscaled.Rdata")
NKIscaledfits<-scalefitsall2030sig
nkinonsig<-unique(NKIscaledfits$outcome[NKIscaledfits$sig!="1"])
acclatmeasures<-NKIscaledfits %>% dplyr::group_by(outcome) %>% dplyr::summarize(type=unique(type))
acclatmeasures$var<-acclatmeasures$outcome
NKIinterpolatedfits<-interpolatebyage(NKIscaledfits,agegrid=agegridforinterp,agegridvar="ages",agevar="pred",vars=unique(NKIscaledfits$outcome),longformat = TRUE,
                                         datanamefromlongformat="fit",varnameforlongformat="outcome",namemodifier="fit",returnlongformat=TRUE)
NKIinterpolatedfitsse<-interpolatebyage(NKIscaledfits,agegrid=agegridforinterp,agegridvar="ages",agevar="pred",vars=unique(NKIscaledfits$outcome),longformat = TRUE,
                                           datanamefromlongformat="se",varnameforlongformat="outcome",namemodifier="se",returnlongformat=TRUE)
NKIallinterpfits<-merge(NKIinterpolatedfits,NKIinterpolatedfitsse,by=c("ages","var"))
NKIallinterpfits$dataset<-"NKI"
NKIallinterpfits<-merge(NKIallinterpfits,acclatmeasures,by=c("var"))
###

load("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/AdultScaledFits/PNCallfitsadultscaled.Rdata")
PNCscaledfits<-scalefitsall2030sig
pncnonsig<-unique(PNCscaledfits$outcome[PNCscaledfits$sig!="1"])
acclatmeasures<-PNCscaledfits %>% dplyr::group_by(outcome) %>% dplyr::summarize(type=unique(type))
acclatmeasures$var<-acclatmeasures$outcome
PNCinterpolatedfits<-interpolatebyage(PNCscaledfits,agegrid=agegridforinterp,agegridvar="ages",agevar="pred",vars=unique(PNCscaledfits$outcome),longformat = TRUE,
                                      datanamefromlongformat="fit",varnameforlongformat="outcome",namemodifier="fit",returnlongformat=TRUE)
PNCinterpolatedfitsse<-interpolatebyage(PNCscaledfits,agegrid=agegridforinterp,agegridvar="ages",agevar="pred",vars=unique(PNCscaledfits$outcome),longformat = TRUE,
                                        datanamefromlongformat="se",varnameforlongformat="outcome",namemodifier="se",returnlongformat=TRUE)
PNCallinterpolatedfits<-merge(PNCinterpolatedfits,PNCinterpolatedfitsse,by=c("ages","var"))
PNCallinterpolatedfits$dataset<-"PNC"
PNCallinterpolatedfits<-merge(PNCallinterpolatedfits,acclatmeasures,by=c("var"))
#############
############
allinterpolatedfits<-plyr::rbind.fill(LUNAallinterpfits,NCANDAallinterpfits) %>% plyr::rbind.fill(.,NKIallinterpfits) %>% plyr::rbind.fill(.,PNCallinterpolatedfits)

allinterpolatedfits$varglobal<-paste(allinterpolatedfits$var,allinterpolatedfits$dataset,sep="_")

allfits<-plyr::rbind.fill(LUNAscaledfits,NCANDAscaledfits) %>% plyr::rbind.fill(.,NKIscaledfits) %>% plyr::rbind.fill(.,PNCscaledfits)

#plyr::rbind.fill(LUNAallinterpfits,NCANDAscaledfits) %>% plyr::rbind.fill(.,NKIscaledfits) %>% plyr::rbind.fill(.,PNCscaledfits)

###link Penn tasks for meta basis
##acc
PNBKacc<-c("SLNB2_SLNB_MCR","LNB_MCR","cnp_sfnb2_sfnb_mcr")
PNBKacc<-paste0(PNBKacc,"_20_30")
PCETacc<-c("cnp_pcet_pcet_acc2","PCET_PCET_ACC2","PCET_ACC2")
PCETacc<-paste0(PCETacc,"_20_30")
PCPTacc<-c("cnp_spcptnl_scpt_tp","PCPT_T_TP","SPCPTNL_SCPT_TP")
PCPTacc<-paste0(PCPTacc,"_20_30")

allinterpolatedfits$varglobal[allinterpolatedfits$var %in% PNBKacc]<-"PNBKacc"
allinterpolatedfits$varglobal[allinterpolatedfits$var %in% PCETacc]<-"PCETacc"
allinterpolatedfits$varglobal[allinterpolatedfits$var %in% PCPTacc]<-"PCPTacc"
##lat###
PNBKlat<-c("cnp_sfnb2_sfnb_mrtc","LNB_MRTC","SLNB2_SLNB_MRTC")
PNBKlat<-paste0(PNBKlat,"_20_30")
PCETlat<-c("cnp_pcet_pcetrtcr","PCET_RTCR","PCET_PCETRTCR")
PCETlat<-paste0(PCETlat,"_20_30")
PCPTlat<-c("cnp_spcptnl_scpt_tprt","PCPT_T_TPRT","SPCPTNL_SCPT_TPRT")
PCPTlat<-paste0(PCPTlat,"_20_30")

allinterpolatedfits$varglobal[allinterpolatedfits$var %in% PNBKlat]<-"PNBKlat"
allinterpolatedfits$varglobal[allinterpolatedfits$var %in% PCETlat]<-"PCETlat"
allinterpolatedfits$varglobal[allinterpolatedfits$var %in% PCPTlat]<-"PCPTlat"

####plot percent of maximum  for those significant age-related change 
###calculating this metric in measures with no change would  bias estimates to appear to mature earlier
allinterpolatedfitssig<-allinterpolatedfits[!(allinterpolatedfits$outcome %in% c(lunanonsig,ncandanonsig)),] ###all sig in NKI and PNC
#allinterpolatedfitssig<-allinterpolatedfitssig[!(allinterpolatedfitssig$dataset %in% c("PNC")),] ###sensitivity removing PNC dataset

accmetafitsig<-metabyagethreelevel(interpolatedagedfs=allinterpolatedfitssig[allinterpolatedfitssig$type=="acc",],agevar="ages",valcol="fit_value",secol ="se_value",variablenestvar="varglobal",datasetvar = "dataset") ###three level meta to define basis function
accmetafitsig$percentmax<-scales::rescale(accmetafitsig$estimate)
latmetafitsig<-metabyagethreelevel(interpolatedagedfs=allinterpolatedfitssig[allinterpolatedfitssig$type=="lat",],agevar="ages",valcol="fit_value",secol ="se_value",variablenestvar="varglobal",datasetvar = "dataset") ###three level meta to define basis function
latmetafitsig$percentmax<-scales::rescale(latmetafitsig$estimate)

max(accmetafitsig$percentmax[accmetafitsig$age<=18],na.rm=TRUE)
max(accmetafitsig$percentmax[accmetafitsig$age<=20],na.rm=TRUE)

1-min(latmetafitsig$percentmax[latmetafitsig$age<=18],na.rm=TRUE) ###subtract by 1 because decreasing
1-min(latmetafitsig$percentmax[latmetafitsig$age<=20],na.rm=TRUE) ###substract by 1 because decreasing 

allfits$varbydataset<-paste(allfits$outcome,allfits$dataset)
allfitsacc<-allfits[allfits$outcome %in% unique(allinterpolatedfitssig$outcome) & allfits$type=="acc",]

accmetafitsigsmooth<-mgcvscalefits(accmetafitsig,outcomevars = c("percentmax"),predvars = "age",mformula = as.formula("outcome~s(pred)"),scale=FALSE,interval_inc = .1) ###fit smooth

ggderivpercentfit<-ggplot(allfitsacc,aes(x=pred,y=fitscale,colour=varbydataset))+geom_line(size=.5,alpha=.35)+
  theme(legend.position = "none")+
  geom_line(data=accmetafitsigsmooth,aes(x=pred,y=fit),colour="black",size=2)+
  scale_x_continuous(limits = c(8, 36),breaks=(c(10,15,20,25,30,35)))+
  scale_colour_manual(values=c("black","black",rep("black",length(unique(allfitsacc$varbydataset)))))+
  scale_linetype_manual(values=c("solid"))
ggderivpercentfit<-LNCDR::lunaize(ggderivpercentfit)+xlab("Age (years)")+ylab("% of Total")+theme(legend.position = "none",legend.title = element_blank())+
  theme(strip.text.y= element_text(angle = 360))

ggsave(ggderivpercentfit,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/MatRasters/percentchangesigacc.fit.plot.pdf",height=4,width=6)


allfitslat<-allfits[allfits$outcome %in% unique(allinterpolatedfitssig$outcome) & allfits$type=="lat",]

latmetafitsigsmooth<-mgcvscalefits(latmetafitsig,outcomevars = c("percentmax"),predvars = "age",mformula = as.formula("outcome~s(pred)"),scale=FALSE,interval_inc = .1) ###fit smooth

ggderivpercentfit<-ggplot(allfitslat,aes(x=pred,y=fitscale,colour=varbydataset))+geom_line(size=.5,alpha=.35)+
  theme(legend.position = "none")+
  geom_line(data=latmetafitsigsmooth,aes(x=pred,y=fit),colour="black",size=2)+
  scale_x_continuous(limits = c(8, 36),breaks=(c(10,15,20,25,30,35)))+
  scale_colour_manual(values=c("black","black",rep("black",length(unique(allfitslat$varbydataset)))))+
  scale_linetype_manual(values=c("solid"))
ggderivpercentfit<-LNCDR::lunaize(ggderivpercentfit)+xlab("Age (years)")+ylab("% of Total")+theme(legend.position = "none",legend.title = element_blank())+
  theme(strip.text.y= element_text(angle = 360))

ggsave(ggderivpercentfit,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/MatRasters/percentchangesiglat.fit.plot.pdf",height=4,width=6)


####supporting data save#######
metafit<-accmetafitsigsmooth[,c("pred","fitscale")]
metafit$varbydataset<-"meta acc"
metafit$dataset<-"meta"
allfitsacc_r<-allfitsacc[,c("pred","fitscale","dataset","varbydataset")]

allaccfits<-plyr::rbind.fill(metafit,allfitsacc_r)
allaccfits$type<-"accuracy"


metafitlat<-latmetafitsigsmooth[,c("pred","fitscale")]
metafitlat$varbydataset<-"meta lat"
metafitlat$dataset<-"meta"
allfitslat_r<-allfitslat[,c("pred","fitscale","dataset","varbydataset")]

alllatfits<-plyr::rbind.fill(metafitlat,allfitslat_r)
alllatfits$type<-"latency"

allsave<-rbind(allaccfits,alllatfits)

write.csv(allsave,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Data/SupportingData/Sup5.csv")


