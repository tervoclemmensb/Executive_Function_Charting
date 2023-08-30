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
###Read full data###
###########READ Full DATA#########
#############################
LUNA<-read.csv("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/data/btc_R03cleaneddata_20220306.outlierremoved.compositeacclat.csv")
LUNA$id<-as.factor(LUNA$id)
LUNA$metaage<-LUNA$Ageatvisit
LUNA$dataset<-"LUNA"
####
NCANDA<-read.csv("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/NCANDA/data/btc_NCANDAscoredmeasures_20220306.outlierremoved.compositeacclat.csv")
NCANDA$id<-as.factor(NCANDA$subject)
NCANDA$metaage<-NCANDA$cnp_age
NCANDA$visitnum<-as.numeric(as.factor(NCANDA$visit))
NCANDA$dataset<-"NCANDA"
####
NKI<-read.csv("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/NKI/Data/btc_NKIscoredmeasures_20220306.outlierremoved.compositeacclat.csv")
NKI$id<-as.factor(NKI$subject)
NKI$metaage<-NKI$age
NKI$dataset<-"NKI"
NKI$visitnum<-1 ##cross-sectional
####
PNC<-read.csv("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/PNC/Data/btc_PNCscoredmeasures_20220306.outlierremoved.compositeacclat.csv")
PNC$id<-as.factor(PNC$dbGaP_Subject_ID)
PNC$metaage<-PNC$age
PNC$dataset<-"PNC"
PNC$visitnum<-1 ##cross-sectional
###fixed age grid###
agegridforinterp<-data.frame(ages=as.numeric(seq(8,35,by=.1)))
####Read common adult scale fits######## 
load("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/AdultScaledFits/LUNAallfitsadultscaled.Rdata")
LUNAscaledfits<-scalefitsall2030sig
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


####all data######
commonvars<-c("id","metaage","visitnum","dataset","Accuracycomposite","Latencycomposite")

Allorigdata<-plyr::rbind.fill(LUNA[,c(commonvars,unlist(strsplit(unique(LUNAscaledfits$outcome),"_20_30")))],NCANDA[,c(commonvars,unlist(strsplit(unique(NCANDAscaledfits$outcome),"_20_30")))]) %>% 
  plyr::rbind.fill(.,NKI[,c(commonvars,unlist(strsplit(unique(NKIscaledfits$outcome),"_20_30")))]) %>% plyr::rbind.fill(.,PNC[,c(commonvars,unlist(strsplit(unique(PNCscaledfits$outcome),"_20_30")))])

Allorigdata$id<-as.factor(Allorigdata$id)

##########
###leave one dataset out for basis function regression####
agegridforbasisregression<-agegridforinterp

outofsamplebasismodelcomp<-lapply(unique(allinterpolatedfits$dataset),function(di){
  print(di)
  
  fitswithoutdiacc<-allinterpolatedfits[!(allinterpolatedfits$dataset %in% di) & (allinterpolatedfits$type=="acc"),] ###select fits that DO NOT include current dataset (leave one out)
  fitswithoutdilat<-allinterpolatedfits[!(allinterpolatedfits$dataset %in% di) & (allinterpolatedfits$type=="lat"),] ###select fits that DO NOT include current dataset (leave one out)
  
  accvarsfortest<-unique(allinterpolatedfits$outcome[allinterpolatedfits$dataset==di & (allinterpolatedfits$type=="acc")]) ###indentify vars to test basis on in left out dataset
  accvarsfortest<-as.character(strsplit(accvarsfortest,"_20_30"))
  latvarsfortest<-unique(allinterpolatedfits$outcome[allinterpolatedfits$dataset==di & (allinterpolatedfits$type=="lat")]) ###indentify vars to test basis on in left out dataset
  latvarsfortest<-as.character(strsplit(latvarsfortest,"_20_30"))
  
  didata<-Allorigdata[Allorigdata$dataset==di,] ### select original data for validation#####
  
  accmetafitwithoutdi<-metabyagethreelevel(interpolatedagedfs=fitswithoutdiacc,agevar="ages",valcol="fit_value",secol ="se_value",variablenestvar="varglobal",datasetvar = "dataset") ###three level meta to define basis function
  accmetafitwithoutdismooth<-mgcvscalefits(accmetafitwithoutdi,outcomevars = c("estimate"),predvars = "age",mformula = as.formula("outcome~s(pred)"),scale=FALSE,interval_inc = .1) ###fit smooth
  
  didata$basisageacc<-approx(x=round(accmetafitwithoutdismooth$pred,1), y = accmetafitwithoutdismooth$fit, xout=round(didata$metaage,1),method="linear")[[2]] ### acc basis function for data set i 
  #print(hist(didata$metaage[is.na(didata$basisageacc)])) ###print missing ages to ensure only where basis function didn't cover (outside of age range)
  
  latmetafitwithoutdi<-metabyagethreelevel(interpolatedagedfs=fitswithoutdilat,agevar="ages",valcol="fit_value",secol ="se_value",variablenestvar="varglobal",datasetvar = "dataset") ###three level meta to define basis function
  latmetafitwithoutdismooth<-mgcvscalefits(latmetafitwithoutdi,outcomevars = c("estimate"),predvars = "age",mformula = as.formula("outcome~s(pred)"),scale=FALSE,interval_inc = .1) ###fit smooth
  
  didata$basisagelat<-approx(x=round(latmetafitwithoutdismooth$pred,1), y = latmetafitwithoutdismooth$fit, xout=round(didata$metaage,1),method="linear")[[2]] ### acc basis function for data set i 
  
  didata<-didata[didata$metaage > min(c(latmetafitwithoutdismooth$pred,accmetafitwithoutdismooth$pred)) & didata$metaage < max(c(latmetafitwithoutdismooth$pred,accmetafitwithoutdismooth$pred)),]
  
  #print(hist(didata$metaage[is.na(didata$basisagelat)])) ###print missing ages to ensure only where basis function didn't cover (outside of age range)
  
  ####standard basis functions used in dev research
  didata$linearage<-didata$metaage 
  didata$invage<-1/didata$metaage
  didata$agesq<-didata$metaage^2
  
  ##plot to monitor basis function change across left out datasets
  plot(didata$metaage,didata$basisageacc)
  plot(didata$metaage,didata$basisagelat)
  
  didata[,c(accvarsfortest,latvarsfortest,"basisageacc","basisagelat","linearage","invage","agesq")]<-lapply(didata[,c(accvarsfortest,latvarsfortest,"basisageacc","basisagelat","linearage","invage","agesq")],function(x){scale(x)}) ###z scale all di vars to help with model convergence
  require(performance)
  require(see)
  #https://easystats.github.io/performance/reference/compare_performance.html
  
  ####plots for visualizing workflow#########
  
  # accviz<-ggplot()+
  #   geom_line(data=accmetafitwithoutdismooth,aes(x=pred,y=fit),linetype="solid",colour="#c9270e",size=1.5)+
  #   geom_line(data=accmetafitwithoutdi,aes(x=age,y=estimate),colour="black",size=1)
  # accviz<-LNCDR::lunaize(accviz)+xlab("Age(years)")+theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())
  # latviz<-ggplot()+
  #   geom_line(data=latmetafitwithoutdismooth,aes(x=pred,y=fit),linetype="solid",colour="#2f42bd",size=1.5)+
  #   geom_line(data=latmetafitwithoutdi,aes(x=age,y=estimate),colour="black",size=1)
  # latviz<-LNCDR::lunaize(latviz)+xlab("Age(years)")+theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())

  if(di %in% c("LUNA","NCANDA")){
    print("longitudinal")
    allbasismodsout<-lapply(c(accvarsfortest,latvarsfortest),function(vari){
      print(vari)
      didata$outcome<-didata[,vari]
      if(vari %in% accvarsfortest){
        basistype<-"acc"
        didata$basisfunctiontest<-didata$basisageacc
        "acc basis"
      }else if (vari %in% latvarsfortest){
        basistype<-"lat"
        didata$basisfunctiontest<-didata$basisagelat
        "lat basis"
      }
      dimodelbasis<-lme4::lmer(outcome~basisfunctiontest+visitnum+(1+basisfunctiontest|id),control=lmerControl(optCtrl=list(maxfun=10000)),data=didata)
      dimodellin<-lme4::lmer(outcome~linearage+visitnum+(1+linearage|id),control=lmerControl(optCtrl=list(maxfun=10000)),data=didata)
      dimodelinv<-lme4::lmer(outcome~invage+visitnum+(1+invage|id),control=lmerControl(optCtrl=list(maxfun=10000)),data=didata)
      dimodelagesq<-lme4::lmer(outcome~linearage+agesq+visitnum+(1+linearage|id),control=lmerControl(optCtrl=list(maxfun=10000)),data=didata)
      dimodelintonly<-lme4::lmer(outcome~visitnum+(1|id),control=lmerControl(optCtrl=list(maxfun=10000)),data=didata)
      modcompare<-performance::compare_performance(dimodelbasis, dimodellin, dimodelinv,dimodelagesq,dimodelintonly,rank=TRUE)
      modcompare$Performance_Scorerank<-rev(rank(modcompare$Performance_Score)) ##higher score better
      modcompare$RMSErank<-rank(modcompare$RMSE)
      modcompareAIC<-performance::compare_performance(dimodelbasis, dimodellin, dimodelinv,dimodelagesq,dimodelintonly,metrics="AIC")
      modcompareAIC$AICrank<-rank(modcompareAIC$AIC)
      allmodcompare<-merge(modcompare,modcompareAIC,by="Name")
      ###exclude age squared so can limit when not significant#######
      modcomparenoagesq<-performance::compare_performance(dimodelbasis, dimodellin, dimodelinv,dimodelintonly,rank=TRUE)
      modcomparenoagesq$Performance_Scorerank<-rev(rank(modcomparenoagesq$Performance_Score)) ##higher score better
      modcomparenoagesq$RMSErank<-rank(modcomparenoagesq$RMSE)
      modcomparenoagesqAIC<-performance::compare_performance(dimodelbasis, dimodellin, dimodelinv,dimodelintonly,metrics="AIC")
      modcomparenoagesqAIC$AICrank<-rank(modcomparenoagesqAIC$AIC)
      allmodcomparenoagesq<-merge(modcomparenoagesq,modcomparenoagesqAIC,by="Name")
      names(allmodcomparenoagesq)[names(allmodcomparenoagesq)!="Name"]<-paste0("noagesq_",names(allmodcomparenoagesq)[names(allmodcomparenoagesq)!="Name"])
      
      allmodcomparewithnoagesq<-merge(allmodcompare,allmodcomparenoagesq,by="Name",all.x=TRUE)
      
      modcompareout<-data.frame(name=allmodcomparewithnoagesq$Name,Performancescore=allmodcomparewithnoagesq$Performance_Score,
                                Performance_Scorerank=allmodcomparewithnoagesq$Performance_Scorerank,AICrank=allmodcomparewithnoagesq$AICrank,
                                RMSErank=allmodcomparewithnoagesq$RMSErank,Performancescorenorsq=allmodcomparewithnoagesq$noagesq_Performance_Score,
                                Performancescore_ranknorsq=allmodcomparewithnoagesq$noagesq_Performance_Scorerank,AICranknorsq=allmodcomparewithnoagesq$noagesq_AICrank,
                                RMSEranknorsq=allmodcomparewithnoagesq$noagesq_RMSErank)
      
      modcompareout$quadtermp<-car::Anova(dimodelagesq)$`Pr(>Chisq)`[row.names(car::Anova(dimodelagesq))=="agesq"]
      
      modcompareout$var<-vari
      modcompareout$type<-basistype
      return(modcompareout)
    })
  }
  if(!(di %in% c("LUNA","NCANDA"))){
    print("cross-sectional")
    allbasismodsout<-lapply(c(accvarsfortest,latvarsfortest),function(vari){
      print(vari)
      didata$outcome<-didata[,vari]
      if(vari %in% accvarsfortest){
        basistype<-"acc"
        didata$basisfunctiontest<-didata$basisageacc
        "acc basis"
      }else if (vari %in% latvarsfortest){
        basistype<-"lat"
        didata$basisfunctiontest<-didata$basisagelat
        "lat basis"
      }
    dimodelbasis<-lm(outcome~basisageacc,data=didata)
    dimodellin<-lm(outcome~linearage,data=didata)
    dimodelinv<-lm(outcome~invage,data=didata)
    dimodelagesq<-lm(outcome~linearage+agesq,data=didata)
    dimodelintonly<-lm(outcome~1,data=didata)
    modcompare<-performance::compare_performance(dimodelbasis, dimodellin, dimodelinv,dimodelagesq,dimodelintonly,rank=TRUE)
    modcompare$Performance_Scorerank<-rev(rank(modcompare$Performance_Score)) ##higher score better
    modcompare$RMSErank<-rank(modcompare$RMSE)
    modcompareAIC<-performance::compare_performance(dimodelbasis, dimodellin, dimodelinv,dimodelagesq,dimodelintonly,metrics="AIC")
    modcompareAIC$AICrank<-rank(modcompareAIC$AIC)
    allmodcompare<-merge(modcompare,modcompareAIC,by="Name")
    ###exclude age squared so can limit when not significant#######
    modcomparenoagesq<-performance::compare_performance(dimodelbasis, dimodellin, dimodelinv,dimodelintonly,rank=TRUE)
    modcomparenoagesq$Performance_Scorerank<-rev(rank(modcomparenoagesq$Performance_Score)) ##higher score better
    modcomparenoagesq$RMSErank<-rank(modcomparenoagesq$RMSE)
    modcomparenoagesqAIC<-performance::compare_performance(dimodelbasis, dimodellin, dimodelinv,dimodelintonly,metrics="AIC")
    modcomparenoagesqAIC$AICrank<-rank(modcomparenoagesqAIC$AIC)
    allmodcomparenoagesq<-merge(modcomparenoagesq,modcomparenoagesqAIC,by="Name")
    names(allmodcomparenoagesq)[names(allmodcomparenoagesq)!="Name"]<-paste0("noagesq_",names(allmodcomparenoagesq)[names(allmodcomparenoagesq)!="Name"])
    
    allmodcomparewithnoagesq<-merge(allmodcompare,allmodcomparenoagesq,by="Name",all.x=TRUE)
    
    modcompareout<-data.frame(name=allmodcomparewithnoagesq$Name,Performancescore=allmodcomparewithnoagesq$Performance_Score,
                              Performance_Scorerank=allmodcomparewithnoagesq$Performance_Scorerank,AICrank=allmodcomparewithnoagesq$AICrank,
                              RMSErank=allmodcomparewithnoagesq$RMSErank,Performancescorenorsq=allmodcomparewithnoagesq$noagesq_Performance_Score,
                              Performancescore_ranknorsq=allmodcomparewithnoagesq$noagesq_Performance_Scorerank,AICranknorsq=allmodcomparewithnoagesq$noagesq_AICrank,
                              RMSEranknorsq=allmodcomparewithnoagesq$noagesq_RMSErank)
    
    modcompareout$quadtermp<-summary(dimodelagesq)$coefficients[row.names(summary(dimodelagesq)$coefficients)=="agesq",4]
    modcompareout$var<-vari
    modcompareout$type<-basistype
    return(modcompareout)
    })
  }
  allbasismodsoutdf<-do.call(rbind,allbasismodsout)
  allbasismodsoutdf$dataset<-di
  return(allbasismodsoutdf)
})

outofsamplebasismodelcompdf<-do.call(rbind,outofsamplebasismodelcomp)
outofsamplebasismodelcompdf$simpname<-outofsamplebasismodelcompdf$name
outofsamplebasismodelcompdf$simpname<-unlist(lapply(strsplit(outofsamplebasismodelcompdf$simpname,"dimodel"),"[[",2))
outofsamplebasismodelcompdf$simpnamef<-factor(outofsamplebasismodelcompdf$simpname,levels=c("basis","agesq","inv","lin","intonly"))
outofsamplebasismodelcompdf$simpnamefbytype<-paste(outofsamplebasismodelcompdf$simpnamef,outofsamplebasismodelcompdf$type,sep="_")
outofsamplebasismodelcompdf$simpnamefbytypef<-factor(outofsamplebasismodelcompdf$simpnamefbytype,
                                                     c(paste(levels(outofsamplebasismodelcompdf$simpnamef),c("acc"),sep="_"),
                                                       paste(levels(outofsamplebasismodelcompdf$simpnamef),c("lat"),sep="_")))
outofsamplebasismodelcompdf$typefull<-dplyr::if_else(outofsamplebasismodelcompdf$type=="acc","Accuracy","Latency")
modelcomparstructformerge<-outofsamplebasismodelcompdf %>% group_by(simpnamefbytypef,simpnamef,typefull) %>% summarize()
###Figure 5C#################

blues<-brewer.pal(n=9,"Blues")
blues<-c(blues,"#2f42bd")
blues<-blues[c(10,9,8,5,3)]
reds<-brewer.pal(n=9,"Reds")
reds<-c(reds,"#c9270e")
reds<-reds[c(10,9,8,6,3)]
mycolors<-c(reds,blues)
ggbasisallvarsmodcompare<-ggplot(outofsamplebasismodelcompdf,aes(x=simpnamef,y=Performancescore,colour=simpnamefbytypef,fill=simpnamefbytypef))+geom_boxplot(alpha=.65,outlier.shape=NA,colour="black")+geom_jitter(height=0,alpha=.85)+facet_grid(~typefull,scales = "free")+
  scale_colour_manual(values=mycolors)+scale_fill_manual(values=mycolors)
ggbasisallvarsmodcompare<-LNCDR::lunaize(ggbasisallvarsmodcompare)+ylab("Age Model Performance Score\n (Percentile Among Models)")+xlab("")+theme(legend.position = "none",strip.background = element_blank(),strip.text.x = (element_text(size = 20)))
save(ggbasisallvarsmodcompare,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/Metabasis/allmeasuresboxplot.plot.Rdata")

###acconly#####
ggbasisallvarsmodcompareacc<-ggplot(outofsamplebasismodelcompdf[outofsamplebasismodelcompdf$type=="acc",],aes(x=simpnamef,y=Performancescore,colour=simpnamefbytypef,fill=simpnamefbytypef))+geom_boxplot(aes(alpha=simpnamefbytypef),outlier.shape=NA,colour="black")+
  geom_jitter(aes(fill=simpnamefbytypef),height=0,colour="black",shape=21,width=.225)+
  scale_fill_manual(values=reds)+scale_alpha_manual(values=c(.9,.75,.6,.5,.3))
ggbasisallvarsmodcompareacc<-LNCDR::lunaize(ggbasisallvarsmodcompareacc)+ylab("Age Model Performance Score\n (Percentile Among Models)")+xlab("")+theme(legend.position = "none",strip.background = element_blank(),strip.text.x = (element_text(size = 20)))
save(ggbasisallvarsmodcompareacc,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/Metabasis/allmeasuresaccboxplot.plot.Rdata")

##latonly#####
ggbasisallvarsmodcomparelat<-ggplot(outofsamplebasismodelcompdf[outofsamplebasismodelcompdf$type=="lat",],aes(x=simpnamef,y=Performancescore,colour=simpnamefbytypef,fill=simpnamefbytypef))+geom_boxplot(aes(alpha=simpnamefbytypef),outlier.shape=NA,colour="black")+
  geom_jitter(aes(fill=simpnamefbytypef),height=0,colour="black",shape=21,width=.225)+
  scale_fill_manual(values=blues)+scale_alpha_manual(values=c(.9,.75,.6,.5,.3))
ggbasisallvarsmodcomparelat<-LNCDR::lunaize(ggbasisallvarsmodcomparelat)+ylab("Age Model Performance Score\n (Percentile Among Models)")+xlab("")+theme(legend.position = "none",strip.background = element_blank(),strip.text.x = (element_text(size = 20)))
save(ggbasisallvarsmodcomparelat,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/Metabasis/allmeasureslatboxplot.plot.Rdata")


###supportdata###

outofsamplebasismodelcompdf_r<-outofsamplebasismodelcompdf[,c("type","simpnamef","Performancescore","simpnamefbytypef")]
write.csv(outofsamplebasismodelcompdf_r,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Data/SupportingData/Figure5B.csv")
#######

####set up CP for pie charts
cp <- coord_polar(theta = "y")
cp$is_free <- function() TRUE
###PerfomanceRANK
outofsamplebasismodelcompdfsummary_PR<-outofsamplebasismodelcompdf %>% count(simpnamefbytypef,Performance_Scorerank) %>% 
  group_by(simpnamefbytypef) %>% mutate(prop=prop.table(n))

outofsamplebasismodelcompdfsummary_PR1<-outofsamplebasismodelcompdfsummary_PR[outofsamplebasismodelcompdfsummary_PR$Performance_Scorerank==1,]

outofsamplebasismodelcompdfsummary_PR1<-merge(outofsamplebasismodelcompdfsummary_PR1,modelcomparstructformerge,by=c("simpnamefbytypef"))
outofsamplebasismodelcompdfsummary_PR1$proplabel<-paste0(round(outofsamplebasismodelcompdfsummary_PR1$prop,2)*100,"%")
mycolourspie<-mycolors[levels(outofsamplebasismodelcompdfsummary_PR1$simpnamefbytypef) %in% unique(outofsamplebasismodelcompdfsummary_PR1$simpnamefbytypef)]
myalphaspie<-c(.9,.75,.6,.9,.75,.6)
ggpiePR<-ggplot(outofsamplebasismodelcompdfsummary_PR1,aes(x="", y=prop, fill=as.factor(simpnamefbytypef))) +
  geom_bar(aes(fill=as.factor(simpnamefbytypef),alpha=as.factor(simpnamefbytypef)),stat="identity", width=1, color="white")+cp+facet_wrap(~typefull,scales="free")+theme(aspect.ratio=1)+theme_void()+
  scale_colour_manual(values=mycolourspie)+scale_fill_manual(values=mycolourspie)+scale_alpha_manual(values=myalphaspie)+
  geom_text(aes(label = proplabel),size=6,
            position = position_stack(vjust = 0.5))+theme_void()

ggpiePR<-LNCDR::lunaize(ggpiePR)+theme_void()+theme(strip.text.x = (element_text(size = 20)))
ggsave(ggpiePR,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/Metabasis/allmeasuresPRrankpie.plot.pdf",height=4,width=12)
save(ggpiePR,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/Metabasis/allmeasuresPRrankpie.plot.Rdata")

outofsamplebasismodelcompdfsummary_PR1<-outofsamplebasismodelcompdfsummary_PR1 %>% group_by(typefull) %>% mutate(totalbytype=sum(n))
outofsamplebasismodelcompdfsummary_PR1$nnotop<-outofsamplebasismodelcompdfsummary_PR1$totalbytype-outofsamplebasismodelcompdfsummary_PR1$n
###overall tests (across accuracy and latency)####
overallsumtest<-outofsamplebasismodelcompdfsummary_PR1 %>% dplyr::group_by(simpnamef) %>% dplyr::summarize(mean=mean(prop),n=sum(n),nnotop=sum(nnotop),totalbytype=sum(totalbytype))
compares<-expand.grid(basis="basis",compare=c("agesq","inv","lin","intonly"))
chisqsoutoverall<-lapply(1:nrow(compares),function(cri){
  cridf<-overallsumtest[overallsumtest$simpnamef %in% c(compares$basis[cri],compares$compare[cri]),]
  if(!(compares$compare[cri] %in% cridf$simpnamef)){
    cridf<-plyr::rbind.fill(cridf,data.frame(simpnamef=compares$compare[cri],n=0,nnotop=cridf$totalbytype[1]))## if no top ranked models for compare manually specify 0
  }
  print(cridf)
  criout<-broom::tidy(chisq.test(cridf[,c("n","nnotop")]))
  criout$basis=compares$basis[cri]
  criout$compare=compares$compare[cri]
  return(criout)
})

chisqsoutoveralldf<-do.call(rbind,chisqsoutoverall)

###test by accuracy and latency separetly
comparesacc<-expand.grid(basis="basis_acc",compare=c("agesq_acc","inv_acc","lin_acc","intonly_acc"))
compareslat<-expand.grid(basis="basis_lat",compare=c("agesq_lat","inv_lat","lin_lat","intonly_lat"))
allcompares<-rbind(comparesacc,compareslat)
chisqsout<-lapply(1:nrow(allcompares),function(cri){
  cridf<-outofsamplebasismodelcompdfsummary_PR1[outofsamplebasismodelcompdfsummary_PR1$simpnamefbytypef %in% c(allcompares$basis[cri],allcompares$compare[cri]),]
  if(!(allcompares$compare[cri] %in% cridf$simpnamefbytypef)){
    cridf<-plyr::rbind.fill(cridf,data.frame(simpnamefbytypef=allcompares$compare[cri],n=0,nnotop=cridf$totalbytype[1]))## if no top ranked models for compare manually specify 0
  }
  print(cridf)
  criout<-broom::tidy(chisq.test(cridf[,c("n","nnotop")]))
  criout$basis=allcompares$basis[cri]
  criout$compare=allcompares$compare[cri]
  return(criout)
})

chisqsoutdf<-do.call(rbind,chisqsout)


####ggPRAcc

ggpiePRacc<-ggplot(outofsamplebasismodelcompdfsummary_PR1[outofsamplebasismodelcompdfsummary_PR1$typefull=="Accuracy",],aes(x="", y=prop, fill=as.factor(simpnamefbytypef))) +
  geom_bar(aes(fill=as.factor(simpnamefbytypef),alpha=as.factor(simpnamefbytypef)),stat="identity", width=1, color="white")+cp+theme(aspect.ratio=1)+theme_void()+
  scale_colour_manual(values=mycolourspie)+scale_fill_manual(values=mycolourspie)+scale_alpha_manual(values=myalphaspie)+
  geom_text(aes(label = proplabel),size=6,
            position = position_stack(vjust = 0.5))+theme_void()
ggpiePRacc<-LNCDR::lunaize(ggpiePRacc)+theme_void()+theme(legend.position = "none")
ggsave(ggpiePRacc,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/Metabasis/allmeasuresPRrankpie.accplot.pdf",height=4,width=4)

###ggPRlat
ggpiePRlat<-ggplot(outofsamplebasismodelcompdfsummary_PR1[outofsamplebasismodelcompdfsummary_PR1$typefull=="Latency",],aes(x="", y=prop, fill=as.factor(simpnamefbytypef))) +
  geom_bar(aes(fill=as.factor(simpnamefbytypef),alpha=as.factor(simpnamefbytypef)),stat="identity", width=1, color="white")+cp+theme(aspect.ratio=1)+theme_void()+
  scale_colour_manual(values=mycolourspie[4:6])+scale_fill_manual(values=mycolourspie[4:6])+scale_alpha_manual(values=myalphaspie[4:6])+
  geom_text(aes(label = proplabel),size=6,
            position = position_stack(vjust = 0.5))+theme_void()
ggpiePRlat<-LNCDR::lunaize(ggpiePRlat)+theme_void()+theme(legend.position = "none")
ggsave(ggpiePRlat,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/Metabasis/allmeasuresPRrankpie.latencyplot.pdf",height=4,width=4)

#####Sensitivity just LUNA and NKI (no overlapping measures between them)

allinterpolatedfits_lunaNKI<-allinterpolatedfits[allinterpolatedfits$dataset %in% c("LUNA","NKI"),]

outofsamplebasismodelcomp_sens<-lapply(c("LUNA","NKI"),function(di){
  print(di)
  
  fitswithoutdiacc<-allinterpolatedfits_lunaNKI[!(allinterpolatedfits_lunaNKI$dataset %in% di) & (allinterpolatedfits_lunaNKI$type=="acc"),] ###select fits that DO NOT include current dataset (leave one out)
  fitswithoutdilat<-allinterpolatedfits_lunaNKI[!(allinterpolatedfits_lunaNKI$dataset %in% di) & (allinterpolatedfits_lunaNKI$type=="lat"),] ###select fits that DO NOT include current dataset (leave one out)
  
  accvarsfortest<-unique(allinterpolatedfits$outcome[allinterpolatedfits$dataset==di & (allinterpolatedfits$type=="acc")]) ###indentify vars to test basis on in left out dataset
  accvarsfortest<-as.character(strsplit(accvarsfortest,"_20_30"))
  latvarsfortest<-unique(allinterpolatedfits$outcome[allinterpolatedfits$dataset==di & (allinterpolatedfits$type=="lat")]) ###indentify vars to test basis on in left out dataset
  latvarsfortest<-as.character(strsplit(latvarsfortest,"_20_30"))
  
  didata<-Allorigdata[Allorigdata$dataset==di,] ### select original data for validation#####
  
  accmetafitwithoutdi<-fitswithoutdiacc %>% dplyr::group_by(ages) %>% dplyr::summarize(estimate=mean(fit_value)) ###only one dataset here so just average
  accmetafitwithoutdismooth<-mgcvscalefits(accmetafitwithoutdi,outcomevars = c("estimate"),predvars = "ages",mformula = as.formula("outcome~s(pred)"),scale=FALSE,interval_inc = .1) ###fit smooth
  
  didata$basisageacc<-approx(x=round(accmetafitwithoutdismooth$pred,1), y = accmetafitwithoutdismooth$fit, xout=round(didata$metaage,1),method="linear")[[2]] ### acc basis function for data set i 
  #print(hist(didata$metaage[is.na(didata$basisageacc)])) ###print missing ages to ensure only where basis function didn't cover (outside of age range)
  
  latmetafitwithoutdi<-fitswithoutdilat %>% dplyr::group_by(ages) %>% dplyr::summarize(estimate=mean(fit_value)) ###only one dataset here so just average
  latmetafitwithoutdismooth<-mgcvscalefits(latmetafitwithoutdi,outcomevars = c("estimate"),predvars = "ages",mformula = as.formula("outcome~s(pred)"),scale=FALSE,interval_inc = .1) ###fit smooth
  
  didata$basisagelat<-approx(x=round(latmetafitwithoutdismooth$pred,1), y = latmetafitwithoutdismooth$fit, xout=round(didata$metaage,1),method="linear")[[2]] ### acc basis function for data set i 
  
  didata<-didata[didata$metaage > min(c(latmetafitwithoutdismooth$pred,accmetafitwithoutdismooth$pred)) & didata$metaage < max(c(latmetafitwithoutdismooth$pred,accmetafitwithoutdismooth$pred)),]
  
  #print(hist(didata$metaage[is.na(didata$basisagelat)])) ###print missing ages to ensure only where basis function didn't cover (outside of age range)
  
  ####standard basis functions used in dev research
  didata$linearage<-didata$metaage 
  didata$invage<-1/didata$metaage
  didata$agesq<-didata$metaage^2
  
  ##plot to monitor basis function change across left out datasets
  plot(didata$metaage,didata$basisageacc)
  plot(didata$metaage,didata$basisagelat)
  
  didata[,c(accvarsfortest,latvarsfortest,"basisageacc","basisagelat","linearage","invage","agesq")]<-lapply(didata[,c(accvarsfortest,latvarsfortest,"basisageacc","basisagelat","linearage","invage","agesq")],function(x){scale(x)}) ###z scale all di vars to help with model convergence
  require(performance)
  require(see)
  #https://easystats.github.io/performance/reference/compare_performance.html
  
  ####plots for visualizing workflow#########
  
  # accviz<-ggplot()+
  #   geom_line(data=accmetafitwithoutdismooth,aes(x=pred,y=fit),linetype="solid",colour="#c9270e",size=1.5)+
  #   geom_line(data=accmetafitwithoutdi,aes(x=age,y=estimate),colour="black",size=1)
  # accviz<-LNCDR::lunaize(accviz)+xlab("Age(years)")+theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())
  # latviz<-ggplot()+
  #   geom_line(data=latmetafitwithoutdismooth,aes(x=pred,y=fit),linetype="solid",colour="#2f42bd",size=1.5)+
  #   geom_line(data=latmetafitwithoutdi,aes(x=age,y=estimate),colour="black",size=1)
  # latviz<-LNCDR::lunaize(latviz)+xlab("Age(years)")+theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())
  # 
  if(di %in% c("LUNA","NCANDA")){
    print("longitudinal")
    allbasismodsout<-lapply(c(accvarsfortest,latvarsfortest),function(vari){
      print(vari)
      didata$outcome<-didata[,vari]
      if(vari %in% accvarsfortest){
        basistype<-"acc"
        didata$basisfunctiontest<-didata$basisageacc
        "acc basis"
      }else if (vari %in% latvarsfortest){
        basistype<-"lat"
        didata$basisfunctiontest<-didata$basisagelat
        "lat basis"
      }
      dimodelbasis<-lme4::lmer(outcome~basisfunctiontest+visitnum+(1+basisfunctiontest|id),control=lmerControl(optCtrl=list(maxfun=10000)),data=didata)
      dimodellin<-lme4::lmer(outcome~linearage+visitnum+(1+linearage|id),control=lmerControl(optCtrl=list(maxfun=10000)),data=didata)
      dimodelinv<-lme4::lmer(outcome~invage+visitnum+(1+invage|id),control=lmerControl(optCtrl=list(maxfun=10000)),data=didata)
      dimodelagesq<-lme4::lmer(outcome~linearage+agesq+visitnum+(1+linearage|id),control=lmerControl(optCtrl=list(maxfun=10000)),data=didata)
      dimodelintonly<-lme4::lmer(outcome~visitnum+(1|id),control=lmerControl(optCtrl=list(maxfun=10000)),data=didata)
      modcompare<-performance::compare_performance(dimodelbasis, dimodellin, dimodelinv,dimodelagesq,dimodelintonly,rank=TRUE)
      modcompare$Performance_Scorerank<-rev(rank(modcompare$Performance_Score)) ##higher score better
      modcompare$RMSErank<-rank(modcompare$RMSE)
      modcompareAIC<-performance::compare_performance(dimodelbasis, dimodellin, dimodelinv,dimodelagesq,dimodelintonly,metrics="AIC")
      modcompareAIC$AICrank<-rank(modcompareAIC$AIC)
      allmodcompare<-merge(modcompare,modcompareAIC,by="Name")
      ###exclude age squared so can limit when not significant#######
      modcomparenoagesq<-performance::compare_performance(dimodelbasis, dimodellin, dimodelinv,dimodelintonly,rank=TRUE)
      modcomparenoagesq$Performance_Scorerank<-rev(rank(modcomparenoagesq$Performance_Score)) ##higher score better
      modcomparenoagesq$RMSErank<-rank(modcomparenoagesq$RMSE)
      modcomparenoagesqAIC<-performance::compare_performance(dimodelbasis, dimodellin, dimodelinv,dimodelintonly,metrics="AIC")
      modcomparenoagesqAIC$AICrank<-rank(modcomparenoagesqAIC$AIC)
      allmodcomparenoagesq<-merge(modcomparenoagesq,modcomparenoagesqAIC,by="Name")
      names(allmodcomparenoagesq)[names(allmodcomparenoagesq)!="Name"]<-paste0("noagesq_",names(allmodcomparenoagesq)[names(allmodcomparenoagesq)!="Name"])
      
      allmodcomparewithnoagesq<-merge(allmodcompare,allmodcomparenoagesq,by="Name",all.x=TRUE)
      
      modcompareout<-data.frame(name=allmodcomparewithnoagesq$Name,Performancescore=allmodcomparewithnoagesq$Performance_Score,
                                Performance_Scorerank=allmodcomparewithnoagesq$Performance_Scorerank,AICrank=allmodcomparewithnoagesq$AICrank,
                                RMSErank=allmodcomparewithnoagesq$RMSErank,Performancescorenorsq=allmodcomparewithnoagesq$noagesq_Performance_Score,
                                Performancescore_ranknorsq=allmodcomparewithnoagesq$noagesq_Performance_Scorerank,AICranknorsq=allmodcomparewithnoagesq$noagesq_AICrank,
                                RMSEranknorsq=allmodcomparewithnoagesq$noagesq_RMSErank)
      
      modcompareout$quadtermp<-car::Anova(dimodelagesq)$`Pr(>Chisq)`[row.names(car::Anova(dimodelagesq))=="agesq"]
      
      modcompareout$var<-vari
      modcompareout$type<-basistype
      return(modcompareout)
    })
  }
  if(!(di %in% c("LUNA","NCANDA"))){
    print("cross-sectional")
    allbasismodsout<-lapply(c(accvarsfortest,latvarsfortest),function(vari){
      print(vari)
      didata$outcome<-didata[,vari]
      if(vari %in% accvarsfortest){
        basistype<-"acc"
        didata$basisfunctiontest<-didata$basisageacc
        "acc basis"
      }else if (vari %in% latvarsfortest){
        basistype<-"lat"
        didata$basisfunctiontest<-didata$basisagelat
        "lat basis"
      }
      dimodelbasis<-lm(outcome~basisageacc,data=didata)
      dimodellin<-lm(outcome~linearage,data=didata)
      dimodelinv<-lm(outcome~invage,data=didata)
      dimodelagesq<-lm(outcome~linearage+agesq,data=didata)
      dimodelintonly<-lm(outcome~1,data=didata)
      modcompare<-performance::compare_performance(dimodelbasis, dimodellin, dimodelinv,dimodelagesq,dimodelintonly,rank=TRUE)
      modcompare$Performance_Scorerank<-rev(rank(modcompare$Performance_Score)) ##higher score better
      modcompare$RMSErank<-rank(modcompare$RMSE)
      modcompareAIC<-performance::compare_performance(dimodelbasis, dimodellin, dimodelinv,dimodelagesq,dimodelintonly,metrics="AIC")
      modcompareAIC$AICrank<-rank(modcompareAIC$AIC)
      allmodcompare<-merge(modcompare,modcompareAIC,by="Name")
      ###exclude age squared so can limit when not significant#######
      modcomparenoagesq<-performance::compare_performance(dimodelbasis, dimodellin, dimodelinv,dimodelintonly,rank=TRUE)
      modcomparenoagesq$Performance_Scorerank<-rev(rank(modcomparenoagesq$Performance_Score)) ##higher score better
      modcomparenoagesq$RMSErank<-rank(modcomparenoagesq$RMSE)
      modcomparenoagesqAIC<-performance::compare_performance(dimodelbasis, dimodellin, dimodelinv,dimodelintonly,metrics="AIC")
      modcomparenoagesqAIC$AICrank<-rank(modcomparenoagesqAIC$AIC)
      allmodcomparenoagesq<-merge(modcomparenoagesq,modcomparenoagesqAIC,by="Name")
      names(allmodcomparenoagesq)[names(allmodcomparenoagesq)!="Name"]<-paste0("noagesq_",names(allmodcomparenoagesq)[names(allmodcomparenoagesq)!="Name"])
      
      allmodcomparewithnoagesq<-merge(allmodcompare,allmodcomparenoagesq,by="Name",all.x=TRUE)
      
      modcompareout<-data.frame(name=allmodcomparewithnoagesq$Name,Performancescore=allmodcomparewithnoagesq$Performance_Score,
                                Performance_Scorerank=allmodcomparewithnoagesq$Performance_Scorerank,AICrank=allmodcomparewithnoagesq$AICrank,
                                RMSErank=allmodcomparewithnoagesq$RMSErank,Performancescorenorsq=allmodcomparewithnoagesq$noagesq_Performance_Score,
                                Performancescore_ranknorsq=allmodcomparewithnoagesq$noagesq_Performance_Scorerank,AICranknorsq=allmodcomparewithnoagesq$noagesq_AICrank,
                                RMSEranknorsq=allmodcomparewithnoagesq$noagesq_RMSErank)
      
      modcompareout$quadtermp<-summary(dimodelagesq)$coefficients[row.names(summary(dimodelagesq)$coefficients)=="agesq",4]
      modcompareout$var<-vari
      modcompareout$type<-basistype
      return(modcompareout)
    })
  }
  allbasismodsoutdf<-do.call(rbind,allbasismodsout)
  allbasismodsoutdf$dataset<-di
  return(allbasismodsoutdf)
})

outofsamplebasismodelcompdf_sens<-do.call(rbind,outofsamplebasismodelcomp_sens)
outofsamplebasismodelcompdf_sens$simpname<-outofsamplebasismodelcompdf_sens$name
outofsamplebasismodelcompdf_sens$simpname<-unlist(lapply(strsplit(outofsamplebasismodelcompdf_sens$simpname,"dimodel"),"[[",2))
outofsamplebasismodelcompdf_sens$simpnamef<-factor(outofsamplebasismodelcompdf_sens$simpname,levels=c("basis","agesq","inv","lin","intonly"))
outofsamplebasismodelcompdf_sens$simpnamefbytype<-paste(outofsamplebasismodelcompdf_sens$simpnamef,outofsamplebasismodelcompdf_sens$type,sep="_")
outofsamplebasismodelcompdf_sens$simpnamefbytypef<-factor(outofsamplebasismodelcompdf_sens$simpnamefbytype,
                                                     c(paste(levels(outofsamplebasismodelcompdf_sens$simpnamef),c("acc"),sep="_"),
                                                       paste(levels(outofsamplebasismodelcompdf_sens$simpnamef),c("lat"),sep="_")))
outofsamplebasismodelcompdf_sens$typefull<-dplyr::if_else(outofsamplebasismodelcompdf_sens$type=="acc","Accuracy","Latency")
modelcomparstructformerge<-outofsamplebasismodelcompdf_sens %>% group_by(simpnamefbytypef,simpnamef,typefull) %>% summarize()
###Figure 5Csens#################

ggbasisallvarsmodcompare_sens<-ggplot(outofsamplebasismodelcompdf_sens,aes(x=simpnamef,y=Performancescore,colour=simpnamefbytypef,fill=simpnamefbytypef))+geom_boxplot(alpha=.65,outlier.shape=NA,colour="black")+geom_jitter(height=0,alpha=.85)+facet_grid(~typefull,scales = "free")+
  scale_colour_manual(values=mycolors)+scale_fill_manual(values=mycolors)
ggbasisallvarsmodcompare_sens<-LNCDR::lunaize(ggbasisallvarsmodcompare_sens)+ylab("Age Model Performance Score\n (Percentile Among Models)")+xlab("")+theme(legend.position = "none",strip.background = element_blank(),strip.text.x = (element_text(size = 20)))

outofsamplebasismodelcompdfsummarysens_PR<-outofsamplebasismodelcompdf_sens %>% count(simpnamefbytypef,Performance_Scorerank) %>% 
  group_by(simpnamefbytypef) %>% mutate(prop=prop.table(n))

outofsamplebasismodelcompdfsummarysens_PR1<-outofsamplebasismodelcompdfsummarysens_PR[outofsamplebasismodelcompdfsummarysens_PR$Performance_Scorerank==1,]

outofsamplebasismodelcompdfsummarysens_PR1 %>% dplyr::summarize(mean=mean(prop))





