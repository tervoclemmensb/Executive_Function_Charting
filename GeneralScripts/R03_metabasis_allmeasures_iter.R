library(dplyr)
library(tidyr)
library(ggplot2)
library(mgcv)
library(parallel)
library(lme4)
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
####offset visualization for workflow#########
allinterpolatedfits_lunaNKI<-allinterpolatedfits[allinterpolatedfits$dataset %in% c("LUNA","NKI"),]

ageoffsets<-c(seq(-5,5,2)*(10),0) ###offsets in years --> convert to 1/10th of years for age grid (multiply by 10)
byageoffsetsLunaNKI<-lapply(ageoffsets,function(thisoffset){
  message_parallel(thisoffset)
  if(sign(thisoffset)==-1){
    leadlagfunction<-function(x,n){
      n<-n*-1
      dplyr::lead(x,n)
    }
  }else{
    leadlagfunction<-function(x,n){
      dplyr::lag(x,n)
    }
  }
  
outofsamplebasismodelcomp_lunanki<-lapply(c("LUNA","NKI"),function(di){
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

  accmetafitwithoutdismoothoffset<-accmetafitwithoutdismooth %>% arrange(pred) %>% mutate(offsetestimate=leadlagfunction(fit,thisoffset))
  
  didata$basisageacc<-approx(x=round(accmetafitwithoutdismoothoffset$pred,1), y = accmetafitwithoutdismoothoffset$offsetestimate, xout=round(didata$metaage,1),method="linear")[[2]] ### acc basis function for data set i 
  
  #print(hist(didata$metaage[is.na(didata$basisageacc)])) ###print missing ages to ensure only where basis function didn't cover (outside of age range)

  latmetafitwithoutdi<-fitswithoutdilat %>% dplyr::group_by(ages) %>% dplyr::summarize(estimate=mean(fit_value)) ###only one dataset here so just average
  latmetafitwithoutdismooth<-mgcvscalefits(latmetafitwithoutdi,outcomevars = c("estimate"),predvars = "ages",mformula = as.formula("outcome~s(pred)"),scale=FALSE,interval_inc = .1) ###fit smooth
  latmetafitwithoutdismoothoffset<-latmetafitwithoutdismooth %>% arrange(pred) %>% mutate(offsetestimate=leadlagfunction(fit,thisoffset))
  
  didata$basisagelat<-approx(x=round(latmetafitwithoutdismoothoffset$pred,1), y = latmetafitwithoutdismoothoffset$offsetestimate, xout=round(didata$metaage,1),method="linear")[[2]] ### acc basis function for data set i 
  
  didata<-didata[didata$metaage > min(c(latmetafitwithoutdismoothoffset$pred[!is.na(latmetafitwithoutdismoothoffset$offsetestimate)],accmetafitwithoutdismoothoffset$pred[!is.na(accmetafitwithoutdismoothoffset$offsetestimate)])) & didata$metaage < max(c(latmetafitwithoutdismoothoffset$pred[!is.na(latmetafitwithoutdismoothoffset$offsetestimate)],accmetafitwithoutdismoothoffset$pred[!is.na(accmetafitwithoutdismoothoffset$offsetestimate)])),]

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
outofsamplebasismodelcompdf<-do.call(rbind,outofsamplebasismodelcomp_lunanki)
outofsamplebasismodelcompdf$ageoffset<-thisoffset
outofsamplebasismodelcompdf$offsetfunction<-deparse1(leadlagfunction)
return(outofsamplebasismodelcompdf)
})

outofsamplebasismodelcompdf<-do.call(rbind,byageoffsetsLunaNKI)
outofsamplebasismodelcompdf$simpname<-outofsamplebasismodelcompdf$name
outofsamplebasismodelcompdf$simpname<-unlist(lapply(strsplit(outofsamplebasismodelcompdf$simpname,"dimodel"),"[[",2))
outofsamplebasismodelcompdf$simpnamef<-factor(outofsamplebasismodelcompdf$simpname,levels=c("basis","agesq","inv","lin","intonly"))
outofsamplebasismodelcompdf$simpnamefbytype<-paste(outofsamplebasismodelcompdf$simpnamef,outofsamplebasismodelcompdf$type,sep="_")
outofsamplebasismodelcompdf$simpnamefbytypef<-factor(outofsamplebasismodelcompdf$simpnamefbytype,
                                                     c(paste(levels(outofsamplebasismodelcompdf$simpnamef),c("acc"),sep="_"),
                                                       paste(levels(outofsamplebasismodelcompdf$simpnamef),c("lat"),sep="_")))
outofsamplebasismodelcompdf$typefull<-dplyr::if_else(outofsamplebasismodelcompdf$type=="acc","Accuracy","Latency")
modelcomparstructformerge<-outofsamplebasismodelcompdf %>% group_by(ageoffset,simpnamefbytypef,simpnamef,typefull) %>% summarize()


ggbasisallvaroffsetcompare<-ggplot(outofsamplebasismodelcompdf[outofsamplebasismodelcompdf$simpnamef=="basis",],aes(x=as.factor(ageoffset/10),y=Performancescore))+geom_boxplot(alpha=.65,outlier.shape=NA,colour="grey55")+geom_jitter(height=0,alpha=.85,colour="grey55")+facet_grid(~typefull,scales = "free")+
  geom_smooth(method = "gam",formula=y~s(x,k=6),se=FALSE, aes(group=1,colour=typefull),size=2)+scale_colour_manual(values=c("#c9270e","#2f42bde6"))
ggbasisallvaroffsetcompare<-LNCDR::lunaize(ggbasisallvaroffsetcompare)+ylab("Age Basis Model Performance Score\n (Percentile Among Models)")+xlab("Age Basis Offset (years)")+theme(legend.position = "none",strip.background = element_blank(),strip.text.x = (element_text(size = 20)))

###supporting data

offsetboxplot_s_r<-outofsamplebasismodelcompdf[outofsamplebasismodelcompdf$simpnamef=="basis",c("ageoffset","Performancescore","typefull")]
write.csv(offsetboxplot_s_r,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Data/SupportingData/Sup12B.csv")

metafitwithoutdismoothoffsets<-lapply(ageoffsets,function(offsetforvis){
  print(offsetforvis)
  if(sign(offsetforvis)==-1){
    visleadlagfunction<-function(x,n,...){
      n<-n*-1
      dplyr::lead(x,n)
    }
  }else{
    visleadlagfunction<-function(x,n,..){
      dplyr::lag(x,n)
    }
  }
  accmetafitwithoutdi<-fitswithoutdiacc %>% dplyr::group_by(ages) %>% dplyr::summarize(estimate=mean(fit_value)) ###only one dataset here so just average
  accmetafitwithoutdismooth<-mgcvscalefits(accmetafitwithoutdi,outcomevars = c("estimate"),predvars = "ages",mformula = as.formula("outcome~s(pred)"),scale=FALSE,interval_inc = .1) ###fit smooth
  
  offsetvistempacc<-accmetafitwithoutdismooth %>% arrange(pred) %>% mutate(offsetestimate=visleadlagfunction(fit,offsetforvis))
  offsetvistempacc$offset<-offsetforvis
  return(offsetvistempacc)
})
accmetafitwithoutdismoothoffsetsdf<-do.call(rbind,metafitwithoutdismoothoffsets)
ggofsetviz<-ggplot(accmetafitwithoutdismoothoffsetsdf,aes(x=pred,y=offsetestimate,colour=as.factor(offset/10)))+geom_point()+scale_colour_brewer(palette = "Reds")
ggofsetviz<-LNCDR::lunaize(ggofsetviz)+labs(col=guide_legend(title="Basis Offset\n(years)"))+theme(legend.title =element_text(size=15),legend.position = c(0.8, 0.39))+ylab("Basis Value (A.U.)")+xlab("Age (years)")
metafitwithoutdismoothoffsetslatency<-lapply(ageoffsets,function(offsetforvis){
  print(offsetforvis)
  if(sign(offsetforvis)==-1){
    visleadlagfunction<-function(x,n,...){
      n<-n*-1
      dplyr::lead(x,n)
    }
  }else{
    visleadlagfunction<-function(x,n,..){
      dplyr::lag(x,n)
    }
  }
  latmetafitwithoutdi<-fitswithoutdilat %>% dplyr::group_by(ages) %>% dplyr::summarize(estimate=mean(fit_value)) ###only one dataset here so just average
  latmetafitwithoutdismooth<-mgcvscalefits(latmetafitwithoutdi,outcomevars = c("estimate"),predvars = "ages",mformula = as.formula("outcome~s(pred)"),scale=FALSE,interval_inc = .1) ###fit smooth
  
  offsetvistemplat<-latmetafitwithoutdismooth %>% arrange(pred) %>% mutate(offsetestimate=visleadlagfunction(fit,offsetforvis))
  offsetvistemplat$offset<-offsetforvis
  return(offsetvistemplat)
})
latmetafitwithoutdismoothoffsetsdf<-do.call(rbind,metafitwithoutdismoothoffsetslatency)
ggofsetvizlat<-ggplot(latmetafitwithoutdismoothoffsetsdf,aes(x=pred,y=offsetestimate,colour=as.factor(offset/10)))+geom_point()+scale_colour_brewer(palette = "Blues")
ggofsetvizlat<-LNCDR::lunaize(ggofsetvizlat)+labs(col=guide_legend(title="Basis Offset\n(years)"))+theme(legend.title =element_text(size=15),legend.position = c(0.8, .60))+ylab("Basis Value (A.U.)")+xlab("Age (years)")


accmetafitwithoutdismoothoffsetsdf_r<-accmetafitwithoutdismoothoffsetsdf[,c("pred","offsetestimate","offset")]
accmetafitwithoutdismoothoffsetsdf_r$type<-"Accuracy"

latmetafitwithoutdismoothoffsetsdf_r<-latmetafitwithoutdismoothoffsetsdf[,c("pred","offsetestimate","offset")]
latmetafitwithoutdismoothoffsetsdf_r$type<-"Latency"


alloffsetvis<-rbind(accmetafitwithoutdismoothoffsetsdf_r,latmetafitwithoutdismoothoffsetsdf_r)
write.csv(alloffsetvis,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Data/SupportingData/Sup12A.csv")

