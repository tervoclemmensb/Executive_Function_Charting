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
#########
source("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/GeneralScripts/growthrate_mgcvgam.R")
source("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/GeneralScripts/metabyage.R")
###########READ DATA#########
#############################
LUNA<-read.csv("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/data/btc_R03cleaneddata_20220306.outlierremoved.compositeacclat.csv")
LUNA$id<-as.factor(LUNA$id)
LUNA$metaage<-LUNA$Ageatvisit
LUNA$dataset<-"LUNA"
####
NCANDA<-read.csv("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/NCANDA/Data/btc_NCANDAscoredmeasures_20220306.outlierremoved.compositeacclat.csv")
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
#################age range######
minmetaage<-range(LUNA$metaage,NCANDA$metaage,NKI$metaage,PNC$metaage)[1]
maxmetaage<-range(LUNA$metaage,NCANDA$metaage,NKI$metage,PNC$metaage)[2]
################################
agegridforinterp<-data.frame(ages=as.numeric(seq(8,35,by=.1)))
###individual fits#####
###Luna####
lunaformula<-as.formula('outcome~s(pred)+s(visitnum)')
luan2030idx<-which(LUNA$metaage > 20 & LUNA$metaage < 30 & LUNA$metaage < max(LUNA$metaage))
LUNA[,c("Accuracycomposite","Latencycomposite")]<-lapply(LUNA[,c("Accuracycomposite","Latencycomposite")],scalebyidx,idx=luan2030idx)
LUNAfits<-mgcvscalefits(LUNA,outcomevars = c("Accuracycomposite","Latencycomposite"),predvars = "Ageatvisit",idvar="id",mformula = lunaformula,scale=FALSE) ##already z scaled within-sample
LUNAfits$dataset<-"Luna"
LUNAfitsfordots<-mgcvscalefits(LUNA,outcomevars = c("Accuracycomposite","Latencycomposite"),predvars = "Ageatvisit",idvar="id",mformula = lunaformula,scale=FALSE,interval_inc = 3) ##already z scaled within-sample
LUNAfitsfordots$dataset<-"Luna"
LUNAinterpolatedfits<-interpolatebyage(LUNAfits,agegrid=agegridforinterp,agegridvar="ages",agevar="pred",vars=c("Accuracycomposite","Latencycomposite"),longformat = TRUE,
                                       datanamefromlongformat="fit",varnameforlongformat="outcome")
LUNAinterpolatedfitsse<-interpolatebyage(LUNAfits,agegrid=agegridforinterp,agegridvar="ages",agevar="pred",vars=c("Accuracycomposite","Latencycomposite"),longformat = TRUE,
                                       datanamefromlongformat="se",varnameforlongformat="outcome",namemodifier="se")
LUNAinterpolatedfits<-merge(LUNAinterpolatedfits,LUNAinterpolatedfitsse,by="ages")
LUNAinterpolatedfits$dataset<-"LUNA"
###NCANDA#####
ncandaformula<-as.formula('outcome~s(pred)+s(visitnum,k=5)')
NCANDA2030idx<-which(NCANDA$metaage > 20 & NCANDA$metaage < 30 & NCANDA$metaage < max(NCANDA$metaage))
NCANDA[,c("Accuracycomposite","Latencycomposite")]<-lapply(NCANDA[,c("Accuracycomposite","Latencycomposite")],scalebyidx,idx=NCANDA2030idx)
NCANDAfits<-mgcvscalefits(NCANDA,outcomevars = c("Accuracycomposite","Latencycomposite"),predvars = "cnp_age",idvar="id",mformula = ncandaformula,scale=FALSE)
NCANDAfits$dataset<-"NCANDA"
NCANDAfitsfordots<-mgcvscalefits(NCANDA,outcomevars = c("Accuracycomposite","Latencycomposite"),predvars = "cnp_age",idvar="id",mformula = ncandaformula,scale=FALSE,interval_inc = 3)
NCANDAfitsfordots$dataset<-"NCANDA"
NCANDAfitsinterpolatedfits<-interpolatebyage(NCANDAfits,agegrid=agegridforinterp,agegridvar="ages",agevar="pred",vars=c("Accuracycomposite","Latencycomposite"),longformat = TRUE,
                                       datanamefromlongformat="fit",varnameforlongformat="outcome")
NCANDAfitsinterpolatedfitsse<-interpolatebyage(NCANDAfits,agegrid=agegridforinterp,agegridvar="ages",agevar="pred",vars=c("Accuracycomposite","Latencycomposite"),longformat = TRUE,
                                         datanamefromlongformat="se",varnameforlongformat="outcome",namemodifier="se")
NCANDAfitsinterpolatedfits<-merge(NCANDAfitsinterpolatedfits,NCANDAfitsinterpolatedfitsse,by="ages")
NCANDAfitsinterpolatedfits$dataset<-"NCANDA"
####NKI#####
nkiformula<-as.formula('outcome~s(pred)')
NKI2030idx<-which(NKI$metaage > 20 & NKI$metaage < 30 & NKI$metaage < max(NKI$metaage))
NKI[,c("Accuracycomposite","Latencycomposite")]<-lapply(NKI[,c("Accuracycomposite","Latencycomposite")],scalebyidx,idx=NKI2030idx)
NKIfits<-mgcvscalefits(NKI,outcomevars = c("Accuracycomposite","Latencycomposite"),predvars = "age",mformula = nkiformula,scale=FALSE)
NKIfits$dataset<-"NKI"
NKIfitsfordots<-mgcvscalefits(NKI,outcomevars = c("Accuracycomposite","Latencycomposite"),predvars = "age",mformula = nkiformula,scale=FALSE,interval_inc = 3)
NKIfitsfordots$dataset<-"NKI"
NKIfitsinterpolatedfits<-interpolatebyage(NKIfits,agegrid=agegridforinterp,agegridvar="ages",agevar="pred",vars=c("Accuracycomposite","Latencycomposite"),longformat = TRUE,
                                             datanamefromlongformat="fit",varnameforlongformat="outcome")
NKIfitsinterpolatedfitsse<-interpolatebyage(NKIfits,agegrid=agegridforinterp,agegridvar="ages",agevar="pred",vars=c("Accuracycomposite","Latencycomposite"),longformat = TRUE,
                                             datanamefromlongformat="se",varnameforlongformat="outcome",namemodifier="se")
NKIfitsinterpolatedfits<-merge(NKIfitsinterpolatedfits,NKIfitsinterpolatedfitsse,by="ages")
NKIfitsinterpolatedfits$dataset<-"NKI"
###PNC
pncformula<-as.formula('outcome~s(pred)')
PNC2030idx<-which(PNC$metaage > 20 & PNC$metaage < 30 & PNC$metaage < max(PNC$metaage))
PNC[,c("Accuracycomposite","Latencycomposite")]<-lapply(PNC[,c("Accuracycomposite","Latencycomposite")],scalebyidx,idx=PNC2030idx)
PNCfits<-mgcvscalefits(PNC,outcomevars = c("Accuracycomposite","Latencycomposite"),predvars = "age",mformula = pncformula,scale=FALSE)
PNCfits$dataset<-"PNC"
PNCfitsfordots<-mgcvscalefits(PNC,outcomevars = c("Accuracycomposite","Latencycomposite"),predvars = "age",mformula = pncformula,scale=FALSE,interval_inc = 3)
PNCfitsfordots$dataset<-"PNC"
PNCfitsinterpolatedfits<-interpolatebyage(PNCfits,agegrid=agegridforinterp,agegridvar="ages",agevar="pred",vars=c("Accuracycomposite","Latencycomposite"),longformat = TRUE,
                                          datanamefromlongformat="fit",varnameforlongformat="outcome")
PNCfitsinterpolatedfitsse<-interpolatebyage(PNCfits,agegrid=agegridforinterp,agegridvar="ages",agevar="pred",vars=c("Accuracycomposite","Latencycomposite"),longformat = TRUE,
                                            datanamefromlongformat="se",varnameforlongformat="outcome",namemodifier="se")
PNCfitsinterpolatedfits<-merge(PNCfitsinterpolatedfits,PNCfitsinterpolatedfitsse,by="ages")
PNCfitsinterpolatedfits$dataset<-"PNC"
##all fits####
allfits<-plyr::rbind.fill(LUNAfits,NCANDAfits) %>% plyr::rbind.fill(.,NKIfits) %>% plyr::rbind.fill(.,PNCfits)
allinterpolatedfits<-plyr::rbind.fill(LUNAinterpolatedfits,NCANDAfitsinterpolatedfits) %>% plyr::rbind.fill(.,NKIfitsinterpolatedfits) %>% plyr::rbind.fill(.,PNCfitsinterpolatedfits)

###acc meta###
metabyageAcccomposite<-metabyage(allinterpolatedfits,agevar="ages",valcol="Accuracycomposite",secol ="se_Accuracycomposite",datasetvar = "dataset")
metabyageAcccomposite$dataset<-"meta"
metacompform<-as.formula('outcome~s(pred)')
metabyageAcccomposite_smooth<-mgcvscalefits(metabyageAcccomposite,outcomevars = c("estimate"),predvars = "age",mformula = metacompform,scale=FALSE)
metabyageAcccomposite_smooth$dataset<-"meta"

metabyageLatencycomposite<-metabyage(allinterpolatedfits,agevar="ages",valcol="Latencycomposite",secol ="se_Latencycomposite",datasetvar = "dataset")
metabyageLatencycomposite$dataset<-"meta"
metabyageLatencycomposite_smooth<-mgcvscalefits(metabyageLatencycomposite,outcomevars = c("estimate"),predvars = "age",mformula = metacompform,scale=FALSE)
metabyageLatencycomposite_smooth$dataset<-"meta"

###Figure 5A,B
ggallaccfits<-ggplot(allfits[allfits$outcome=="Accuracycomposite",],aes(x=pred,y=fit,fill=dataset))+geom_ribbon(aes(ymin=fit-2*se,ymax=fit+2*se,fill=dataset,alpha=dataset),colour="black")+
  geom_line(colour="black",se=FALSE)+
  scale_fill_manual(values=c("black","black","black","black","black"))+scale_alpha_manual(values=c(.1,.35,.70,1,.80))+
  geom_line(data=metabyageAcccomposite_smooth,aes(x=pred,y=fit),colour="#c9270e",linetype="solid",size=2)
  #+
  #geom_line(data=metabyageAcccomposite,aes(x=age,y=estimate),colour="#c9270e",linetype="dashed",size=.75) ##visualization switched to separate supplemental
ggallaccfits<-LNCDR::lunaize(ggallaccfits)+xlab("Age (years)")+ylab("Executive Function Accuracy\n(adult z-score)")+
  theme(legend.position = "none")

save(ggallaccfits,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/Metabasis/allaccMetabasis.plot.Rdata")
save(metabyageAcccomposite_smooth,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/Metabasis/Metabasis.acc.data.Rdata")


ggalllatfits<-ggplot(allfits[allfits$outcome=="Latencycomposite",],aes(x=pred,y=fit,fill=dataset))+geom_ribbon(aes(ymin=fit-2*se,ymax=fit+2*se,fill=dataset,alpha=dataset),colour="black")+
  geom_line(colour="black",se=FALSE)+
  scale_fill_manual(values=c("black","black","black","black","black"))+scale_alpha_manual(values=c(.1,.35,.70,1,.80))+
  geom_line(data=metabyageLatencycomposite_smooth,aes(x=pred,y=fit),colour="#2f42bd",linetype="solid",size=2)
  #+
  #geom_line(data=metabyageLatencycomposite,aes(x=age,y=estimate),colour="#2f42bd",linetype="dashed",size=.75) ##visualization switched to separate supplemental
latlegend<-cowplot::get_legend(ggalllatfits)
ggalllatfits<-LNCDR::lunaize(ggalllatfits)+xlab("Age (years)")+ylab("Executive Function Latency\n(adult z-score)")+
  theme(legend.position = "none")

save(ggalllatfits,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/Metabasis/alllatMetabasis.plot.Rdata")
save(metabyageLatencycomposite_smooth,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/Metabasis/Metabasis.lat.data.Rdata")

#####save trajectories
metabyageAcccomposite_smoothsave<-metabyageAcccomposite_smooth[,c("pred","fit")]
names(metabyageAcccomposite_smoothsave)<-c("age","accuracy_trajectory_adultzscale")
metabyageAcccomposite_smoothsave$accuracy_trajectory_minmaxscale<-scales::rescale(metabyageAcccomposite_smoothsave$accuracy_trajectory_adultzscale)

metabyageLatencycomposite_smoothsave<-metabyageLatencycomposite_smooth[,c("pred","fit")]
names(metabyageLatencycomposite_smoothsave)<-c("age","latency_trajectory_adultzscale")
metabyageLatencycomposite_smoothsave$latency_trajectory_minmaxscale<-scales::rescale(metabyageLatencycomposite_smoothsave$latency_trajectory_adultzscale)

metabasissaveall<-merge(metabyageAcccomposite_smoothsave,metabyageLatencycomposite_smoothsave,by=c("age"))
write.csv(metabasissaveall,"~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/Metabasis/Metabasis.acclat.data.csv")

####visualize against inverse, linear, quadratic for talk

###lmers
allfits$invpred<-1/allfits$pred
acctempvis<-allfits[allfits$outcome=="Accuracycomposite",]
preddata<-data.frame(invpred=seq(min(acctempvis$invpred),max(acctempvis$invpred),by=.001))
preddata$pred<-1/preddata$invpred
preddata$predsq<-preddata$pred^2
acctempvis$predsq<-acctempvis$pred^2
lmerinversevisfit<-lme4::lmer(fit~invpred + (invpred|dataset),data=acctempvis)
preddata$invpred<-predict(lmerinversevisfit,newdata=preddata,re.form=NA)
lmerlinfit<-lme4::lmer(fit~pred + (pred|dataset),data=acctempvis)
preddata$linpred<-predict(lmerlinfit,newdata=preddata,re.form=NA)
lmerquadfit<-lme4::lmer(fit~pred+predsq+(pred|dataset),data=acctempvis)
preddata$sqpred<-predict(lmerquadfit,newdata=preddata,re.form=NA)



ggotherfitvis<-ggplot()+
  geom_line(colour="black")+
  scale_fill_manual(values=c("black","black","black","black","black"))+scale_alpha_manual(values=c(.1,.35,.70,1,.80))+
  geom_line(data=metabyageAcccomposite_smooth,aes(x=pred,y=fit,linetype="basis"),colour="#c9270e",size=2)+
  geom_line(data=preddata,aes(x=pred,y=linpred,linetype="linear"),colour="#c9270e",size=1)+
  geom_line(data=preddata,aes(x=pred,y=invpred,linetype="inverse"),colour="#c9270e",size=1)+
  geom_line(data=preddata,aes(x=pred,y=sqpred,linetype="quadratic"),colour="#c9270e",size=1)+
  scale_linetype_manual(name = "Age Model", values = c("basis" = "solid", "linear" = "dashed","inverse"="dotted","quadratic"="twodash"))

ggotherfitvis<-LNCDR::lunaize(ggotherfitvis)+xlab("Age (years)")+ylab("Executive Function \n (z-score)")+theme(legend.position = "top",legend.title = element_blank())
  
 # geom_smooth(data=preddata,aes(x=pred,y=fit),method="lm",formula = y=1/x,colour="red",linetype="dashed",se=FALSE)








