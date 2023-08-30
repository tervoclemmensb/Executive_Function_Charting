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
LUNA<-read.csv("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/data/btc_R03cleaneddata_20220306.outlierremoved.compositeacclat.2030scale.csv")
LUNA$id<-as.factor(LUNA$id)
LUNA$metaage<-LUNA$Ageatvisit
LUNA$dataset<-"LUNA"
LUNAsex<-read.csv("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Data/gender_20220214.csv")
LUNAsex$sex[LUNAsex$Gender=="Female"]<-"F" ###coded as "gender" but choices were limited to male female; framed as sex
LUNAsex$sex[LUNAsex$Gender=="Male"]<-"M"
Lunasexadd<-read.csv("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Data/gender_20220214_manualadd.csv",sep=" ",header=FALSE)
names(Lunasexadd)<-c("id","X","Y","sex","Z")
LUNAsex<-plyr::rbind.fill(LUNAsex,Lunasexadd[c("id","sex")])
LUNAsex<-LUNAsex %>% dplyr::group_by(id) %>% dplyr::summarise(sex=unique(sex),sexn=length(unique(sex)))
LUNAsex<-LUNAsex[LUNAsex$id %in% unique(LUNA$id),]
length(unique(LUNAsex$id[LUNAsex$sexn==2]))
LUNAsex<-LUNAsex[LUNAsex$sexn==1,] ###limited for analysis but reported in manuscript
LUNA<-merge(LUNA,LUNAsex[,c("id","sex")],by=c("id"),all.x=TRUE)
LUNA<-LUNA[!duplicated(LUNA),]
####
NCANDA<-read.csv("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/NCANDA/data/btc_R03cleaneddata_20220306.outlierremoved.compositeacclat.2030scale.csv")
NCANDA$id<-as.factor(NCANDA$subject)
NCANDA$metaage<-NCANDA$cnp_age
NCANDA$visitnum<-as.numeric(as.factor(NCANDA$visit))
NCANDA$dataset<-"NCANDA"
NCANDAsex<-read.csv("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/NCANDA_Behavioral/Data/NCANDA_RELEASE_4Y_REDCAP_MEASUREMENTS_V02/summaries/redcap/demographics.csv")
NCANDAsex<-NCANDAsex %>% dplyr::group_by(subject) %>% dplyr::summarise(sex=unique(sex))
NCANDA<-merge(NCANDA,NCANDAsex,by=("subject"),all.x=TRUE)
####
NKI<-read.csv("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/NKI/Data/btc_NKIscoredmeasures_20220306.outlierremoved.compositeacclat.2030scale.csv")
NKI$id<-as.factor(NKI$subject)
NKI$metaage<-NKI$age
NKI$dataset<-"NKI"
NKI$visitnum<-1 ##cross-sectional
NKIsex<-read.csv("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/NKI/data/NKIdemodata-2022-03-03T18_51_26.161Z.csv")
#dem_002 (Male=1,Female=2)
NKIsex$sex<-NKIsex$dem_002
NKIsex$sex[NKIsex$sex=="1"]<-"M"
NKIsex$sex[NKIsex$sex=="2"]<-"F"
NKIsex$sex[NKIsex$sex=="MD"]<-NA
NKIsex$subject<-unlist(lapply(strsplit(NKIsex$Identifiers,"[,]"),"[[",1))
NKIsex<-NKIsex %>% dplyr::group_by(subject) %>% dplyr::summarise(sex=unique(sex))
NKI<-merge(NKI,NKIsex[,c("subject","sex")],by=c("subject"),all.x=TRUE)
####
PNC<-read.csv("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/PNC/Data/btc_R03cleaneddata_20220306.outlierremoved.compositeacclat.2030scale.csv")
PNC$id<-as.factor(PNC$dbGaP_Subject_ID)
PNC$metaage<-PNC$age
PNC$dataset<-"PNC"
PNC$visitnum<-1 ##cross-sectional
PNCsex<-read.csv("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/PNC/Data/PNCsex.csv")
names(PNCsex)[names(PNCsex)=="Sex"]<-"sex"
PNCsex<-PNCsex %>% dplyr::group_by(dbGaP_Subject_ID) %>% dplyr::summarise(sex=unique(sex))
PNC<-merge(PNC,PNCsex,by="dbGaP_Subject_ID",all.x=TRUE)
PNC$sex[PNC$sex %in% c("+")]<-NA ###set + to NA to not bin with M/F (NA and + reported separately in manuscript)
####################
###note all internal scaling set to false due to subsetting by sex and want to keep on original, full sample z units of composite metrics
###without internal scaling derivative plots will revert to T-scaling; saturating this scale here such that significant always shows same colour to minimize plotting complexity (i.e., is significant male, is significant female )
#######models#######
###LUNA
lunaformula<-as.formula('outcome~s(pred)+s(visitnum)')
lunaformula_m<-as.formula('outcome~s(pred)+s(visitnum,k=9)') ###9 visit max for male participants
lunascaledfitM<-mgcvscalefits(LUNA[LUNA$sex=="M",],outcomevars = c("Accuracycomposite","Latencycomposite"),predvars = "Ageatvisit",idvar="id",mformula = lunaformula_m,scale=FALSE)
lunascaledfitM$sex<-"M"
lunaderivMacc<-multi_mgcvgam_growthrate_multiplot_rasteronly(df=LUNA[LUNA$sex=="M",],outcomevars = c("Accuracycomposite"),idvar="id",predvars='Ageatvisit',mformula=lunaformula_m,datarangegrey=TRUE,derivrangemanual=c(-.2,.2),derivcolourmanual=c("grey33","#00BFC4"),devageguides=FALSE,zscale=FALSE) ###"M" colour key corresponds to expected direction for acc/lat 
lunaderivMlat<-multi_mgcvgam_growthrate_multiplot_rasteronly(df=LUNA[LUNA$sex=="M",],outcomevars = c("Latencycomposite"),idvar="id",predvars='Ageatvisit',mformula=lunaformula_m,datarangegrey=TRUE,derivrangemanual=c(-.2,.2),derivcolourmanual=c("#00BFC4","grey33"),devageguides=FALSE,zscale=FALSE) ###"M" colour key corresponds to expected direction for acc/lat 

lunascaledfitsF<-mgcvscalefits(LUNA[LUNA$sex=="F",],outcomevars = c("Accuracycomposite","Latencycomposite"),predvars = "Ageatvisit",idvar="id",mformula = lunaformula,scale=FALSE)
lunascaledfitsF$sex<-"F"
lunaderivFacc<-multi_mgcvgam_growthrate_multiplot_rasteronly(df=LUNA[LUNA$sex=="F",],outcomevars = c("Accuracycomposite"),idvar="id",predvars='Ageatvisit',mformula=lunaformula,datarangegrey=TRUE,derivrangemanual=c(-.2,.2),derivcolourmanual=c("grey33","#F8766D"),devageguides=FALSE,zscale=FALSE) ###"F" colour key corresponds to expected direction for acc/lat 
lunaderivFlat<-multi_mgcvgam_growthrate_multiplot_rasteronly(df=LUNA[LUNA$sex=="F",],outcomevars = c("Latencycomposite"),idvar="id",predvars='Ageatvisit',mformula=lunaformula,datarangegrey=TRUE,derivrangemanual=c(-.2,.2),derivcolourmanual=c("#F8766D","grey33"),devageguides=FALSE,zscale=FALSE) ###"F" colour key corresponds to expected direction for acc/lat 

alllunafits<-plyr::rbind.fill(lunascaledfitM,lunascaledfitsF)
alllunafits$outcomelabel<-dplyr::if_else(alllunafits$outcome=="Accuracycomposite","Accuracy","Latency")

gglunaalllatfitsacc<-ggplot(alllunafits,aes(x=pred,y=fit,colour=sex,fill=sex))+geom_line(colour="black")+geom_ribbon(aes(ymin=fit-2*se,ymax=fit+2*se,fill=sex),colour="black",alpha=.5)+
  theme(legend.position = "top")+
  facet_grid(cols=vars(outcomelabel),scales="free")+theme(strip.text.y= element_text(angle = 180))+
  scale_x_continuous(limits = c(8, 36),breaks=(c(10,15,20,25,30,35)))
gglunaalllatfitsacc<-LNCDR::lunaize(gglunaalllatfitsacc)+xlab("Age (years)")+ylab("Executive Function\n (z-score)")+theme(legend.position = "none",legend.title = element_blank())+
  theme(strip.text.y= element_text(angle = 360))+theme(strip.background = element_blank())

lunaall<-(gglunaalllatfitsacc)/(lunaderivFacc$returnplot+lunaderivFlat$returnplot)/(lunaderivMacc$returnplot+lunaderivMlat$returnplot)+patchwork::plot_layout(heights = c(3,1,1))
ggsave(lunaall,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/Sensplots/lunasexsens.pdf",height=6,width=8)
#####NCANDA######
NCANDAformula<-as.formula('outcome~s(pred)+s(visitnum,k=5)')
NCANDAscaledfitM<-mgcvscalefits(NCANDA[NCANDA$sex=="M",],outcomevars = c("Accuracycomposite","Latencycomposite"),predvars = "cnp_age",idvar="id",mformula = NCANDAformula,scale=FALSE)
NCANDAscaledfitM$sex<-"M"
NCANDAderivMacc<-multi_mgcvgam_growthrate_multiplot_rasteronly(df=NCANDA[NCANDA$sex=="M",],outcomevars = c("Accuracycomposite"),idvar="id",predvars='cnp_age',mformula=NCANDAformula,datarangegrey=TRUE,derivrangemanual=c(-.2,.2),derivcolourmanual=c("grey33","#00BFC4"),devageguides=FALSE,zscale=FALSE) ###"M" colour key corresponds to expected direction for acc/lat 
NCANDAderivMlat<-multi_mgcvgam_growthrate_multiplot_rasteronly(df=NCANDA[NCANDA$sex=="M",],outcomevars = c("Latencycomposite"),idvar="id",predvars='cnp_age',mformula=NCANDAformula,datarangegrey=TRUE,derivrangemanual=c(-.2,.2),derivcolourmanual=c("#00BFC4","grey33"),devageguides=FALSE,zscale=FALSE) ###"M" colour key corresponds to expected direction for acc/lat 

NCANDAscaledfitsF<-mgcvscalefits(NCANDA[NCANDA$sex=="F",],outcomevars = c("Accuracycomposite","Latencycomposite"),predvars = "cnp_age",idvar="id",mformula = NCANDAformula,scale=FALSE)
NCANDAscaledfitsF$sex<-"F"
NCANDAderivFacc<-multi_mgcvgam_growthrate_multiplot_rasteronly(df=NCANDA[NCANDA$sex=="F",],outcomevars = c("Accuracycomposite"),idvar="id",predvars='cnp_age',mformula=NCANDAformula,datarangegrey=TRUE,derivrangemanual=c(-.2,.2),derivcolourmanual=c("grey33","#F8766D"),devageguides=FALSE,zscale=FALSE) ###"F" colour key corresponds to expected direction for acc/lat 
NCANDAderivFlat<-multi_mgcvgam_growthrate_multiplot_rasteronly(df=NCANDA[NCANDA$sex=="F",],outcomevars = c("Latencycomposite"),idvar="id",predvars='cnp_age',mformula=NCANDAformula,datarangegrey=TRUE,derivrangemanual=c(-.2,.2),derivcolourmanual=c("#F8766D","grey33"),devageguides=FALSE,zscale=FALSE) ###"F" colour key corresponds to expected direction for acc/lat 

allNCANDAfits<-plyr::rbind.fill(NCANDAscaledfitM,NCANDAscaledfitsF)
allNCANDAfits$outcomelabel<-dplyr::if_else(allNCANDAfits$outcome=="Accuracycomposite","Accuracy","Latency")

ggNCANDAalllatfitsacc<-ggplot(allNCANDAfits,aes(x=pred,y=fit,colour=sex,fill=sex))+geom_line(colour="black")+geom_ribbon(aes(ymin=fit-2*se,ymax=fit+2*se,fill=sex),colour="black",alpha=.5)+
  theme(legend.position = "top")+
  facet_grid(cols=vars(outcomelabel),scales="free")+theme(strip.text.y= element_text(angle = 180))+
  scale_x_continuous(limits = c(8, 36),breaks=(c(10,15,20,25,30,35)))
ggNCANDAalllatfitsacc<-LNCDR::lunaize(ggNCANDAalllatfitsacc)+xlab("Age (years)")+ylab("Executive Function\n (z-score)")+theme(legend.position = "none",legend.title = element_blank())+
  theme(strip.text.y= element_text(angle = 360))+theme(strip.background = element_blank())

NCANDAall<-(ggNCANDAalllatfitsacc)/(NCANDAderivFacc$returnplot+NCANDAderivFlat$returnplot)/(NCANDAderivMacc$returnplot+NCANDAderivMlat$returnplot)+patchwork::plot_layout(heights = c(3,1,1))
ggsave(NCANDAall,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/Sensplots/NCANDAsexsens.pdf",height=6,width=8)

####NKI####
NKIformula<-as.formula('outcome~s(pred)')
NKIscaledfitM<-mgcvscalefits(NKI[NKI$sex=="M",],outcomevars = c("Accuracycomposite","Latencycomposite"),predvars = "age",mformula = NKIformula,scale=FALSE)
NKIscaledfitM$sex<-"M"
NKIderivMacc<-multi_mgcvgam_growthrate_multiplot_rasteronly(df=NKI[NKI$sex=="M",],outcomevars = c("Accuracycomposite"),predvars='age',mformula=NKIformula,datarangegrey=TRUE,derivrangemanual=c(-.2,.2),derivcolourmanual=c("grey33","#00BFC4"),devageguides=FALSE,zscale=FALSE) ###"M" colour key corresponds to expected direction for acc/lat 
NKIderivMlat<-multi_mgcvgam_growthrate_multiplot_rasteronly(df=NKI[NKI$sex=="M",],outcomevars = c("Latencycomposite"),predvars='age',mformula=NKIformula,datarangegrey=TRUE,derivrangemanual=c(-.2,.2),derivcolourmanual=c("#00BFC4","grey33"),devageguides=FALSE,zscale=FALSE) ###"M" colour key corresponds to expected direction for acc/lat 

NKIscaledfitsF<-mgcvscalefits(NKI[NKI$sex=="F",],outcomevars = c("Accuracycomposite","Latencycomposite"),predvars = "age",mformula = NKIformula,scale=FALSE)
NKIscaledfitsF$sex<-"F"
NKIderivFacc<-multi_mgcvgam_growthrate_multiplot_rasteronly(df=NKI[NKI$sex=="F",],outcomevars = c("Accuracycomposite"),predvars='age',mformula=NKIformula,datarangegrey=TRUE,derivrangemanual=c(-.2,.2),derivcolourmanual=c("grey33","#F8766D"),devageguides=FALSE,zscale=FALSE) ###"F" colour key corresponds to expected direction for acc/lat 
NKIderivFlat<-multi_mgcvgam_growthrate_multiplot_rasteronly(df=NKI[NKI$sex=="F",],outcomevars = c("Latencycomposite"),predvars='age',mformula=NKIformula,datarangegrey=TRUE,derivrangemanual=c(-.2,.2),derivcolourmanual=c("#F8766D","grey33"),devageguides=FALSE,zscale=FALSE) ###"F" colour key corresponds to expected direction for acc/lat 

allNKIfits<-plyr::rbind.fill(NKIscaledfitM,NKIscaledfitsF)
allNKIfits$outcomelabel<-dplyr::if_else(allNKIfits$outcome=="Accuracycomposite","Accuracy","Latency")
ggNKIalllatfitsacc<-ggplot(allNKIfits,aes(x=pred,y=fit,colour=sex,fill=sex))+geom_line(colour="black")+geom_ribbon(aes(ymin=fit-2*se,ymax=fit+2*se,fill=sex),colour="black",alpha=.5)+
  theme(legend.position = "top")+
  facet_grid(cols=vars(outcomelabel),scales="free")+theme(strip.text.y= element_text(angle = 180))+
  scale_x_continuous(limits = c(8, 36),breaks=(c(10,15,20,25,30,35)))
ggNKIalllatfitsacc<-LNCDR::lunaize(ggNKIalllatfitsacc)+xlab("Age (years)")+ylab("Executive Function\n (z-score)")+theme(legend.position = "none",legend.title = element_blank())+
  theme(strip.text.y= element_text(angle = 360))+theme(strip.background = element_blank())

NKIall<-(ggNKIalllatfitsacc)/(NKIderivFacc$returnplot+NKIderivFlat$returnplot)/(NKIderivMacc$returnplot+NKIderivMlat$returnplot)+patchwork::plot_layout(heights = c(3,1,1))
ggsave(NKIall,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/Sensplots/NKIsexsens.pdf",height=6,width=8)

####PNC####
PNCformula<-as.formula('outcome~s(pred)')
PNCscaledfitM<-mgcvscalefits(PNC[PNC$sex=="M",],outcomevars = c("Accuracycomposite","Latencycomposite"),predvars = "age",mformula = PNCformula,scale=FALSE)
PNCscaledfitM$sex<-"M"
PNCderivMacc<-multi_mgcvgam_growthrate_multiplot_rasteronly(df=PNC[PNC$sex=="M",],outcomevars = c("Accuracycomposite"),predvars='age',mformula=PNCformula,datarangegrey=TRUE,derivrangemanual=c(-.2,.2),derivcolourmanual=c("grey33","#00BFC4"),devageguides=FALSE,zscale=FALSE) ###"M" colour key corresponds to expected direction for acc/lat 
PNCderivMlat<-multi_mgcvgam_growthrate_multiplot_rasteronly(df=PNC[PNC$sex=="M",],outcomevars = c("Latencycomposite"),predvars='age',mformula=PNCformula,datarangegrey=TRUE,derivrangemanual=c(-.2,.2),derivcolourmanual=c("#00BFC4","grey33"),devageguides=FALSE,zscale=FALSE) ###"M" colour key corresponds to expected direction for acc/lat 

PNCscaledfitsF<-mgcvscalefits(PNC[PNC$sex=="F",],outcomevars = c("Accuracycomposite","Latencycomposite"),predvars = "age",mformula = PNCformula,scale=FALSE)
PNCscaledfitsF$sex<-"F"
PNCderivFacc<-multi_mgcvgam_growthrate_multiplot_rasteronly(df=PNC[PNC$sex=="F",],outcomevars = c("Accuracycomposite"),predvars='age',mformula=PNCformula,datarangegrey=TRUE,derivrangemanual=c(-.2,.2),derivcolourmanual=c("grey33","#F8766D"),devageguides=FALSE,zscale=FALSE) ###"F" colour key corresponds to expected direction for acc/lat 
PNCderivFlat<-multi_mgcvgam_growthrate_multiplot_rasteronly(df=PNC[PNC$sex=="F",],outcomevars = c("Latencycomposite"),predvars='age',mformula=PNCformula,datarangegrey=TRUE,derivrangemanual=c(-.2,.2),derivcolourmanual=c("#F8766D","grey33"),devageguides=FALSE,zscale=FALSE) ###"F" colour key corresponds to expected direction for acc/lat 

allPNCfits<-plyr::rbind.fill(PNCscaledfitM,PNCscaledfitsF)
allPNCfits$outcomelabel<-dplyr::if_else(allPNCfits$outcome=="Accuracycomposite","Accuracy","Latency")
ggPNCalllatfitsacc<-ggplot(allPNCfits,aes(x=pred,y=fit,colour=sex,fill=sex))+geom_line(colour="black")+geom_ribbon(aes(ymin=fit-2*se,ymax=fit+2*se,fill=sex),colour="black",alpha=.5)+
  theme(legend.position = "top")+
  facet_grid(cols=vars(outcomelabel),scales="free")+theme(strip.text.y= element_text(angle = 180))+
  scale_x_continuous(limits = c(8, 36),breaks=(c(10,15,20,25,30,35)))
ggPNCalllatfitsacc<-LNCDR::lunaize(ggPNCalllatfitsacc)+xlab("Age (years)")+ylab("Executive Function\n (z-score)")+theme(legend.position = "none",legend.title = element_blank())+
  theme(strip.text.y= element_text(angle = 360))+theme(strip.background = element_blank())

library(patchwork)
PNCall<-(ggPNCalllatfitsacc)/(PNCderivFacc$returnplot+PNCderivFlat$returnplot)/(PNCderivMacc$returnplot+PNCderivMlat$returnplot)+patchwork::plot_layout(heights = c(3,1,1))
ggsave(PNCall,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/Sensplots/PNCsexsens.pdf",height=6,width=8)

#########all panels#########
lunaallwithNCANDAall<-(gglunaalllatfitsacc+ggNCANDAalllatfitsacc)/(lunaderivFacc$returnplot+lunaderivFlat$returnplot+NCANDAderivFacc$returnplot+NCANDAderivFlat$returnplot)/(lunaderivMacc$returnplot+lunaderivMlat$returnplot+NCANDAderivMacc$returnplot+NCANDAderivMlat$returnplot)+patchwork::plot_layout(heights = c(3,1,1))
#####

###supporting data
alllunafits$dataset<-"Luna"
allNCANDAfits$dataset<-"NCANDA"
allNKIfits$dataset<-"NKI"
allPNCfits$dataset<-"PNC"

outvars<-c("pred","fit","sex","se","outcomelabel","dataset")

supdatafits<-rbind(alllunafits[,outvars],allNCANDAfits[,outvars]) %>% rbind(.,allNKIfits[,outvars]) %>%
  rbind(.,allPNCfits[,outvars])

write.csv(supdatafits,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Data/SupportingData/Sup13.csv")


