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
####
NCANDA<-read.csv("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/NCANDA/data/btc_R03cleaneddata_20220306.outlierremoved.compositeacclat.2030scale.csv")
NCANDA$id<-as.factor(NCANDA$subject)
NCANDA$metaage<-NCANDA$cnp_age
NCANDA$visitnum<-as.numeric(as.factor(NCANDA$visit))
NCANDA$dataset<-"NCANDA"
NCANDAdemo<-read.csv("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/NCANDA_Behavioral/Data/NCANDA_RELEASE_4Y_REDCAP_MEASUREMENTS_V02/summaries/redcap/demographics.csv")
NCANDAclinical<-read.csv("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/NCANDA_Behavioral/Data/NCANDA_RELEASE_4Y_REDCAP_MEASUREMENTS_V02/summaries/redcap/clinical.csv")
diagvars<-grep("diag|abuse",grep("youth",grep("dsm",names(NCANDAclinical),value=TRUE),value=TRUE),value=TRUE)
diagvars<-diagvars[!grepl("label",diagvars)]
diagvars<-diagvars[!diagvars %in% c("lssaga_dsm4_youth_d03_diag","lssaga_dsm5_youth_d03_diag")] ###exclude nicotine dependence
NCANDAclinicaldatadict<-read.csv("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/NCANDA_Behavioral/Data/NCANDA_RELEASE_4Y_REDCAP_MEASUREMENTS_V02/datadict/redcap/clinical_datadict.csv")
NCANDAclinicaldatadict$Field.Label[NCANDAclinicaldatadict$Variable...Field.Name %in% diagvars]

NCANDAclinical$anydxcount<-unlist(lapply(1:nrow(NCANDAclinical),function(cli){
  clindatacli<-NCANDAclinical[cli,diagvars]
  return(length(which(clindatacli==1)))
}))
NCANDA<-merge(NCANDA,NCANDAclinical[,c("subject","anydxcount","visit")],by=c("subject","visit"))
NCANDA$anypsydx<-dplyr::if_else(NCANDA$anydxcount>0,1,0)
####
NKI<-read.csv("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/NKI/Data/btc_NKIscoredmeasures_20220306.outlierremoved.compositeacclat.2030scale.csv")
NKI$id<-as.factor(NKI$subject)
NKI$metaage<-NKI$age
NKI$dataset<-"NKI"
NKI$visitnum<-1 ##cross-sectional
NKIpsy<-read.csv("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/NKI/Data/PSYNKIdata-2022-02-14T22_08_40.494Z.csv")
NKIpsy$id<-unlist(lapply(strsplit(NKIpsy$Identifiers,","),"[[",1))
NKIpsy$Identifiers<-NULL
NKIpsy$anypsydx<-dplyr::if_else(NKIpsy$diag_01.desc=="No Diagnosis or Condition on Axis I",0,1)
NKIpsysum<-NKIpsy %>% group_by(id) %>% dplyr::summarise(anypsydx=sum(anypsydx,na.rm=TRUE)) ###sum so any diagnosis also counted
NKIpsysum$anypsydx<-dplyr::if_else(NKIpsysum$anypsydx>0,1,0)
NKI<-merge(NKI,NKIpsysum,by=c("id"),all.x=TRUE)
####
PNC<-read.csv("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/PNC/Data/btc_R03cleaneddata_20220306.outlierremoved.compositeacclat.2030scale.csv")
PNC$id<-as.factor(PNC$dbGaP_Subject_ID)
PNC$metaage<-PNC$age
PNC$dataset<-"PNC"
PNC$visitnum<-1 ##cross-sectional
PNCpsy<-read.csv("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/PNC/Data/PNCpsy.csv")
PNCpsy$SCR007_bin<-dplyr::if_else(PNCpsy$SCR007==1,1,0)
PNCpsy$SCR006_bin<-dplyr::if_else(PNCpsy$SCR006==1,1,0)
PNCpsy$LTNsum<-rowSums(PNCpsy[,c("SCR007_bin","SCR006_bin")],na.rm=TRUE)
PNCpsy$LTN<-dplyr::if_else(PNCpsy$LTNsum>0,1,0)
PNCpsy<-PNCpsy %>% dplyr::group_by(dbGaP_Subject_ID) %>% dplyr::summarise(LTN=sum(LTN,na.rm=TRUE))###sum so any diagnosis also counted
PNCpsy$LTN<-dplyr::if_else(PNCpsy$LTN>0,1,0)
PNC<-merge(PNC,PNCpsy,by="dbGaP_Subject_ID",all.x=TRUE)
####################
###note all internal scaling set to false due to subsetting by sex and want to keep on original, full sample z units of composite metrics
###without internal scaling derivative plots will revert to T-scaling; saturating this scale here such that significant always shows same colour to minimize plotting complexity (i.e., is significant male, is significant female )
#######models#######
###LUNA
lunaformula<-as.formula('outcome~s(pred)+s(visitnum)')
lunascaledfitfull<-mgcvscalefits(LUNA,outcomevars = c("Accuracycomposite","Latencycomposite"),predvars = "Ageatvisit",idvar="id",mformula = lunaformula,scale=FALSE)
lunascaledfitfull$sampletyp<-"full"
lunaderivfullacc<-multi_mgcvgam_growthrate_multiplot_rasteronly(df=LUNA,outcomevars = c("Accuracycomposite"),idvar="id",predvars='Ageatvisit',mformula=lunaformula,datarangegrey=TRUE,derivrangemanual=c(-.2,.2),derivcolourmanual=c("grey33","grey33"),devageguides=FALSE,zscale=FALSE) ###"M" colour key corresponds to expected direction for acc/lat
lunaderivfulllat<-multi_mgcvgam_growthrate_multiplot_rasteronly(df=LUNA,outcomevars = c("Latencycomposite"),idvar="id",predvars='Ageatvisit',mformula=lunaformula,datarangegrey=TRUE,derivrangemanual=c(-.2,.2),derivcolourmanual=c("grey33","grey33"),devageguides=FALSE,zscale=FALSE) ###"M" colour key corresponds to expected direction for acc/lat

alllunafits<-lunascaledfitfull
alllunafits$outcomelabel<-dplyr::if_else(alllunafits$outcome=="Accuracycomposite","Accuracy","Latency")

gglunaalllatfitsacc<-ggplot(alllunafits,aes(x=pred,y=fit))+geom_line(colour="black")+geom_ribbon(aes(ymin=fit-2*se,ymax=fit+2*se),colour="black",alpha=.5)+
  theme(legend.position = "top")+
  facet_grid(cols=vars(outcomelabel),scales="free")+theme(strip.text.y= element_text(angle = 180))+
  scale_x_continuous(limits = c(8, 36),breaks=(c(10,15,20,25,30,35)))
gglunaalllatfitsacc<-LNCDR::lunaize(gglunaalllatfitsacc)+xlab("Age (years)")+ylab("Executive Function\n (z-score)")+theme(legend.position = "none",legend.title = element_blank())+
  theme(strip.text.y= element_text(angle = 360))+theme(strip.background = element_blank())

lunaall<-(gglunaalllatfitsacc)/(lunaderivfullacc$returnplot+lunaderivfulllat$returnplot)/(lunaderivfullacc$returnplot+lunaderivfulllat$returnplot)+patchwork::plot_layout(heights = c(3,1,1))
ggsave(lunaall,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/Sensplots/lunapsysens.pdf",height=6,width=8)
#####NCANDA######
NCANDAformula<-as.formula('outcome~s(pred)+s(visitnum,k=5)')
NCANDAscaledfitnoclin<-mgcvscalefits(NCANDA[NCANDA$anypsydx!=1,],outcomevars = c("Accuracycomposite","Latencycomposite"),predvars = "cnp_age",idvar="id",mformula = NCANDAformula,scale=FALSE)
NCANDAscaledfitnoclin$sampletyp<-"noclin"
NCANDAscaledfitnoclinacc<-multi_mgcvgam_growthrate_multiplot_rasteronly(df=NCANDA[NCANDA$anypsydx!=1,],outcomevars = c("Accuracycomposite"),idvar="id",predvars='cnp_age',mformula=NCANDAformula,datarangegrey=TRUE,derivrangemanual=c(-.2,.2),derivcolourmanual=c("grey33","#6a50a3"),devageguides=FALSE,zscale=FALSE) ###"M" colour key corresponds to expected direction for acc/lat 
NCANDAscaledfitnoclinlat<-multi_mgcvgam_growthrate_multiplot_rasteronly(df=NCANDA[NCANDA$anypsydx!=1,],outcomevars = c("Latencycomposite"),idvar="id",predvars='cnp_age',mformula=NCANDAformula,datarangegrey=TRUE,derivrangemanual=c(-.2,.2),derivcolourmanual=c("#6a50a3","grey33"),devageguides=FALSE,zscale=FALSE) ###"M" colour key corresponds to expected direction for acc/lat 

NCANDAscaledfitsF<-mgcvscalefits(NCANDA,outcomevars = c("Accuracycomposite","Latencycomposite"),predvars = "cnp_age",idvar="id",mformula = NCANDAformula,scale=FALSE)
NCANDAscaledfitsF$sampletyp<-"full"
NCANDAderivFacc<-multi_mgcvgam_growthrate_multiplot_rasteronly(df=NCANDA,outcomevars = c("Accuracycomposite"),idvar="id",predvars='cnp_age',mformula=NCANDAformula,datarangegrey=TRUE,derivrangemanual=c(-.2,.2),derivcolourmanual=c("grey33","grey33"),devageguides=FALSE,zscale=FALSE) ###"F" colour key corresponds to expected direction for acc/lat 
NCANDAderivFlat<-multi_mgcvgam_growthrate_multiplot_rasteronly(df=NCANDA,outcomevars = c("Latencycomposite"),idvar="id",predvars='cnp_age',mformula=NCANDAformula,datarangegrey=TRUE,derivrangemanual=c(-.2,.2),derivcolourmanual=c("grey33","grey33"),devageguides=FALSE,zscale=FALSE) ###"F" colour key corresponds to expected direction for acc/lat 

allNCANDAfits<-plyr::rbind.fill(NCANDAscaledfitnoclin,NCANDAscaledfitsF)
allNCANDAfits$outcomelabel<-dplyr::if_else(allNCANDAfits$outcome=="Accuracycomposite","Accuracy","Latency")

ggNCANDAalllatfitsacc<-ggplot(allNCANDAfits,aes(x=pred,y=fit,colour=sampletyp,fill=sampletyp))+geom_line(colour="black")+geom_ribbon(aes(ymin=fit-2*se,ymax=fit+2*se,fill=sampletyp),colour="black",alpha=.5)+
  theme(legend.position = "top")+
  facet_grid(cols=vars(outcomelabel),scales="free")+theme(strip.text.y= element_text(angle = 180))+
  scale_x_continuous(limits = c(8, 36),breaks=(c(10,15,20,25,30,35)))+scale_colour_manual(values=c("grey33","#6a50a3"))+scale_fill_manual(values=c("grey33","#6a50a3"))
ggNCANDAalllatfitsacc<-LNCDR::lunaize(ggNCANDAalllatfitsacc)+xlab("Age (years)")+ylab("Executive Function\n (z-score)")+theme(legend.position = "none",legend.title = element_blank())+
  theme(strip.text.y= element_text(angle = 360))+theme(strip.background = element_blank())

NCANDAall<-(ggNCANDAalllatfitsacc)/(NCANDAderivFacc$returnplot+NCANDAderivFlat$returnplot)/(NCANDAscaledfitnoclinacc$returnplot+NCANDAscaledfitnoclinlat$returnplot)+patchwork::plot_layout(heights = c(3,1,1))
ggsave(NCANDAall,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/Sensplots/NCANDApsysens.pdf",height=6,width=8)

####NKI####
NKIformula<-as.formula('outcome~s(pred)')
NKIscaledfitnoclin<-mgcvscalefits(NKI[NKI$anypsydx!=1,],outcomevars = c("Accuracycomposite","Latencycomposite"),predvars = "metaage",mformula = NKIformula,scale=FALSE)
NKIscaledfitnoclin$sampletyp<-"noclin"
NKIscaledfitnoclinacc<-multi_mgcvgam_growthrate_multiplot_rasteronly(df=NKI[NKI$anypsydx!=1,],outcomevars = c("Accuracycomposite"),predvars='metaage',mformula=NKIformula,datarangegrey=TRUE,derivrangemanual=c(-.2,.2),derivcolourmanual=c("grey33","#6a50a3"),devageguides=FALSE,zscale=FALSE) ###"M" colour key corresponds to expected direction for acc/lat 
NKIscaledfitnoclinlat<-multi_mgcvgam_growthrate_multiplot_rasteronly(df=NKI[NKI$anypsydx!=1,],outcomevars = c("Latencycomposite"),predvars='metaage',mformula=NKIformula,datarangegrey=TRUE,derivrangemanual=c(-.2,.2),derivcolourmanual=c("#6a50a3","grey33"),devageguides=FALSE,zscale=FALSE) ###"M" colour key corresponds to expected direction for acc/lat 

NKIscaledfitsF<-mgcvscalefits(NKI,outcomevars = c("Accuracycomposite","Latencycomposite"),predvars = "metaage",mformula = NKIformula,scale=FALSE)
NKIscaledfitsF$sampletyp<-"full"
NKIderivFacc<-multi_mgcvgam_growthrate_multiplot_rasteronly(df=NKI,outcomevars = c("Accuracycomposite"),predvars='metaage',mformula=NKIformula,datarangegrey=TRUE,derivrangemanual=c(-.2,.2),derivcolourmanual=c("grey33","grey33"),devageguides=FALSE,zscale=FALSE) ###"F" colour key corresponds to expected direction for acc/lat 
NKIderivFlat<-multi_mgcvgam_growthrate_multiplot_rasteronly(df=NKI,outcomevars = c("Latencycomposite"),predvars='metaage',mformula=NKIformula,datarangegrey=TRUE,derivrangemanual=c(-.2,.2),derivcolourmanual=c("grey33","grey33"),devageguides=FALSE,zscale=FALSE) ###"F" colour key corresponds to expected direction for acc/lat 

allNKIfits<-plyr::rbind.fill(NKIscaledfitnoclin,NKIscaledfitsF)
allNKIfits$outcomelabel<-dplyr::if_else(allNKIfits$outcome=="Accuracycomposite","Accuracy","Latency")
ggNKIalllatfitsacc<-ggplot(allNKIfits,aes(x=pred,y=fit,colour=sampletyp,fill=sampletyp))+geom_line(colour="black")+geom_ribbon(aes(ymin=fit-2*se,ymax=fit+2*se,fill=sampletyp),colour="black",alpha=.5)+
  theme(legend.position = "top")+
  facet_grid(cols=vars(outcomelabel),scales="free")+theme(strip.text.y= element_text(angle = 180))+
  scale_x_continuous(limits = c(8, 36),breaks=(c(10,15,20,25,30,35)))+scale_colour_manual(values=c("grey33","#6a50a3"))+scale_fill_manual(values=c("grey33","#6a50a3"))
ggNKIalllatfitsacc<-LNCDR::lunaize(ggNKIalllatfitsacc)+xlab("Age (years)")+ylab("Executive Function\n (z-score)")+theme(legend.position = "none",legend.title = element_blank())+
  theme(strip.text.y= element_text(angle = 360))+theme(strip.background = element_blank())

NKIall<-(ggNKIalllatfitsacc)/(NKIderivFacc$returnplot+NKIderivFlat$returnplot)/(NKIscaledfitnoclinacc$returnplot+NKIscaledfitnoclinlat$returnplot)+patchwork::plot_layout(heights = c(3,1,1))
ggsave(NKIall,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/Sensplots/NKIpsysens.pdf",height=6,width=8)

####PNC####
PNCformula<-as.formula('outcome~s(pred)')
PNCscaledfitnoclin<-mgcvscalefits(PNC[PNC$LTN!=1,],outcomevars = c("Accuracycomposite","Latencycomposite"),predvars = "metaage",mformula = PNCformula,scale=FALSE)
PNCscaledfitnoclin$sampletyp<-"noclin"
PNCscaledfitnoclinacc<-multi_mgcvgam_growthrate_multiplot_rasteronly(df=PNC[PNC$LTN!=1,],outcomevars = c("Accuracycomposite"),predvars='metaage',mformula=PNCformula,datarangegrey=TRUE,derivrangemanual=c(-.2,.2),derivcolourmanual=c("grey33","#6a50a3"),devageguides=FALSE,zscale=FALSE) ###"M" colour key corresponds to expected direction for acc/lat 
PNCscaledfitnoclinlat<-multi_mgcvgam_growthrate_multiplot_rasteronly(df=PNC[PNC$LTN!=1,],outcomevars = c("Latencycomposite"),predvars='metaage',mformula=PNCformula,datarangegrey=TRUE,derivrangemanual=c(-.2,.2),derivcolourmanual=c("#6a50a3","grey33"),devageguides=FALSE,zscale=FALSE) ###"M" colour key corresponds to expected direction for acc/lat 

PNCscaledfitsF<-mgcvscalefits(PNC,outcomevars = c("Accuracycomposite","Latencycomposite"),predvars = "metaage",mformula = PNCformula,scale=FALSE)
PNCscaledfitsF$sampletyp<-"full"
PNCderivFacc<-multi_mgcvgam_growthrate_multiplot_rasteronly(df=PNC,outcomevars = c("Accuracycomposite"),predvars='metaage',mformula=PNCformula,datarangegrey=TRUE,derivrangemanual=c(-.2,.2),derivcolourmanual=c("grey33","grey33"),devageguides=FALSE,zscale=FALSE) ###"F" colour key corresponds to expected direction for acc/lat 
PNCderivFlat<-multi_mgcvgam_growthrate_multiplot_rasteronly(df=PNC,outcomevars = c("Latencycomposite"),predvars='metaage',mformula=PNCformula,datarangegrey=TRUE,derivrangemanual=c(-.2,.2),derivcolourmanual=c("grey33","grey33"),devageguides=FALSE,zscale=FALSE) ###"F" colour key corresponds to expected direction for acc/lat 

allPNCfits<-plyr::rbind.fill(PNCscaledfitnoclin,PNCscaledfitsF)
allPNCfits$outcomelabel<-dplyr::if_else(allPNCfits$outcome=="Accuracycomposite","Accuracy","Latency")
ggPNCalllatfitsacc<-ggplot(allPNCfits,aes(x=pred,y=fit,colour=sampletyp,fill=sampletyp))+geom_line(colour="black")+geom_ribbon(aes(ymin=fit-2*se,ymax=fit+2*se,fill=sampletyp),colour="black",alpha=.5)+
  theme(legend.position = "top")+
  facet_grid(cols=vars(outcomelabel),scales="free")+theme(strip.text.y= element_text(angle = 180))+
  scale_x_continuous(limits = c(8, 36),breaks=(c(10,15,20,25,30,35)))+scale_colour_manual(values=c("grey33","#6a50a3"))+scale_fill_manual(values=c("grey33","#6a50a3"))
ggPNCalllatfitsacc<-LNCDR::lunaize(ggPNCalllatfitsacc)+xlab("Age (years)")+ylab("Executive Function\n (z-score)")+theme(legend.position = "none",legend.title = element_blank())+
  theme(strip.text.y= element_text(angle = 360))+theme(strip.background = element_blank())

PNCall<-(ggPNCalllatfitsacc)/(PNCderivFacc$returnplot+PNCderivFlat$returnplot)/(PNCscaledfitnoclinacc$returnplot+PNCscaledfitnoclinlat$returnplot)+patchwork::plot_layout(heights = c(3,1,1))
ggsave(PNCall,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/Sensplots/PNCpsysens.pdf",height=6,width=8)

#########all panels#########
lunaallwithNCANDAall<-(gglunaalllatfitsacc+ggNCANDAalllatfitsacc)/(lunaderivfullacc$returnplot+lunaderivfulllat$returnplot+NCA$returnplot+NCANDAderivFlat$returnplot)/(lunaderivMacc$returnplot+lunaderivMlat$returnplot+NCANDAderivMacc$returnplot+NCANDAderivMlat$returnplot)+patchwork::plot_layout(heights = c(3,1,1))

alllunafits$dataset<-"Luna"
allNCANDAfits$dataset<-"NCANDA"
allNKIfits$dataset<-"NKI"
allPNCfits$dataset<-"PNC"

outvars<-c("pred","fit","sampletyp","se","outcomelabel","dataset")

supdatafits<-rbind(alllunafits[,outvars],allNCANDAfits[,outvars]) %>% rbind(.,allNKIfits[,outvars]) %>%
  rbind(.,allPNCfits[,outvars])

write.csv(supdatafits,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Data/SupportingData/Sup17.csv")

supdatafits$dataset_sampletyp<-paste(supdatafits$dataset,supdatafits$sampletyp)
ggplot(supdatafits,aes(x=pred,y=fit,colour=sampletyp))+geom_line()+
  geom_ribbon(aes(ymin=fit-2*se,ymax=fit+2*se,fill=sampletyp),colour="black",alpha=.5)+
  facet_grid(rows=vars(dataset_sampletyp),cols=vars(outcomelabel),scales = "free")




