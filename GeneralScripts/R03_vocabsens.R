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
####
NCANDA<-read.csv("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/NCANDA/data/btc_R03cleaneddata_20220306.outlierremoved.compositeacclat.2030scale.csv")
NCANDA$id<-as.factor(NCANDA$subject)
NCANDA$metaage<-NCANDA$cnp_age
NCANDA$visitnum<-as.numeric(as.factor(NCANDA$visit))
NCANDA$dataset<-"NCANDA"
NCANDAparentreport<-read.csv("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/NCANDA_Behavioral/Data/NCANDA_RELEASE_4Y_REDCAP_MEASUREMENTS_V02/summaries/redcap/parentreport.csv")
NCANDAparentreport<-NCANDAparentreport[NCANDAparentreport$visit=="baseline",]
NCANDAparentreport$parentedaverage<-compositecols(c("parentreport_sesp3","parentreport_sesp3b"),data=NCANDAparentreport,mincols=1)##composite but allow for only one parental score
NCANDA<-merge(NCANDA,NCANDAparentreport[,c("subject","parentreport_sesp3","parentreport_sesp3b","parentedaverage")],by=("subject"),all.x=TRUE)
NCANDA$vocab<-NCANDA$cnp_pvoc_pvoccr
NCANDA$PVR_CR<-NCANDA$cnp_pvrt_pvrtcr
cor(NCANDA$vocab,NCANDA$cnp_pvrt_pvrtcr,use="complete")
####
###
NKI<-read.csv("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/NKI/Data/btc_NKIscoredmeasures_20220306.outlierremoved.compositeacclat.2030scale.csv")
NKI$id<-as.factor(NKI$subject)
NKI$metaage<-NKI$age
NKI$dataset<-"NKI"
NKI$visitnum<-1 ##cross-sectional
NKI$PVR_CR<-NKI$penncnp_0195
NKI$vocab<-NKI$int_04 ###using t scores because of upperbound that differs by age for Raw scores
cor(NKI$PVR_CR,NKI$vocab,use="complete")
####
#PVRT_CR#
PNC<-read.csv("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/PNC/Data/btc_R03cleaneddata_20220306.outlierremoved.compositeacclat.2030scale.csv")
PNC$id<-as.factor(PNC$dbGaP_Subject_ID)
PNC$metaage<-PNC$age
PNC$dataset<-"PNC"
PNC$visitnum<-1 ##cross-sectional
PNCPVRT<-read.csv("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/PNC/Data/PNCPVRT.csv")
PNCPVRT$id<-PNCPVRT$dbGaP_Subject_ID
PNCPVRT<-PNCPVRT[PNCPVRT$id %in% PNC$id,]
PNCPVRT$PVRT_CR<-as.numeric(PNCPVRT$PVRT_CR)
PNCPVRT<-PNCPVRT %>% group_by(id) %>% summarize(PVRT_CR=mean(PVRT_CR,na.rm=TRUE)) ##reduce to 1 measurement per participant; mean in case inconsistent report
PNCPVRT$PVR_CR<-PNCPVRT$PVRT_CR ##match name to other datasets
PNC<-merge(PNC,PNCPVRT,by=c("id"),all.x=TRUE)
####################
###note all internal scaling set to false due want to keep on original, full sample z units of composite metrics
###without internal scaling derivative plots will revert to T-scaling; saturating this scale here such that significant always shows same colour to minimize plotting complexity (i.e., is significant basemodel, is significant covmodel )
#####NCANDA######
NCANDAformula<-as.formula('outcome~s(pred)+s(visitnum,k=5)')
NCANDAscaledfitbase<-mgcvscalefits(NCANDA,outcomevars = c("Accuracycomposite","Latencycomposite"),predvars = "metaage",idvar="id",mformula = NCANDAformula,scale=FALSE)
NCANDAscaledfitbase$modeltype<-"base"
NCANDAscaledfitbaseacc<-multi_mgcvgam_growthrate_multiplot_rasteronly(df=NCANDA,outcomevars = c("Accuracycomposite"),idvar="id",predvars='metaage',mformula=NCANDAformula,datarangegrey=TRUE,derivrangemanual=c(-.2,.2),derivcolourmanual=c("grey33","grey33"),devageguides=FALSE,zscale=FALSE)
NCANDAscaledfitbaseacclat<-multi_mgcvgam_growthrate_multiplot_rasteronly(df=NCANDA,outcomevars = c("Latencycomposite"),idvar="id",predvars='metaage',mformula=NCANDAformula,datarangegrey=TRUE,derivrangemanual=c(-.2,.2),derivcolourmanual=c("grey33","grey33"),devageguides=FALSE,zscale=FALSE) 

NCANDAformulacov<-as.formula('outcome~s(pred)+s(visitnum,k=5)+s(PVR_CR)')
NCANDAscaledfitscov<-mgcvscalefits(NCANDA,outcomevars = c("Accuracycomposite","Latencycomposite"),predvars = "metaage",idvar="id",mformula = NCANDAformulacov,scale=FALSE)
NCANDAscaledfitscov$modeltype<-"coved"
NCANDAscaledfitscovacc<-multi_mgcvgam_growthrate_multiplot_rasteronly(df=NCANDA,outcomevars = c("Accuracycomposite"),idvar="id",predvars='metaage',mformula=NCANDAformulacov,datarangegrey=TRUE,derivrangemanual=c(-.2,.2),derivcolourmanual=c("grey33","#238b45"),devageguides=FALSE,zscale=FALSE)###colour key corresponds to expected direction for acc/lat 
NCANDAscaledfitscovlat<-multi_mgcvgam_growthrate_multiplot_rasteronly(df=NCANDA,outcomevars = c("Latencycomposite"),idvar="id",predvars='metaage',mformula=NCANDAformulacov,datarangegrey=TRUE,derivrangemanual=c(-.2,.2),derivcolourmanual=c("#238b45","grey33"),devageguides=FALSE,zscale=FALSE) ###colour key corresponds to expected direction for acc/lat  

NCANDAformulacovvocab<-as.formula('outcome~s(pred)+s(visitnum,k=5)+s(vocab)')
NCANDAscaledfitscovvocab<-mgcvscalefits(NCANDA,outcomevars = c("Accuracycomposite","Latencycomposite"),predvars = "metaage",idvar="id",mformula = NCANDAformulacovvocab,scale=FALSE)
NCANDAscaledfitscovvocab$modeltype<-"covedvocab"
NCANDAscaledfitscovvocabacc<-multi_mgcvgam_growthrate_multiplot_rasteronly(df=NCANDA,outcomevars = c("Accuracycomposite"),idvar="id",predvars='metaage',mformula=NCANDAformulacovvocab,datarangegrey=TRUE,derivrangemanual=c(-.2,.2),derivcolourmanual=c("grey33","#00441b"),devageguides=FALSE,zscale=FALSE)###colour key corresponds to expected direction for acc/lat 
NCANDAscaledfitscovvocablat<-multi_mgcvgam_growthrate_multiplot_rasteronly(df=NCANDA,outcomevars = c("Latencycomposite"),idvar="id",predvars='metaage',mformula=NCANDAformulacovvocab,datarangegrey=TRUE,derivrangemanual=c(-.2,.2),derivcolourmanual=c("#00441b","grey33"),devageguides=FALSE,zscale=FALSE) ###colour key corresponds to expected direction for acc/lat  


allNCANDAfits<-plyr::rbind.fill(NCANDAscaledfitbase,NCANDAscaledfitscov) %>% plyr::rbind.fill(.,NCANDAscaledfitscovvocab)
allNCANDAfits$outcomelabel<-dplyr::if_else(allNCANDAfits$outcome=="Accuracycomposite","Accuracy","Latency")

ggNCANDAalllatfitsacc<-ggplot(allNCANDAfits,aes(x=pred,y=fit,colour=modeltype,fill=modeltype))+geom_line(colour="black")+geom_ribbon(aes(ymin=fit-2*se,ymax=fit+2*se,fill=modeltype),colour="black",alpha=.5)+
  theme(legend.position = "top")+
  facet_grid(cols=vars(outcomelabel),scales="free")+theme(strip.text.y= element_text(angle = 180))+
  scale_x_continuous(limits = c(8, 36),breaks=(c(10,15,20,25,30,35)))+scale_colour_manual(values=c("grey11","#238B45","#00441B"))+scale_fill_manual(values=c("grey11","#238B45","#00441B"))
ggNCANDAalllatfitsacc<-LNCDR::lunaize(ggNCANDAalllatfitsacc)+xlab("Age (years)")+ylab("Executive Function\n (z-score)")+theme(legend.position = "none",legend.title = element_blank())+
  theme(strip.text.y= element_text(angle = 360))+theme(strip.background = element_blank())

NCANDAall<-(ggNCANDAalllatfitsacc)/(NCANDAscaledfitbaseacc$returnplot+NCANDAscaledfitbaseacclat$returnplot)/(NCANDAscaledfitscovacc$returnplot+NCANDAscaledfitscovlat$returnplot)/(NCANDAscaledfitscovvocabacc$returnplot+NCANDAscaledfitscovvocablat$returnplot)+patchwork::plot_layout(heights = c(3,1,1,1))
ggsave(NCANDAall,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/Sensplots/NCANDAvocabsens.pdf",height=6,width=8)
####NKI####
NKIformula<-as.formula('outcome~s(pred)')
NKIscaledfitbase<-mgcvscalefits(NKI,outcomevars = c("Accuracycomposite","Latencycomposite"),predvars = "metaage",mformula = NKIformula,scale=FALSE)
NKIscaledfitbase$modeltype<-"base"
NKIscaledfitbaseacc<-multi_mgcvgam_growthrate_multiplot_rasteronly(df=NKI,outcomevars = c("Accuracycomposite"),predvars='metaage',mformula=NKIformula,datarangegrey=TRUE,derivrangemanual=c(-.2,.2),derivcolourmanual=c("grey33","grey33"),devageguides=FALSE,zscale=FALSE) 
NKIscaledfitbaseacclat<-multi_mgcvgam_growthrate_multiplot_rasteronly(df=NKI,outcomevars = c("Latencycomposite"),predvars='metaage',mformula=NKIformula,datarangegrey=TRUE,derivrangemanual=c(-.2,.2),derivcolourmanual=c("grey33","grey33"),devageguides=FALSE,zscale=FALSE) 

NKIformulacov<-as.formula('outcome~s(pred)+s(PVR_CR)')
NKIscaledfitscov<-mgcvscalefits(NKI,outcomevars = c("Accuracycomposite","Latencycomposite"),predvars = "metaage",mformula = NKIformulacov,scale=FALSE)
NKIscaledfitscov$modeltype<-"coved"
NKIscaledfitscovacc<-multi_mgcvgam_growthrate_multiplot_rasteronly(df=NKI,outcomevars = c("Accuracycomposite"),predvars='metaage',mformula=NKIformulacov,datarangegrey=TRUE,derivrangemanual=c(-.2,.2),derivcolourmanual=c("grey33","#238b45"),devageguides=FALSE,zscale=FALSE) ###colour key corresponds to expected direction for acc/lat 
NKIscaledfitscovlat<-multi_mgcvgam_growthrate_multiplot_rasteronly(df=NKI,outcomevars = c("Latencycomposite"),predvars='metaage',mformula=NKIformulacov,datarangegrey=TRUE,derivrangemanual=c(-.2,.2),derivcolourmanual=c("#238b45","grey33"),devageguides=FALSE,zscale=FALSE) ###colour key corresponds to expected direction for acc/lat 

NIformulacovvocab<-as.formula('outcome~s(pred)+s(vocab)')
NKIscaledfitscovvocab<-mgcvscalefits(NKI,outcomevars = c("Accuracycomposite","Latencycomposite"),predvars = "metaage",mformula = NIformulacovvocab,scale=FALSE)
NKIscaledfitscovvocab$modeltype<-"covedvocab"
NKIscaledfitscovvocabacc<-multi_mgcvgam_growthrate_multiplot_rasteronly(df=NKI,outcomevars = c("Accuracycomposite"),predvars='metaage',mformula=NIformulacovvocab,datarangegrey=TRUE,derivrangemanual=c(-.2,.2),derivcolourmanual=c("grey33","#00441b"),devageguides=FALSE,zscale=FALSE)###colour key corresponds to expected direction for acc/lat 
NKIscaledfitscovvocablat<-multi_mgcvgam_growthrate_multiplot_rasteronly(df=NKI,outcomevars = c("Latencycomposite"),predvars='metaage',mformula=NIformulacovvocab,datarangegrey=TRUE,derivrangemanual=c(-.2,.2),derivcolourmanual=c("#00441b","grey33"),devageguides=FALSE,zscale=FALSE) ###colour key corresponds to expected direction for acc/lat  

allNKIfits<-plyr::rbind.fill(NKIscaledfitbase,NKIscaledfitscov) %>% plyr::rbind.fill(.,NKIscaledfitscovvocab)
allNKIfits$outcomelabel<-dplyr::if_else(allNKIfits$outcome=="Accuracycomposite","Accuracy","Latency")

ggNKIalllatfitsacc<-ggplot(allNKIfits,aes(x=pred,y=fit,colour=modeltype,fill=modeltype))+geom_line(colour="black")+geom_ribbon(aes(ymin=fit-2*se,ymax=fit+2*se,fill=modeltype),colour="black",alpha=.5)+
  theme(legend.position = "top")+
  facet_grid(cols=vars(outcomelabel),scales="free")+theme(strip.text.y= element_text(angle = 180))+
  scale_x_continuous(limits = c(8, 36),breaks=(c(10,15,20,25,30,35)))+scale_colour_manual(values=c("grey11","#238B45","#00441B"))+scale_fill_manual(values=c("grey11","#238B45","#00441B"))
ggNKIalllatfitsacc<-LNCDR::lunaize(ggNKIalllatfitsacc)+xlab("Age (years)")+ylab("Executive Function\n (z-score)")+theme(legend.position = "none",legend.title = element_blank())+
  theme(strip.text.y= element_text(angle = 360))+theme(strip.background = element_blank())

NKIall<-(ggNKIalllatfitsacc)/(NKIscaledfitbaseacc$returnplot+NKIscaledfitbaseacclat$returnplot)/(NKIscaledfitscovacc$returnplot+NKIscaledfitscovlat$returnplot)/(NKIscaledfitscovvocabacc$returnplot+NKIscaledfitscovvocablat$returnplot)+patchwork::plot_layout(heights = c(3,1,1,1))
ggsave(NKIall,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/Sensplots/NKIvocabsens.pdf",height=6,width=8)
####PNC####
PNCformula<-as.formula('outcome~s(pred)')
PNCscaledfitbase<-mgcvscalefits(PNC,outcomevars = c("Accuracycomposite","Latencycomposite"),predvars = "metaage",mformula = PNCformula,scale=FALSE)
PNCscaledfitbase$modeltype<-"base"
PNCscaledfitbaseacc<-multi_mgcvgam_growthrate_multiplot_rasteronly(df=PNC,outcomevars = c("Accuracycomposite"),predvars='metaage',mformula=PNCformula,datarangegrey=TRUE,derivrangemanual=c(-.2,.2),derivcolourmanual=c("grey33","grey33"),devageguides=FALSE,zscale=FALSE) 
PNCscaledfitbaseacclat<-multi_mgcvgam_growthrate_multiplot_rasteronly(df=PNC,outcomevars = c("Latencycomposite"),predvars='metaage',mformula=PNCformula,datarangegrey=TRUE,derivrangemanual=c(-.2,.2),derivcolourmanual=c("grey33","grey33"),devageguides=FALSE,zscale=FALSE) 

PNCformulacov<-as.formula('outcome~s(pred)+s(PVR_CR)')
PNCscaledfitscov<-mgcvscalefits(PNC,outcomevars = c("Accuracycomposite","Latencycomposite"),predvars = "metaage",mformula = PNCformulacov,scale=FALSE)
PNCscaledfitscov$modeltype<-"coved"
PNCscaledfitscovacc<-multi_mgcvgam_growthrate_multiplot_rasteronly(df=PNC,outcomevars = c("Accuracycomposite"),predvars='metaage',mformula=PNCformulacov,datarangegrey=TRUE,derivrangemanual=c(-.2,.2),derivcolourmanual=c("grey33","#238b45"),devageguides=FALSE,zscale=FALSE) ###colour key corresponds to expected direction for acc/lat 
PNCscaledfitscovlat<-multi_mgcvgam_growthrate_multiplot_rasteronly(df=PNC,outcomevars = c("Latencycomposite"),predvars='metaage',mformula=PNCformulacov,datarangegrey=TRUE,derivrangemanual=c(-.2,.2),derivcolourmanual=c("#238b45","grey33"),devageguides=FALSE,zscale=FALSE) ###colour key corresponds to expected direction for acc/lat 

allPNCfits<-plyr::rbind.fill(PNCscaledfitbase,PNCscaledfitscov)
allPNCfits$outcomelabel<-dplyr::if_else(allPNCfits$outcome=="Accuracycomposite","Accuracy","Latency")

ggPNCalllatfitsacc<-ggplot(allPNCfits,aes(x=pred,y=fit,colour=modeltype,fill=modeltype))+geom_line(colour="black")+geom_ribbon(aes(ymin=fit-2*se,ymax=fit+2*se,fill=modeltype),colour="black",alpha=.5)+
  theme(legend.position = "top")+
  facet_grid(cols=vars(outcomelabel),scales="free")+theme(strip.text.y= element_text(angle = 180))+
  scale_x_continuous(limits = c(8, 36),breaks=(c(10,15,20,25,30,35)))+scale_colour_manual(values=c("grey11","#238b45"))+scale_fill_manual(values=c("grey11","#238b45"))
ggPNCalllatfitsacc<-LNCDR::lunaize(ggPNCalllatfitsacc)+xlab("Age (years)")+ylab("Executive Function\n (z-score)")+theme(legend.position = "none",legend.title = element_blank())+
  theme(strip.text.y= element_text(angle = 360))+theme(strip.background = element_blank())

PNCall<-(ggPNCalllatfitsacc)/(PNCscaledfitbaseacc$returnplot+PNCscaledfitbaseacclat$returnplot)/(PNCscaledfitscovacc$returnplot+PNCscaledfitscovlat$returnplot)/(PNCscaledfitscovacc$returnplot+PNCscaledfitscovlat$returnplot)+patchwork::plot_layout(heights = c(3,1,1,1))
ggsave(PNCall,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/Sensplots/PNCvocabsens.pdf",height=6,width=8)
#########all panels#########

