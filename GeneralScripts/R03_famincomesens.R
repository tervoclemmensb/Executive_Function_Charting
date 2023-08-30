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
lunaincome<-read.csv("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_Behavioral/Data/lunaincome.csv")
LUNA<-merge(LUNA,lunaincome,by=c("id"),all.x=TRUE)
LUNA$income_recode_impute<-relevel(factor(LUNA$income_recode_impute),ref="75k-99k")
####
NCANDA<-read.csv("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/NCANDA/data/btc_R03cleaneddata_20220306.outlierremoved.compositeacclat.2030scale.csv")
NCANDA$id<-as.factor(NCANDA$subject)
NCANDA$metaage<-NCANDA$cnp_age
NCANDA$visitnum<-as.numeric(as.factor(NCANDA$visit))
NCANDA$dataset<-"NCANDA"
NCANDAincome<-read.csv("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_Behavioral/NCANDA/Data/ncandaincome.csv")
NCANDA<-merge(NCANDA,NCANDAincome,by=("id"),all.x=TRUE)
NCANDA$income_recode_impute<-relevel(factor(NCANDA$income_recode_impute),ref="75k-99k")
####
NKI<-read.csv("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/NKI/Data/btc_NKIscoredmeasures_20220306.outlierremoved.compositeacclat.2030scale.csv")
NKI$id<-as.factor(NKI$subject)
NKI$metaage<-NKI$age
NKI$dataset<-"NKI"
NKI$visitnum<-1 ##cross-sectional
NKIincome<-read.csv("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_Behavioral/NKI/Data/nkiincome.csv")
NKI<-merge(NKI,NKIincome,by=c("id"),all.x=TRUE)
NKI$income_recode_impute<-relevel(factor(NKI$income_recode_impute),ref="75k-99k")
######
#####
####
###no income in PNC####
####
####################
###note all internal scaling set to false due want to keep on original, full sample z units of composite metrics
###without internal scaling derivative plots will revert to T-scaling; saturating this scale here such that significant always shows same colour to minimize plotting complexity (i.e., is significant basemodel, is significant covmodel )
#######models#######
###LUNA
lunaformula<-as.formula('outcome~s(pred)+s(visitnum)')
lunascaledfitbase<-mgcvscalefits(LUNA,outcomevars = c("Accuracycomposite","Latencycomposite"),predvars = "Ageatvisit",idvar="id",mformula = lunaformula,scale=FALSE)
lunascaledfitbase$modeltype<-"base"
lunascaledfitbaseacc<-multi_mgcvgam_growthrate_multiplot_rasteronly(df=LUNA,outcomevars = c("Accuracycomposite"),idvar="id",predvars='Ageatvisit',mformula=lunaformula,datarangegrey=TRUE,derivrangemanual=c(-.2,.2),derivcolourmanual=c("grey33","grey33"),devageguides=FALSE,zscale=FALSE)
lunascaledfitbaseacclat<-multi_mgcvgam_growthrate_multiplot_rasteronly(df=LUNA,outcomevars = c("Latencycomposite"),idvar="id",predvars='Ageatvisit',mformula=lunaformula,datarangegrey=TRUE,derivrangemanual=c(-.2,.2),derivcolourmanual=c("grey33","grey33"),devageguides=FALSE,zscale=FALSE)  

lunaformulacov<-as.formula('outcome~s(pred)+s(visitnum)+income_recode_impute')
lunascaledfitscov<-mgcvscalefits(LUNA,outcomevars = c("Accuracycomposite","Latencycomposite"),predvars = "Ageatvisit",idvar="id",mformula = lunaformulacov,scale=FALSE)
lunascaledfitscov$modeltype<-"coved"
lunascaledfitscovacc<-multi_mgcvgam_growthrate_multiplot_rasteronly(df=LUNA,outcomevars = c("Accuracycomposite"),idvar="id",predvars='Ageatvisit',mformula=lunaformulacov,datarangegrey=TRUE,derivrangemanual=c(-.2,.2),derivcolourmanual=c("grey33","#cb9c3f"),devageguides=FALSE,zscale=FALSE) ###colour key corresponds to expected direction for acc/lat 
lunascaledfitscovlat<-multi_mgcvgam_growthrate_multiplot_rasteronly(df=LUNA,outcomevars = c("Latencycomposite"),idvar="id",predvars='Ageatvisit',mformula=lunaformulacov,datarangegrey=TRUE,derivrangemanual=c(-.2,.2),derivcolourmanual=c("#cb9c3f","grey33"),devageguides=FALSE,zscale=FALSE)###colour key corresponds to expected direction for acc/lat 

alllunafits<-plyr::rbind.fill(lunascaledfitbase,lunascaledfitscov)
alllunafits$outcomelabel<-dplyr::if_else(alllunafits$outcome=="Accuracycomposite","Accuracy","Latency")

gglunaalllatfitsacc<-ggplot(alllunafits,aes(x=pred,y=fit,colour=modeltype,fill=modeltype))+geom_line(colour="black")+geom_ribbon(aes(ymin=fit-2*se,ymax=fit+2*se,fill=modeltype),colour="black",alpha=.5)+
  theme(legend.position = "top")+
  facet_grid(cols=vars(outcomelabel),scales="free")+theme(strip.text.y= element_text(angle = 180))+
  scale_x_continuous(limits = c(8, 36),breaks=(c(10,15,20,25,30,35)))+scale_colour_manual(values=c("grey11","#cb9c3f"))+scale_fill_manual(values=c("grey11","#cb9c3f"))
gglunaalllatfitsacc<-LNCDR::lunaize(gglunaalllatfitsacc)+xlab("Age (years)")+ylab("Executive Function\n (z-score)")+theme(legend.position = "none",legend.title = element_blank())+
  theme(strip.text.y= element_text(angle = 360))+theme(strip.background = element_blank())

lunaall<-(gglunaalllatfitsacc)/(lunascaledfitbaseacc$returnplot+lunascaledfitbaseacclat$returnplot)/(lunascaledfitscovacc$returnplot+lunascaledfitscovlat$returnplot)+patchwork::plot_layout(heights = c(3,1,1))
ggsave(lunaall,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/Sensplots/lunaincomesens.pdf",height=6,width=8)
#####NCANDA######
NCANDAformula<-as.formula('outcome~s(pred)+s(visitnum,k=5)')
NCANDAscaledfitbase<-mgcvscalefits(NCANDA,outcomevars = c("Accuracycomposite","Latencycomposite"),predvars = "metaage",idvar="id",mformula = NCANDAformula,scale=FALSE)
NCANDAscaledfitbase$modeltype<-"base"
NCANDAscaledfitbaseacc<-multi_mgcvgam_growthrate_multiplot_rasteronly(df=NCANDA,outcomevars = c("Accuracycomposite"),idvar="id",predvars='metaage',mformula=NCANDAformula,datarangegrey=TRUE,derivrangemanual=c(-.2,.2),derivcolourmanual=c("grey33","grey33"),devageguides=FALSE,zscale=FALSE)
NCANDAscaledfitbaseacclat<-multi_mgcvgam_growthrate_multiplot_rasteronly(df=NCANDA,outcomevars = c("Latencycomposite"),idvar="id",predvars='metaage',mformula=NCANDAformula,datarangegrey=TRUE,derivrangemanual=c(-.2,.2),derivcolourmanual=c("grey33","grey33"),devageguides=FALSE,zscale=FALSE) 

NCANDAformulacov<-as.formula('outcome~s(pred)+s(visitnum,k=5)+income_recode_impute')
NCANDAscaledfitscov<-mgcvscalefits(NCANDA,outcomevars = c("Accuracycomposite","Latencycomposite"),predvars = "metaage",idvar="id",mformula = NCANDAformulacov,scale=FALSE)
NCANDAscaledfitscov$modeltype<-"coved"
NCANDAscaledfitscovacc<-multi_mgcvgam_growthrate_multiplot_rasteronly(df=NCANDA,outcomevars = c("Accuracycomposite"),idvar="id",predvars='metaage',mformula=NCANDAformulacov,datarangegrey=TRUE,derivrangemanual=c(-.2,.2),derivcolourmanual=c("grey33","#cb9c3f"),devageguides=FALSE,zscale=FALSE)###colour key corresponds to expected direction for acc/lat 
NCANDAscaledfitscovlat<-multi_mgcvgam_growthrate_multiplot_rasteronly(df=NCANDA,outcomevars = c("Latencycomposite"),idvar="id",predvars='metaage',mformula=NCANDAformulacov,datarangegrey=TRUE,derivrangemanual=c(-.2,.2),derivcolourmanual=c("#cb9c3f","grey33"),devageguides=FALSE,zscale=FALSE) ###colour key corresponds to expected direction for acc/lat  

allNCANDAfits<-plyr::rbind.fill(NCANDAscaledfitbase,NCANDAscaledfitscov)
allNCANDAfits$outcomelabel<-dplyr::if_else(allNCANDAfits$outcome=="Accuracycomposite","Accuracy","Latency")

ggNCANDAalllatfitsacc<-ggplot(allNCANDAfits,aes(x=pred,y=fit,colour=modeltype,fill=modeltype))+geom_line(colour="black")+geom_ribbon(aes(ymin=fit-2*se,ymax=fit+2*se,fill=modeltype),colour="black",alpha=.5)+
  theme(legend.position = "top")+
  facet_grid(cols=vars(outcomelabel),scales="free")+theme(strip.text.y= element_text(angle = 180))+
  scale_x_continuous(limits = c(8, 36),breaks=(c(10,15,20,25,30,35)))+scale_colour_manual(values=c("grey11","#cb9c3f"))+scale_fill_manual(values=c("grey11","#cb9c3f"))
ggNCANDAalllatfitsacc<-LNCDR::lunaize(ggNCANDAalllatfitsacc)+xlab("Age (years)")+ylab("Executive Function\n (z-score)")+theme(legend.position = "none",legend.title = element_blank())+
  theme(strip.text.y= element_text(angle = 360))+theme(strip.background = element_blank())

NCANDAall<-(ggNCANDAalllatfitsacc)/(NCANDAscaledfitbaseacc$returnplot+NCANDAscaledfitbaseacclat$returnplot)/(NCANDAscaledfitscovacc$returnplot+NCANDAscaledfitscovlat$returnplot)+patchwork::plot_layout(heights = c(3,1,1))
ggsave(NCANDAall,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/Sensplots/NCANDAincomesens.pdf",height=6,width=8)
####NKI####
NKIformula<-as.formula('outcome~s(pred)')
NKIscaledfitbase<-mgcvscalefits(NKI,outcomevars = c("Accuracycomposite","Latencycomposite"),predvars = "metaage",mformula = NKIformula,scale=FALSE)
NKIscaledfitbase$modeltype<-"base"
NKIscaledfitbaseacc<-multi_mgcvgam_growthrate_multiplot_rasteronly(df=NKI,outcomevars = c("Accuracycomposite"),predvars='metaage',mformula=NKIformula,datarangegrey=TRUE,derivrangemanual=c(-.2,.2),derivcolourmanual=c("grey33","grey33"),devageguides=FALSE,zscale=FALSE) 
NKIscaledfitbaseacclat<-multi_mgcvgam_growthrate_multiplot_rasteronly(df=NKI,outcomevars = c("Latencycomposite"),predvars='metaage',mformula=NKIformula,datarangegrey=TRUE,derivrangemanual=c(-.2,.2),derivcolourmanual=c("grey33","grey33"),devageguides=FALSE,zscale=FALSE) 

NKIformulacov<-as.formula('outcome~s(pred)+income_recode_impute')
NKIscaledfitscov<-mgcvscalefits(NKI,outcomevars = c("Accuracycomposite","Latencycomposite"),predvars = "metaage",mformula = NKIformulacov,scale=FALSE)
NKIscaledfitscov$modeltype<-"coved"
NKIscaledfitscovacc<-multi_mgcvgam_growthrate_multiplot_rasteronly(df=NKI,outcomevars = c("Accuracycomposite"),predvars='metaage',mformula=NKIformulacov,datarangegrey=TRUE,derivrangemanual=c(-.2,.2),derivcolourmanual=c("grey33","#cb9c3f"),devageguides=FALSE,zscale=FALSE) ###colour key corresponds to expected direction for acc/lat 
NKIscaledfitscovlat<-multi_mgcvgam_growthrate_multiplot_rasteronly(df=NKI,outcomevars = c("Latencycomposite"),predvars='metaage',mformula=NKIformulacov,datarangegrey=TRUE,derivrangemanual=c(-.2,.2),derivcolourmanual=c("#cb9c3f","grey33"),devageguides=FALSE,zscale=FALSE) ###colour key corresponds to expected direction for acc/lat 

allNKIfits<-plyr::rbind.fill(NKIscaledfitbase,NKIscaledfitscov)
allNKIfits$outcomelabel<-dplyr::if_else(allNKIfits$outcome=="Accuracycomposite","Accuracy","Latency")

ggNKIalllatfitsacc<-ggplot(allNKIfits,aes(x=pred,y=fit,colour=modeltype,fill=modeltype))+geom_line(colour="black")+geom_ribbon(aes(ymin=fit-2*se,ymax=fit+2*se,fill=modeltype),colour="black",alpha=.5)+
  theme(legend.position = "top")+
  facet_grid(cols=vars(outcomelabel),scales="free")+theme(strip.text.y= element_text(angle = 180))+
  scale_x_continuous(limits = c(8, 36),breaks=(c(10,15,20,25,30,35)))+scale_colour_manual(values=c("grey11","#cb9c3f"))+scale_fill_manual(values=c("grey11","#cb9c3f"))
ggNKIalllatfitsacc<-LNCDR::lunaize(ggNKIalllatfitsacc)+xlab("Age (years)")+ylab("Executive Function\n (z-score)")+theme(legend.position = "none",legend.title = element_blank())+
  theme(strip.text.y= element_text(angle = 360))+theme(strip.background = element_blank())

NKIall<-(ggNKIalllatfitsacc)/(NKIscaledfitbaseacc$returnplot+NKIscaledfitbaseacclat$returnplot)/(NKIscaledfitscovacc$returnplot+NKIscaledfitscovlat$returnplot)+patchwork::plot_layout(heights = c(3,1,1))
ggsave(NKIall,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/Sensplots/NKIincomesens.pdf",height=6,width=8)
#########all panels#########
###supporting data
alllunafits$dataset<-"Luna"
allNCANDAfits$dataset<-"NCANDA"
allNKIfits$dataset<-"NKI"

outvars<-c("pred","fit","modeltype","se","outcomelabel","dataset")

supdatafits<-rbind(alllunafits[,outvars],allNCANDAfits[,outvars]) %>% rbind(.,allNKIfits[,outvars])

write.csv(supdatafits,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Data/SupportingData/Sup15.income.csv")

ggplot(supdatafits,aes(x=pred,y=fit,colour=modeltype))+geom_line()+
  geom_ribbon(aes(ymin=fit-2*se,ymax=fit+2*se,fill=modeltype),colour="black",alpha=.5)+
  facet_grid(rows=vars(dataset),cols=vars(outcomelabel))


