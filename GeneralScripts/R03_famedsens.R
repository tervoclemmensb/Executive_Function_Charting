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
Lunaed<-read.table("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/data/cog_parent_edu.tsv",sep="\t")
names(Lunaed)<-c("id","surveydate","maternaled","paternaled")
Lunaed[,c("maternaled","paternaled")]<-lapply(Lunaed[,c("maternaled","paternaled")],function(x){as.numeric(x)})
Lunaed$maternaled[Lunaed$maternaled %in% c(-8,-9)]<-NA
Lunaed$paternaled[Lunaed$paternaled %in% c(-8,-9,9)]<-NA ###-8,-9 N/A values, 9 miscoded
Lunaednona<-Lunaed[rowSums(is.na(Lunaed[,c("maternaled","paternaled")])) != ncol(Lunaed[,c("maternaled","paternaled")]), ]
##take first visit info
Lunaed<-Lunaed %>% dplyr::group_by(id) %>% mutate(visitrank=rank(surveydate))
Lunaedunique<-Lunaed[Lunaed$visitrank==1,] %>% dplyr::group_by(id) %>% dplyr::summarize(maternaled=unique(maternaled),paternaled=unique(paternaled))
LunaeduniquefromLunacog<-Lunaedunique[Lunaedunique$id %in% LUNA$id,]
LunaeduniquefromLunacognona<-data.frame(LunaeduniquefromLunacog[rowSums(is.na(LunaeduniquefromLunacog[,c("maternaled","paternaled")])) != ncol(LunaeduniquefromLunacog[,c("maternaled","paternaled")]), ])
LunaeduniquefromLunacognona$parentedaverage<-compositecols(c("maternaled","paternaled"),data=LunaeduniquefromLunacognona,mincols=1)##composite but allow for only one parental score
LUNA<-merge(LUNA,LunaeduniquefromLunacognona,by=c("id"),all.x=TRUE)
####for table#####
LUNAtable<-LUNA[LUNA$visitnum==1,]
LUNAtable$guardedfact<-NA
LUNAtable$guardedfact[LUNAtable$maternaled %in% c(3,2,1)]<-"Incomplete Highschool"
LUNAtable$guardedfact[LUNAtable$maternaled %in% c(4)]<-"Highschool"
LUNAtable$guardedfact[LUNAtable$maternaled %in% c(5)]<-"Some College"
LUNAtable$guardedfact[LUNAtable$maternaled %in% c(6)]<-"College"
LUNAtable$guardedfact[LUNAtable$maternaled %in% c(7)]<-"Postgraduate"
LUNAtable$guardedfact<-factor(LUNAtable$guardedfact,levels=c("Incomplete Highschool","Highschool","Some College","College","Postgraduate"))
table(LUNAtable$guardedfact)/nrow(LUNAtable)
prop.table(table(LUNAtable$guardedfact))
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
#####
NCANDAtable<-NCANDA[NCANDA$visitnum==1,]
NCANDAtable$guardedfact<-NA
NCANDAtable$guardedfact[NCANDAtable$parentreport_sesp3 < 12 ]<-"Incomplete Highschool"
NCANDAtable$guardedfact[NCANDAtable$parentreport_sesp3==12]<-"Highschool"
NCANDAtable$guardedfact[NCANDAtable$parentreport_sesp3> 12 & NCANDAtable$parentreport_sesp3 < 16]<-"Some College"
NCANDAtable$guardedfact[NCANDAtable$parentreport_sesp3==16]<-"College"
NCANDAtable$guardedfact[NCANDAtable$parentreport_sesp3>16]<-"Postgraduate"
NCANDAtable$guardedfact<-factor(NCANDAtable$guardedfact,levels=c("Incomplete Highschool","Highschool","Some College","College","Postgraduate"))
table(NCANDAtable$guardedfact)/nrow(NCANDAtable)
prop.table(table(NCANDAtable$guardedfact))
####
NKI<-read.csv("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/NKI/Data/btc_NKIscoredmeasures_20220306.outlierremoved.compositeacclat.2030scale.csv")
NKI$id<-as.factor(NKI$subject)
NKI$metaage<-NKI$age
NKI$dataset<-"NKI"
NKI$visitnum<-1 ##cross-sectional
NKIedadult<-read.csv("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/NKI/data/SESNKIdata-2022-02-13T23_09_33.328Z.csv",skip=1)
NKIedadult$subject<-unlist(lapply(strsplit(NKIedadult$Identifiers,"[,]"),"[[",1))
NKIedadult<-NKIedadult[NKIedadult$subject %in% NKI$subject,]
NKIedadult$maternaled_guard1<-NKIedadult$nkises_05b
NKIedadult$guard2ed<-NKIedadult$nkises_07b
NKIedadultsum <-NKIedadult %>% group_by(subject) %>% summarize(maternaled_guard1=unique(maternaled_guard1),guard2ed=unique(guard2ed)) ##reduce to 1 measurement per participant
NKIedadultsum[,c("maternaled_guard1","guard2ed")]<-lapply(NKIedadultsum[,c("maternaled_guard1","guard2ed")],function(x){as.numeric(x)}) ##convert to numeric

NKIedchild<-read.csv("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/NKI/data/SENSNKICHILD-2022-03-07T18_03_11.143Z.csv")
NKIedchild$subject<-unlist(lapply(strsplit(NKIedchild$Identifiers,"[,]"),"[[",1))
NKIedchild<-NKIedchild[NKIedchild$subject %in% NKI$subject,]
NKIedchild$maternaled_guard1<-NKIedchild$sesc_02
NKIedchild$guard2ed<-NKIedchild$sesc_06
NKIedchildsum <-NKIedchild %>% group_by(subject) %>% summarize(maternaled_guard1=unique(maternaled_guard1),guard2ed=unique(guard2ed)) ##reduce to 1 measurement per participant
NKIedchildsum[,c("maternaled_guard1","guard2ed")]<-lapply(NKIedchildsum[,c("maternaled_guard1","guard2ed")],function(x){as.numeric(x)})

allNKIed<-rbind(NKIedadultsum,NKIedchildsum)
allNKIedsum <-allNKIed %>% group_by(subject) %>% summarize(maternaled_guard1=mean(maternaled_guard1,na.rm=TRUE),maternaled_guard1=mean(guard2ed,na.rm=TRUE)) ##reduce to 1 measurement per participant; mean in case inconsistent report

allNKIedsum$parentedaverage<-compositecols(c("maternaled_guard1","maternaled_guard1"),data=allNKIedsum,mincols=1)##composite but allow for only one parental score

NKI<-merge(NKI,allNKIedsum,by=c("subject"),all.x=TRUE)
######
#####
NKI$guardedfact<-NA
NKI$guardedfact[NKI$maternaled_guard1 < 12 ]<-"Incomplete Highschool"
NKI$guardedfact[NKI$maternaled_guard1==12]<-"Highschool"
NKI$guardedfact[NKI$maternaled_guard1> 12 & NKI$maternaled_guard1 < 16]<-"Some College"
NKI$guardedfact[NKI$maternaled_guard1==16]<-"College"
NKI$guardedfact[NKI$maternaled_guard1>16]<-"Postgraduate"
NKI$guardedfact<-factor(NKI$guardedfact,levels=c("Incomplete Highschool","Highschool","Some College","College","Postgraduate"))
table(NKI$guardedfact)/nrow(NKI)
prop.table(table(NKI$guardedfact))
####
PNC<-read.csv("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/PNC/Data/btc_R03cleaneddata_20220306.outlierremoved.compositeacclat.2030scale.csv")
PNC$id<-as.factor(PNC$dbGaP_Subject_ID)
PNC$metaage<-PNC$age
PNC$dataset<-"PNC"
PNC$visitnum<-1 ##cross-sectional
PNCed<-read.csv("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/PNC/Data/PNCparentaled.csv")
PNCed$id<-PNCed$dbGaP_Subject_ID
PNCed<-PNCed[PNCed$id %in% PNC$id,]

PNCedsum <-PNCed %>% group_by(id) %>% summarize(Mother_Education=mean(Mother_Education,na.rm=TRUE),Father_Education=mean(Father_Education,na.rm=TRUE)) ##reduce to 1 measurement per participant; mean in case inconsistent report

PNCedsum$parentedaverage<-compositecols(c("Mother_Education","Father_Education"),data=PNCedsum,mincols=1)
PNC<-merge(PNC,PNCedsum,by=c("id"),all.x=TRUE)

####
PNC$guardedfact<-NA
PNC$guardedfact[PNC$Mother_Education < 12 ]<-"Incomplete Highschool"
PNC$guardedfact[PNC$Mother_Education==12]<-"Highschool"
PNC$guardedfact[PNC$Mother_Education> 12 & PNC$Mother_Education < 16]<-"Some College"
PNC$guardedfact[PNC$Mother_Education==16]<-"College"
PNC$guardedfact[PNC$Mother_Education>16]<-"Postgraduate"
PNC$guardedfact<-factor(PNC$guardedfact,levels=c("Incomplete Highschool","Highschool","Some College","College","Postgraduate"))
table(PNC$guardedfact)/nrow(PNC)
t(prop.table(table(PNC$guardedfact)))
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

lunaformulacov<-as.formula('outcome~s(pred)+s(visitnum)+s(parentedaverage)')
lunascaledfitscov<-mgcvscalefits(LUNA,outcomevars = c("Accuracycomposite","Latencycomposite"),predvars = "Ageatvisit",idvar="id",mformula = lunaformulacov,scale=FALSE)
lunascaledfitscov$modeltype<-"coved"
lunascaledfitscovacc<-multi_mgcvgam_growthrate_multiplot_rasteronly(df=LUNA,outcomevars = c("Accuracycomposite"),idvar="id",predvars='Ageatvisit',mformula=lunaformulacov,datarangegrey=TRUE,derivrangemanual=c(-.2,.2),derivcolourmanual=c("grey33","#fec44f"),devageguides=FALSE,zscale=FALSE) ###colour key corresponds to expected direction for acc/lat 
lunascaledfitscovlat<-multi_mgcvgam_growthrate_multiplot_rasteronly(df=LUNA,outcomevars = c("Latencycomposite"),idvar="id",predvars='Ageatvisit',mformula=lunaformulacov,datarangegrey=TRUE,derivrangemanual=c(-.2,.2),derivcolourmanual=c("#fec44f","grey33"),devageguides=FALSE,zscale=FALSE)###colour key corresponds to expected direction for acc/lat 

alllunafits<-plyr::rbind.fill(lunascaledfitbase,lunascaledfitscov)
alllunafits$outcomelabel<-dplyr::if_else(alllunafits$outcome=="Accuracycomposite","Accuracy","Latency")

gglunaalllatfitsacc<-ggplot(alllunafits,aes(x=pred,y=fit,colour=modeltype,fill=modeltype))+geom_line(colour="black")+geom_ribbon(aes(ymin=fit-2*se,ymax=fit+2*se,fill=modeltype),colour="black",alpha=.5)+
  theme(legend.position = "top")+
  facet_grid(cols=vars(outcomelabel),scales="free")+theme(strip.text.y= element_text(angle = 180))+
  scale_x_continuous(limits = c(8, 36),breaks=(c(10,15,20,25,30,35)))+scale_colour_manual(values=c("grey11","#fec44f"))+scale_fill_manual(values=c("grey11","#fec44f"))
gglunaalllatfitsacc<-LNCDR::lunaize(gglunaalllatfitsacc)+xlab("Age (years)")+ylab("Executive Function\n (z-score)")+theme(legend.position = "none",legend.title = element_blank())+
  theme(strip.text.y= element_text(angle = 360))+theme(strip.background = element_blank())

lunaall<-(gglunaalllatfitsacc)/(lunascaledfitbaseacc$returnplot+lunascaledfitbaseacclat$returnplot)/(lunascaledfitscovacc$returnplot+lunascaledfitscovlat$returnplot)+patchwork::plot_layout(heights = c(3,1,1))
ggsave(lunaall,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/Sensplots/lunaedsens.pdf",height=6,width=8)
#####NCANDA######
NCANDAformula<-as.formula('outcome~s(pred)+s(visitnum,k=5)')
NCANDAscaledfitbase<-mgcvscalefits(NCANDA,outcomevars = c("Accuracycomposite","Latencycomposite"),predvars = "metaage",idvar="id",mformula = NCANDAformula,scale=FALSE)
NCANDAscaledfitbase$modeltype<-"base"
NCANDAscaledfitbaseacc<-multi_mgcvgam_growthrate_multiplot_rasteronly(df=NCANDA,outcomevars = c("Accuracycomposite"),idvar="id",predvars='metaage',mformula=NCANDAformula,datarangegrey=TRUE,derivrangemanual=c(-.2,.2),derivcolourmanual=c("grey33","grey33"),devageguides=FALSE,zscale=FALSE)
NCANDAscaledfitbaseacclat<-multi_mgcvgam_growthrate_multiplot_rasteronly(df=NCANDA,outcomevars = c("Latencycomposite"),idvar="id",predvars='metaage',mformula=NCANDAformula,datarangegrey=TRUE,derivrangemanual=c(-.2,.2),derivcolourmanual=c("grey33","grey33"),devageguides=FALSE,zscale=FALSE) 

NCANDAformulacov<-as.formula('outcome~s(pred)+s(visitnum,k=5)+s(parentedaverage)')
NCANDAscaledfitscov<-mgcvscalefits(NCANDA,outcomevars = c("Accuracycomposite","Latencycomposite"),predvars = "metaage",idvar="id",mformula = NCANDAformulacov,scale=FALSE)
NCANDAscaledfitscov$modeltype<-"coved"
NCANDAscaledfitscovacc<-multi_mgcvgam_growthrate_multiplot_rasteronly(df=NCANDA,outcomevars = c("Accuracycomposite"),idvar="id",predvars='metaage',mformula=NCANDAformulacov,datarangegrey=TRUE,derivrangemanual=c(-.2,.2),derivcolourmanual=c("grey33","#fec44f"),devageguides=FALSE,zscale=FALSE)###colour key corresponds to expected direction for acc/lat 
NCANDAscaledfitscovlat<-multi_mgcvgam_growthrate_multiplot_rasteronly(df=NCANDA,outcomevars = c("Latencycomposite"),idvar="id",predvars='metaage',mformula=NCANDAformulacov,datarangegrey=TRUE,derivrangemanual=c(-.2,.2),derivcolourmanual=c("#fec44f","grey33"),devageguides=FALSE,zscale=FALSE) ###colour key corresponds to expected direction for acc/lat  

allNCANDAfits<-plyr::rbind.fill(NCANDAscaledfitbase,NCANDAscaledfitscov)
allNCANDAfits$outcomelabel<-dplyr::if_else(allNCANDAfits$outcome=="Accuracycomposite","Accuracy","Latency")

ggNCANDAalllatfitsacc<-ggplot(allNCANDAfits,aes(x=pred,y=fit,colour=modeltype,fill=modeltype))+geom_line(colour="black")+geom_ribbon(aes(ymin=fit-2*se,ymax=fit+2*se,fill=modeltype),colour="black",alpha=.5)+
  theme(legend.position = "top")+
  facet_grid(cols=vars(outcomelabel),scales="free")+theme(strip.text.y= element_text(angle = 180))+
  scale_x_continuous(limits = c(8, 36),breaks=(c(10,15,20,25,30,35)))+scale_colour_manual(values=c("grey11","#fec44f"))+scale_fill_manual(values=c("grey11","#fec44f"))
ggNCANDAalllatfitsacc<-LNCDR::lunaize(ggNCANDAalllatfitsacc)+xlab("Age (years)")+ylab("Executive Function\n (z-score)")+theme(legend.position = "none",legend.title = element_blank())+
  theme(strip.text.y= element_text(angle = 360))+theme(strip.background = element_blank())

NCANDAall<-(ggNCANDAalllatfitsacc)/(NCANDAscaledfitbaseacc$returnplot+NCANDAscaledfitbaseacclat$returnplot)/(NCANDAscaledfitscovacc$returnplot+NCANDAscaledfitscovlat$returnplot)+patchwork::plot_layout(heights = c(3,1,1))
ggsave(NCANDAall,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/Sensplots/NCANDAedsens.pdf",height=6,width=8)
####NKI####
NKIformula<-as.formula('outcome~s(pred)')
NKIscaledfitbase<-mgcvscalefits(NKI,outcomevars = c("Accuracycomposite","Latencycomposite"),predvars = "metaage",mformula = NKIformula,scale=FALSE)
NKIscaledfitbase$modeltype<-"base"
NKIscaledfitbaseacc<-multi_mgcvgam_growthrate_multiplot_rasteronly(df=NKI,outcomevars = c("Accuracycomposite"),predvars='metaage',mformula=NKIformula,datarangegrey=TRUE,derivrangemanual=c(-.2,.2),derivcolourmanual=c("grey33","grey33"),devageguides=FALSE,zscale=FALSE) 
NKIscaledfitbaseacclat<-multi_mgcvgam_growthrate_multiplot_rasteronly(df=NKI,outcomevars = c("Latencycomposite"),predvars='metaage',mformula=NKIformula,datarangegrey=TRUE,derivrangemanual=c(-.2,.2),derivcolourmanual=c("grey33","grey33"),devageguides=FALSE,zscale=FALSE) 

NKIformulacov<-as.formula('outcome~s(pred)+s(parentedaverage)')
NKIscaledfitscov<-mgcvscalefits(NKI,outcomevars = c("Accuracycomposite","Latencycomposite"),predvars = "metaage",mformula = NKIformulacov,scale=FALSE)
NKIscaledfitscov$modeltype<-"coved"
NKIscaledfitscovacc<-multi_mgcvgam_growthrate_multiplot_rasteronly(df=NKI,outcomevars = c("Accuracycomposite"),predvars='metaage',mformula=NKIformulacov,datarangegrey=TRUE,derivrangemanual=c(-.2,.2),derivcolourmanual=c("grey33","#fec44f"),devageguides=FALSE,zscale=FALSE) ###colour key corresponds to expected direction for acc/lat 
NKIscaledfitscovlat<-multi_mgcvgam_growthrate_multiplot_rasteronly(df=NKI,outcomevars = c("Latencycomposite"),predvars='metaage',mformula=NKIformulacov,datarangegrey=TRUE,derivrangemanual=c(-.2,.2),derivcolourmanual=c("#fec44f","grey33"),devageguides=FALSE,zscale=FALSE) ###colour key corresponds to expected direction for acc/lat 

allNKIfits<-plyr::rbind.fill(NKIscaledfitbase,NKIscaledfitscov)
allNKIfits$outcomelabel<-dplyr::if_else(allNKIfits$outcome=="Accuracycomposite","Accuracy","Latency")

ggNKIalllatfitsacc<-ggplot(allNKIfits,aes(x=pred,y=fit,colour=modeltype,fill=modeltype))+geom_line(colour="black")+geom_ribbon(aes(ymin=fit-2*se,ymax=fit+2*se,fill=modeltype),colour="black",alpha=.5)+
  theme(legend.position = "top")+
  facet_grid(cols=vars(outcomelabel),scales="free")+theme(strip.text.y= element_text(angle = 180))+
  scale_x_continuous(limits = c(8, 36),breaks=(c(10,15,20,25,30,35)))+scale_colour_manual(values=c("grey11","#fec44f"))+scale_fill_manual(values=c("grey11","#fec44f"))
ggNKIalllatfitsacc<-LNCDR::lunaize(ggNKIalllatfitsacc)+xlab("Age (years)")+ylab("Executive Function\n (z-score)")+theme(legend.position = "none",legend.title = element_blank())+
  theme(strip.text.y= element_text(angle = 360))+theme(strip.background = element_blank())

NKIall<-(ggNKIalllatfitsacc)/(NKIscaledfitbaseacc$returnplot+NKIscaledfitbaseacclat$returnplot)/(NKIscaledfitscovacc$returnplot+NKIscaledfitscovlat$returnplot)+patchwork::plot_layout(heights = c(3,1,1))
ggsave(NKIall,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/Sensplots/NKIedsens.pdf",height=6,width=8)
####PNC####
PNCformula<-as.formula('outcome~s(pred)')
PNCscaledfitbase<-mgcvscalefits(PNC,outcomevars = c("Accuracycomposite","Latencycomposite"),predvars = "metaage",mformula = PNCformula,scale=FALSE)
PNCscaledfitbase$modeltype<-"base"
PNCscaledfitbaseacc<-multi_mgcvgam_growthrate_multiplot_rasteronly(df=PNC,outcomevars = c("Accuracycomposite"),predvars='metaage',mformula=PNCformula,datarangegrey=TRUE,derivrangemanual=c(-.2,.2),derivcolourmanual=c("grey33","grey33"),devageguides=FALSE,zscale=FALSE) 
PNCscaledfitbaseacclat<-multi_mgcvgam_growthrate_multiplot_rasteronly(df=PNC,outcomevars = c("Latencycomposite"),predvars='metaage',mformula=PNCformula,datarangegrey=TRUE,derivrangemanual=c(-.2,.2),derivcolourmanual=c("grey33","grey33"),devageguides=FALSE,zscale=FALSE) 

PNCformulacov<-as.formula('outcome~s(pred)+s(parentedaverage)')
PNCscaledfitscov<-mgcvscalefits(PNC,outcomevars = c("Accuracycomposite","Latencycomposite"),predvars = "metaage",mformula = PNCformulacov,scale=FALSE)
PNCscaledfitscov$modeltype<-"coved"
PNCscaledfitscovacc<-multi_mgcvgam_growthrate_multiplot_rasteronly(df=PNC,outcomevars = c("Accuracycomposite"),predvars='metaage',mformula=PNCformulacov,datarangegrey=TRUE,derivrangemanual=c(-.2,.2),derivcolourmanual=c("grey33","#fec44f"),devageguides=FALSE,zscale=FALSE) ###colour key corresponds to expected direction for acc/lat 
PNCscaledfitscovlat<-multi_mgcvgam_growthrate_multiplot_rasteronly(df=PNC,outcomevars = c("Latencycomposite"),predvars='metaage',mformula=PNCformulacov,datarangegrey=TRUE,derivrangemanual=c(-.2,.2),derivcolourmanual=c("#fec44f","grey33"),devageguides=FALSE,zscale=FALSE) ###colour key corresponds to expected direction for acc/lat 

allPNCfits<-plyr::rbind.fill(PNCscaledfitbase,PNCscaledfitscov)
allPNCfits$outcomelabel<-dplyr::if_else(allPNCfits$outcome=="Accuracycomposite","Accuracy","Latency")

ggPNCalllatfitsacc<-ggplot(allPNCfits,aes(x=pred,y=fit,colour=modeltype,fill=modeltype))+geom_line(colour="black")+geom_ribbon(aes(ymin=fit-2*se,ymax=fit+2*se,fill=modeltype),colour="black",alpha=.5)+
  theme(legend.position = "top")+
  facet_grid(cols=vars(outcomelabel),scales="free")+theme(strip.text.y= element_text(angle = 180))+
  scale_x_continuous(limits = c(8, 36),breaks=(c(10,15,20,25,30,35)))+scale_colour_manual(values=c("grey11","#fec44f"))+scale_fill_manual(values=c("grey11","#fec44f"))
ggPNCalllatfitsacc<-LNCDR::lunaize(ggPNCalllatfitsacc)+xlab("Age (years)")+ylab("Executive Function\n (z-score)")+theme(legend.position = "none",legend.title = element_blank())+
  theme(strip.text.y= element_text(angle = 360))+theme(strip.background = element_blank())

PNCall<-(ggPNCalllatfitsacc)/(PNCscaledfitbaseacc$returnplot+PNCscaledfitbaseacclat$returnplot)/(PNCscaledfitscovacc$returnplot+PNCscaledfitscovlat$returnplot)+patchwork::plot_layout(heights = c(3,1,1))
ggsave(PNCall,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/Sensplots/PNCedsens.pdf",height=6,width=8)
#########all panels#########

alllunafits$dataset<-"Luna"
allNCANDAfits$dataset<-"NCANDA"
allNKIfits$dataset<-"NKI"
allPNCfits$dataset<-"PNC"

outvars<-c("pred","fit","modeltype","se","outcomelabel","dataset")

supdatafits<-rbind(alllunafits[,outvars],allNCANDAfits[,outvars]) %>% rbind(.,allNKIfits[,outvars]) %>% rbind(.,allPNCfits[,outvars])

write.csv(supdatafits,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Data/SupportingData/Sup14.famed.csv")

ggplot(supdatafits,aes(x=pred,y=fit,colour=modeltype))+geom_line()+
  geom_ribbon(aes(ymin=fit-2*se,ymax=fit+2*se,fill=modeltype),colour="black",alpha=.5)+
  facet_grid(rows=vars(dataset),cols=vars(outcomelabel))
