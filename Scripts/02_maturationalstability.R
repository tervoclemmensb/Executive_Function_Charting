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
library(grid)
######functions#####
source("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/GeneralScripts/growthrate_mgcvgam.R")
######functions#####
############
coglongdata<-read.csv("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/data/btc_R03cleaneddata_20220306.outlierremoved.compositeacclat.csv")
coglongdata$id<-as.factor(coglongdata$id)
####vars############
####################
allefvars<-c("Anti_CRLat","Anti_CRR","Mix_CRLat","Mix_CRR","SOC.Problems.solved.in.minimum.moves",
             "SOC.Overallmeaninitialthinkingtime","DMS.PC","DMS.Median.correct.latency","SSP.Span.length",
             "nfixbreak","best_acc_m_exclude","first_lat_m_exclude")

accvars<-c("Anti_CRR","Mix_CRR","DMS.PC","SSP.Span.length","nfixbreak_fl","SOC.Problems.solved.in.minimum.moves","best_acc_m_exclude_fl")

latvars<-c("Anti_CRLat","Mix_CRLat","SOC.Overallmeaninitialthinkingtime","DMS.Median.correct.latency","first_lat_m_exclude")

allefvars<-c(accvars,latvars)

coglongdata$visitnum<-as.numeric(coglongdata$visit)

groups<-data.frame(acc=c("Accuracycomposite","Anti_CRR","Mix_CRR","DMS.PC","SSP.Span.length","nfixbreak_fl","SOC.Problems.solved.in.minimum.moves","best_acc_m_exclude_fl"),
                   lat=c("Latencycomposite","Anti_CRLat","Mix_CRLat","DMS.Median.correct.latency",NA,NA,"SOC.Overallmeaninitialthinkingtime","first_lat_m_exclude"),
                   group=c("Composite","ANTI","MIX","DMS","SSP","FIX","SOC","MGS"))
groupslong<-gather(groups, type, outcome, acc:lat, factor_key=TRUE)

####maturation points######
baseformula<-as.formula('outcome~s(pred)+s(visitnum)')

####scalefits################
scaledfitsacc<-mgcvscalefits(coglongdata,outcomevars = c("Accuracycomposite",accvars),predvars = "Ageatvisit",idvar="id",mformula = baseformula)
scaledfitsacc$type<-"acc"
scaledfitslat<-mgcvscalefits(coglongdata,outcomevars = c("Latencycomposite",latvars),predvars = "Ageatvisit",idvar="id",mformula = baseformula)
scaledfitslat$type<-"lat"
scalefitsall<-rbind(scaledfitsacc,scaledfitslat)
scalefitsall<-merge(scalefitsall,groupslong,by=c("outcome","type"))
scalefitsall<-scalefitsall %>% group_by(outcome) %>% mutate(adjusted.p=p.adjust(unique(modelp)))
scalefitsallsig<-scalefitsall
scalefitsallsig$sig<-as.character(dplyr::if_else(scalefitsallsig$adjusted.p < .05,1,0))

Fig1statstable<-scalefitsall %>% group_by(group,outcome) %>% summarize(unique(adjusted.p),unique(edf),unique(Fstat))
write.csv(Fig1statstable,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/ScaledFits/LUNAscaledfitsstats.csv")
write.csv(scalefitsallsig,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/ScaledFits/LUNAscaledfits.csv")
###scaled points#######
scaledpointsacc<-mgcvscalefits(coglongdata,outcomevars = c("Accuracycomposite",accvars),predvars = "Ageatvisit",idvar="id",mformula = baseformula,interval_inc = 3)
scaledpointsacc$type<-"acc"

scaledpointslat<-mgcvscalefits(coglongdata,outcomevars = c("Latencycomposite",latvars),predvars = "Ageatvisit",idvar="id",mformula = baseformula,interval_inc = 3)
scaledpointslat$type<-"lat"

scalepointssall<-rbind(scaledpointsacc,scaledpointslat)
scalepointssallsig<-scalepointssall[scalepointssall$outcome %in% unique(scalefitsallsig$outcome),]
scalepointssallsig<-merge(scalepointssallsig,groupslong,by=c("outcome","type"))

ggscaledall<-ggplot()+
  geom_line(data=scalefitsallsig[scalefitsallsig$group=="Composite",],aes(x=pred,y=fit),colour="black",size=1.5,alpha=.5)+
  geom_line(data=scalefitsallsig[scalefitsallsig$group!="Composite",],aes(x=pred,y=fit,colour=group,linetype=sig),size=.5,alpha=1)+
  scale_x_continuous(limits = c(8, 36),breaks=(c(10,15,20,25,30,35)))+
  geom_point(data=scalepointssallsig[scalepointssallsig$group!="Composite",],aes(x=pred,y=fit,shape=group),size=2)+
  facet_grid(rows=vars(type),scales="free")+theme(strip.text.y= element_text(angle = 180))+
  scale_shape_manual(values=c(0,1,2,3,4,5,6,7))+
  scale_colour_manual(values=rep("black",length(unique(scalefitsallsig$group))))+
  scale_linetype_manual(values=c("dashed","solid"))
ggscaledall<-LNCDR::lunaize(ggscaledall)+xlab("Age (years)")+ylab("Executive Function\n (z-score)")+theme(legend.position = "top",legend.title = element_blank())+
  theme(strip.text.y= element_text(angle = 360))
ggscaledall
ggsave(ggscaledall,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/ScaledFits/LUNAscaledfits.pdf",height=6.5,width=8)
ggscaledallLUNA<-ggscaledall
save(ggscaledallLUNA,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/ScaledFits/LUNAscaledfits.Rdata")

#####panel plots######Figure 3
###acc
accvarsfromgrouplong<-groupslong[groupslong$type=="acc",]
accvarsfromgrouplong<-accvarsfromgrouplong[accvarsfromgrouplong$group!="Composite",]
accvarsfromgrouplongorder<-accvarsfromgrouplong[order(accvarsfromgrouplong$group),]
accvarsplotorder<-accvarsfromgrouplongorder$outcome
accvarsplotorder<-accvarsplotorder[!is.na(accvarsplotorder)]

ggderivplotsacc<-multi_mgcvgam_growthrate_multiplot_rasteronly(df=coglongdata,outcomevars = c("Accuracycomposite",accvarsplotorder),idvar="id",predvars='Ageatvisit',mformula=baseformula,datarangegrey=FALSE,derivrangemanual=c(-.2,.2),derivcolourmanual=c("#2f42bde6","#c9270e"))
noaorder<-match(as.character(ggderivplotsacc$pairs$outcome[!is.na(ggderivplotsacc$pairs$outcome)]),accvarsfromgrouplongorder$outcome)+1 ##offset 1 for composite
noaorder<-c(1,noaorder[!is.na(noaorder)]) ##add 1 for composite

finallist<-vector("list", length = length(ggderivplotsacc$pairs$outcome))
for(ni in 1:length(ggderivplotsacc$allplots)){
  print(ni)
  niforfinallist<-noaorder[ni]
  finallist[[niforfinallist]]<-ggderivplotsacc$allplots[[ni]]
}

naorder<-seq(1:length(ggderivplotsacc$pairs$outcome))[!seq(1:length(ggderivplotsacc$pairs$outcome)) %in% noaorder]
finalorderggmatacc<-cowplot::plot_grid(plotlist=c(finallist,list(ggderivplotsacc$tilescaleplot)),ncol = 1)

save(finalorderggmatacc,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/MatRasters/LUNAaccmatrasters.plot.Rdata")
save(ggderivplotsacc,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/MatRasters/LUNAaccmatrasters.fulldata.Rdata")


#####lat#####
latarsfromgrouplong<-groupslong[groupslong$type=="lat",]
latarsfromgrouplong<-latarsfromgrouplong[latarsfromgrouplong$group!="Composite",]
latarsfromgrouplongorder<-latarsfromgrouplong[order(latarsfromgrouplong$group),]
latvarsplotorder<-latarsfromgrouplongorder$outcome
latvarsplotorder<-latvarsplotorder[!is.na(latvarsplotorder)]

ggderivplotslat<-multi_mgcvgam_growthrate_multiplot_rasteronly(df=coglongdata,outcomevars = c("Latencycomposite",latvarsplotorder),idvar="id",predvars='Ageatvisit',mformula=baseformula,totaltiles = 8,datarangegrey=FALSE,derivrangemanual=c(-.2,.2),derivcolourmanual=c("#2f42bde6","#c9270e"))
noaorder<-match(as.character(ggderivplotslat$pairswithfill$outcome[!is.na(ggderivplotslat$pairswithfill$outcome)]),latarsfromgrouplongorder$outcome)+1 ##offset 1 for composite
noaorder<-c(1,noaorder[!is.na(noaorder)]) ##add 1 for composite

finallist<-vector("list", length = length(ggderivplotslat$pairswithfill$outcome))
for(ni in 1:length(ggderivplotslat$allplots)){
  print(ni)
  niforfinallist<-noaorder[ni]
  finallist[[niforfinallist]]<-ggderivplotslat$allplots[[ni]]
}

naorder<-seq(1:length(ggderivplotslat$pairswithfill$outcome))[!seq(1:length(ggderivplotslat$pairswithfill$outcome)) %in% noaorder]
finalorderggmatlat<-cowplot::plot_grid(plotlist=c(finallist,list(ggderivplotslat$tilescaleplot)),ncol = 1)

save(finalorderggmatlat,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/MatRasters/LUNAlatmatrasters.plot.Rdata")
save(ggderivplotslat,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/MatRasters/LUNAlatmatrasters.fulldata.Rdata")

#####scaled to 20-30###For Figure 5
source("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/GeneralScripts/metabyage.R")
coglongdatascale20<-read.csv("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/data/btc_R03cleaneddata_20220306.outlierremoved.compositeacclat.2030scale.csv")
groupslong2030<-groupslong
groupslong2030$outcome<-paste0(groupslong2030$outcome,"_20_30")
accvarsscale20<-paste0(accvars,"_20_30")
scaledfitsacc2030<-mgcvscalefits(coglongdatascale20,outcomevars = c("Accuracycomposite_20_30",accvarsscale20),predvars = "Ageatvisit",idvar="id",mformula = baseformula,scale=FALSE) ###already scaled to adults
scaledfitsacc2030$type<-"acc"
latvarsscale20<-paste0(latvars,"_20_30")
scaledfitslat2030<-mgcvscalefits(coglongdatascale20,outcomevars = c("Latencycomposite_20_30",latvarsscale20),predvars = "Ageatvisit",idvar="id",mformula = baseformula,scale=FALSE) ###already scaled to adults
scaledfitslat2030$type<-"lat"
scalefitsall2030<-rbind(scaledfitsacc2030,scaledfitslat2030)
scalefitsall2030<-merge(scalefitsall2030,groupslong2030,by=c("outcome","type"))
scalefitsall2030<-scalefitsall2030 %>% group_by(outcome) %>% mutate(adjusted.p=p.adjust(unique(modelp)))
scalefitsall2030sig<-scalefitsall2030
scalefitsall2030sig$sig<-as.character(dplyr::if_else(scalefitsall2030sig$adjusted.p < .05,1,0))
scalefitsall2030sig$dataset<-"LUNA"
save(scalefitsall2030sig,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/AdultScaledFits/LUNAallfitsadultscaled.Rdata")


