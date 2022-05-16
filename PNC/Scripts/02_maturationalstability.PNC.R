library(dplyr)
library(tidyr)
library(lubridate)
library(ggplot2)
library(lsmeans)
library(mgcv)
library(itsadug)
library(lme4)
library(stats)
library(psych)
library(LNCDR)
library(FactoMineR)
library(corrplot)
library(mgcv)
library(cowplot)
library(RColorBrewer)
library(wesanderson)
library(grid)
######functions#####
source("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/GeneralScripts/growthrate_mgcvgam.R")
############dataread
coglongdata<-read.csv("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/PNC/Data/btc_PNCscoredmeasures_20220306.outlierremoved.compositeacclat.csv")

allcogtestsacc<-c("PADT_A","PFMT_IFAC_TOT","PEIT_CR","PWMT_KIWRD_TOT","PVRT_CR","PEDT_A","PMAT_CR","VOLT_SVT","LNB_MCR","PCET_ACC2","PCPT_T_TP","PLOT_TC","WRAT_CR_RAW")
allcogtestslatency<-c("PADT_T","PFMT_IFAC_RTC","PEIT_CRT","PWMT_KIWRD_RTC","PVRT_RTCR","PEDT_T","MP_MP2RTCR","PMAT_RTCR","VOLT_SVTCRT","LNB_MRTC","PCET_RTCR","PCPT_T_TPRT","PLOT_TCRT")

allcogvars<-c(allcogtestsacc,allcogtestslatency)
allfactorvars<-allcogvars

testgroups<-c("PADT","PFMT","PEIT","PWMT","PVRT","PEDT","PMAT","VOLT","LNB","PCET","PCPT","PLOT","WRAT","MP")
efficiencytestgroups<-c("PADT","PFMT","PEIT","PWMT","PVRT","PEDT","PMAT","VOLT","LNB","PCPT","PLOT")
eftestgroups<-c("PCET","PCPT","LNB") ###EXECUTIVE-CONTROL#HTTPS://WWW.NCBI.NLM.NIH.GOV/PMC/ARTICLES/PMC3295891/PDF/NIHMS348849.PDF

accvars<-unlist(lapply(eftestgroups,function(vg){grep(vg,allcogtestsacc,value=TRUE)}))
latvars<-unlist(lapply(eftestgroups,function(vg){grep(vg,allcogtestslatency,value=TRUE)}))

allefvars<-c(accvars,latvars) 
groups<-data.frame(acc=c("Accuracycomposite","PCET_ACC2","PCPT_T_TP","LNB_MCR"),
                   lat=c("Latencycomposite","PCET_RTCR","PCPT_T_TPRT","LNB_MRTC"),
                   group=c("Composite","PCET","PCTP","PNBK"))
groupslong<-gather(groups, type, outcome, acc:lat, factor_key=TRUE)
####maturation points######
baseformula<-as.formula("outcome~s(pred)")

####scalefits################
scaledfitsacc<-mgcvscalefits(coglongdata,outcomevars = c("Accuracycomposite",accvars),predvars = "age",mformula = baseformula)
scaledfitsacc$type<-"acc"
scaledfitslat<-mgcvscalefits(coglongdata,outcomevars = c("Latencycomposite",latvars),predvars = "age",mformula = baseformula)
scaledfitslat$type<-"lat"
scalefitsall<-rbind(scaledfitsacc,scaledfitslat)
scalefitsall<-merge(scalefitsall,groupslong,by=c("outcome","type"))
scalefitsall<-scalefitsall %>% group_by(outcome) %>% mutate(adjusted.p=p.adjust(unique(modelp)))
scalefitsallsig<-scalefitsall
scalefitsallsig$sig<-as.character(dplyr::if_else(scalefitsallsig$adjusted.p < .05,1,0))

Fig1statstable<-scalefitsall %>% group_by(group,outcome) %>% summarize(unique(adjusted.p),unique(edf),unique(Fstat))
write.csv(Fig1statstable,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/ScaledFits/PNCcaledfitsstats.csv")
write.csv(scalefitsallsig,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/ScaledFits/PNCscaledfits.csv")

###scaled points#######
scaledpointsacc<-mgcvscalefits(coglongdata,outcomevars = c("Accuracycomposite",accvars),predvars = "age",mformula = baseformula,interval_inc = 3)
scaledpointsacc$type<-"acc"

scaledpointslat<-mgcvscalefits(coglongdata,outcomevars = c("Latencycomposite",latvars),predvars = "age",mformula = baseformula,interval_inc = 3)
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
  scale_shape_manual(values=c(0,1,2))+
  scale_colour_manual(values=rep("black",length(unique(scalefitsallsig$group))))+
  scale_linetype_manual(values=c("solid")) ###all are sig so only one level
ggscaledall<-LNCDR::lunaize(ggscaledall)+xlab("Age (years)")+ylab("Executive Function\n (z-score)")+theme(legend.position = "top",legend.title = element_blank())+
  theme(strip.text.y= element_text(angle = 360))
ggscaledall
ggsave(ggscaledall,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/ScaledFits/PNCscaledfits.pdf",height=6.5,width=8)
ggscaledallPNC<-ggscaledall
save(ggscaledallPNC,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/ScaledFits/PNCscaledfits.Rdata")

#####panel plots######Figure 3
accvarsfromgrouplong<-groupslong[groupslong$type=="acc",]
accvarsfromgrouplong<-accvarsfromgrouplong[accvarsfromgrouplong$group!="Composite",]
accvarsfromgrouplongorder<-accvarsfromgrouplong[order(accvarsfromgrouplong$group),]
accvarsplotorder<-accvarsfromgrouplongorder$outcome
accvarsplotorder<-accvarsplotorder[!is.na(accvarsplotorder)]

ggderivplotsacc<-multi_mgcvgam_growthrate_multiplot_rasteronly(df=coglongdata,outcomevars = c("Accuracycomposite",accvarsplotorder),predvars='age',mformula=baseformula,totaltiles = 8,datarangegrey=TRUE,xrangetickvals=c(8,35),derivrangemanual=c(-.2,.2),derivcolourmanual=c("#2f42bde6","#c9270e"))
noaorder<-match(as.character(ggderivplotsacc$pairs$outcome[!is.na(ggderivplotsacc$pairs$outcome)]),accvarsfromgrouplongorder$outcome)+1 ##offset 1 for composite
noaorder<-c(1,noaorder[!is.na(noaorder)]) ##add 1 for composite

finallist<-vector("list", length = length(ggderivplotsacc$pairs$outcome))
for(ni in 1:length(ggderivplotsacc$allplots)){
  print(ni)
  niforfinallist<-noaorder[ni]
  finallist[[niforfinallist]]<-ggderivplotsacc$allplots[[ni]]
}

naorder<-seq(1:length(ggderivplotsacc$pairs$outcome))[!seq(1:length(ggderivplotsacc$pairs$outcome)) %in% noaorder]
###PNC has four extras and only there for spacing, setting to NA
for (nulli in 8:(max(length(c("Accuracycomposite",accvarsplotorder)))+1)){
  print(nulli)
  finallist[[nulli]]<-NA
}

finalorderggmatacc<-cowplot::plot_grid(plotlist=c(finallist,list(ggderivplotsacc$tilescaleplot)),ncol = 1)

save(finalorderggmatacc,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/MatRasters/PNCaccmatrasters.plot.Rdata")
save(ggderivplotsacc,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/MatRasters/PNCaccmatrasters.fulldata.Rdata")


#####lat#####
latarsfromgrouplong<-groupslong[groupslong$type=="lat",]
latarsfromgrouplong<-latarsfromgrouplong[latarsfromgrouplong$group!="Composite",]
latarsfromgrouplongorder<-latarsfromgrouplong[order(latarsfromgrouplong$group),]
latvarsplotorder<-latarsfromgrouplongorder$outcome
latvarsplotorder<-latvarsplotorder[!is.na(latvarsplotorder)]

ggderivplotslat<-multi_mgcvgam_growthrate_multiplot_rasteronly(df=coglongdata,outcomevars = c("Latencycomposite",latvarsplotorder),predvars='age',mformula=baseformula,totaltiles = 8,datarangegrey=TRUE,xrangetickvals=c(8,35),derivrangemanual=c(-.2,.2),derivcolourmanual=c("#2f42bde6","#c9270e"))
noaorder<-match(as.character(ggderivplotslat$pairswithfill$outcome[!is.na(ggderivplotslat$pairswithfill$outcome)]),latarsfromgrouplongorder$outcome)+1 ##offset 1 for composite
noaorder<-c(1,noaorder[!is.na(noaorder)]) ##add 1 for composite

finallist<-vector("list", length = length(ggderivplotslat$pairswithfill$outcome))
for(ni in 1:length(ggderivplotslat$allplots)){
  print(ni)
  niforfinallist<-noaorder[ni]
  finallist[[niforfinallist]]<-ggderivplotslat$allplots[[ni]]
}

naorder<-seq(1:length(ggderivplotslat$pairswithfill$outcome))[!seq(1:length(ggderivplotslat$pairswithfill$outcome)) %in% noaorder]
###PNC has four extras and only there for spacing, setting to NA
for (nulli in 8:(max(length(c("Accuracycomposite",accvarsplotorder)))+1)){
  print(nulli)
  finallist[[nulli]]<-NA
}

finalorderggmatlat<-cowplot::plot_grid(plotlist=c(finallist,list(ggderivplotslat$tilescaleplot)),ncol = 1)

save(finalorderggmatlat,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/MatRasters/PNClatmatrasters.plot.Rdata")
save(ggderivplotslat,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/MatRasters/PNClatmatrasters.fulldata.Rdata")

#####scaled to 20-30###For Figure 5
source("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/GeneralScripts/metabyage.R")
coglongdatascale20<-read.csv("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/PNC/data/btc_R03cleaneddata_20220306.outlierremoved.compositeacclat.2030scale.csv")
groupslong2030<-groupslong
groupslong2030$outcome<-paste0(groupslong2030$outcome,"_20_30")
accvarsscale20<-paste0(accvars,"_20_30")
scaledfitsacc2030<-mgcvscalefits(coglongdatascale20,outcomevars = c("Accuracycomposite_20_30",accvarsscale20),predvars = "age",mformula = baseformula,scale=FALSE)### already scaled to adults
scaledfitsacc2030$type<-"acc"
latvarsscale20<-paste0(latvars,"_20_30")
scaledfitslat2030<-mgcvscalefits(coglongdatascale20,outcomevars = c("Latencycomposite_20_30",latvarsscale20),predvars = "age",mformula = baseformula,scale=FALSE) ### already scaled to adults
scaledfitslat2030$type<-"lat"
scalefitsall2030<-rbind(scaledfitsacc2030,scaledfitslat2030)
scalefitsall2030<-merge(scalefitsall2030,groupslong2030,by=c("outcome","type"))
scalefitsall2030<-scalefitsall2030 %>% group_by(outcome) %>% mutate(adjusted.p=p.adjust(unique(modelp)))
scalefitsall2030sig<-scalefitsall2030
scalefitsall2030sig$sig<-as.character(dplyr::if_else(scalefitsall2030sig$adjusted.p < .05,1,0))
scalefitsall2030sig$dataset<-"PNC"
save(scalefitsall2030sig,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/AdultScaledFits/PNCallfitsadultscaled.Rdata")
