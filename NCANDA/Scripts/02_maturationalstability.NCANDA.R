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
#########
source("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/GeneralScripts/growthrate_mgcvgam.R")
###########READ DATA#########
#############################
coglongdata<-read.csv("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/NCANDA/Data/btc_NCANDAscoredmeasures_20220306.outlierremoved.compositeacclat.csv")
coglongdata$id<-as.factor(coglongdata$subject)
coglongdata$visitnum<-as.numeric(as.factor(coglongdata$visit))
####vars############
####################
allaccvars<-c("cnp_cpf_ifac_tot","cnp_cpw_iwrd_tot","cnp_spcptnl_scpt_tp","cnp_sfnb2_sfnb_mcr","cnp_pmat24a_pmat24_a_cr","cnp_cpfd_dfac_tot","cnp_cpwd_dwrd_tot",
           "cnp_shortvolt_svt","cnp_er40d_er40_cr","cnp_pcet_pcet_acc2","cnp_medf36_medf36_a","cnp_pvoc_pvoccr","cnp_pvrt_pvrt_pc","cnp_svdelay_svt_ld","np_wais4_rawscore_computed","latentdd","Accuracycomposite")   

alllatvars<-c("cnp_cpf_ifac_rtc","cnp_cpw_iwrd_rtc","cnp_spcptnl_scpt_tprt","cnp_sfnb2_sfnb_mrtc","cnp_pmat24a_pmat24_a_rtcr","cnp_cpfd_dfac_rtc",
           "cnp_cpwd_dwrd_rtc","cnp_shortvolt_svtcrt","cnp_er40d_er40_crt","cnp_pcet_pcetrtcr","cnp_medf36_medf36_t","cnp_pvoc_pvocrtcr","cnp_pvrt_pvrtrtcr","cnp_svdelay_svtldrtc","stroop_total_mean","latentgroove","Latencycomposite")

accvars<-c("cnp_pcet_pcet_acc2","cnp_sfnb2_sfnb_mcr","cnp_spcptnl_scpt_tp")
latvars<-c("cnp_pcet_pcetrtcr","cnp_sfnb2_sfnb_mrtc","cnp_spcptnl_scpt_tprt","stroop_total_mean")

allefvars<-c(accvars,latvars)

groups<-data.frame(acc=c("Accuracycomposite","cnp_pcet_pcet_acc2","cnp_spcptnl_scpt_tp","cnp_sfnb2_sfnb_mcr",NA),
                   lat=c("Latencycomposite","cnp_pcet_pcetrtcr","cnp_spcptnl_scpt_tprt","cnp_sfnb2_sfnb_mrtc","stroop_total_mean"),
                   group=c("Composite","PCET","PCTP","PNBK","STRP"))
groupslong<-gather(groups, type, outcome, acc:lat, factor_key=TRUE)
baseformula<-as.formula('outcome~s(pred)+s(visitnum,k=5)')
#####scalefits################
scaledfitsacc<-mgcvscalefits(coglongdata,outcomevars = c("Accuracycomposite",accvars),predvars = "cnp_age",idvar="id",mformula = baseformula)
scaledfitsacc$type<-"acc"
scaledfitslat<-mgcvscalefits(coglongdata,outcomevars = c("Latencycomposite",latvars),predvars = "cnp_age",idvar="id",mformula = baseformula)
scaledfitslat$type<-"lat"
scalefitsall<-rbind(scaledfitsacc,scaledfitslat)
scalefitsall<-merge(scalefitsall,groupslong,by=c("outcome","type"))
scalefitsall<-scalefitsall %>% group_by(outcome) %>% mutate(adjusted.p=p.adjust(unique(modelp)))
scalefitsallsig<-scalefitsall
scalefitsallsig$sig<-as.character(dplyr::if_else(scalefitsallsig$adjusted.p < .05,1,0))

Fig1statstable<-scalefitsall %>% group_by(group,outcome) %>% summarize(unique(adjusted.p),unique(edf),unique(Fstat))
write.csv(Fig1statstable,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/ScaledFits/NCANDAscaledfitsstats.csv")
write.csv(scalefitsallsig,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/ScaledFits/NCANDAscaledfits.csv")
###scaled points#######
scaledpointsacc<-mgcvscalefits(coglongdata,outcomevars = c("Accuracycomposite",accvars),predvars = "cnp_age",idvar="id",mformula = baseformula,interval_inc = 3)
scaledpointsacc$type<-"acc"

scaledpointslat<-mgcvscalefits(coglongdata,outcomevars = c("Latencycomposite",latvars),predvars = "cnp_age",idvar="id",mformula = baseformula,interval_inc = 3)
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
  scale_shape_manual(values=c(0,1,2,3))+
  scale_colour_manual(values=rep("black",length(unique(scalefitsallsig$group))))+
  scale_linetype_manual(values=c("dashed","solid"))
ggscaledall<-LNCDR::lunaize(ggscaledall)+xlab("Age (years)")+ylab("Executive Function\n (z-score)")+theme(legend.position = "top",legend.title = element_blank())+
  theme(strip.text.y= element_text(angle = 360))
ggscaledall
ggsave(ggscaledall,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/ScaledFits/NCANDAscaledfits.pdf",height=6.5,width=8)
ggscaledallNCANDA<-ggscaledall
save(ggscaledallNCANDA,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/ScaledFits/NCANDAscaledfits.Rdata")

#####panel plots######Figure 3
accvarsfromgrouplong<-groupslong[groupslong$type=="acc",]
accvarsfromgrouplong<-accvarsfromgrouplong[accvarsfromgrouplong$group!="Composite",]
accvarsfromgrouplongorder<-accvarsfromgrouplong[order(accvarsfromgrouplong$group),]
accvarsplotorder<-accvarsfromgrouplongorder$outcome
accvarsplotorder<-accvarsplotorder[!is.na(accvarsplotorder)]

ggderivplotsacc<-multi_mgcvgam_growthrate_multiplot_rasteronly(df=coglongdata,outcomevars = c("Accuracycomposite",accvarsplotorder),idvar="id",predvars='cnp_age',mformula=baseformula,totaltiles = 8,datarangegrey=TRUE,xrangetickvals=c(8,35),derivrangemanual=c(-.2,.2),derivcolourmanual=c("#2f42bde6","#c9270e"))
noaorder<-match(as.character(ggderivplotsacc$pairs$outcome[!is.na(ggderivplotsacc$pairs$outcome)]),accvarsfromgrouplongorder$outcome)+1 ##offset 1 for composite
noaorder<-c(1,noaorder[!is.na(noaorder)]) ##add 1 for composite

finallist<-vector("list", length = length(ggderivplotsacc$pairs$outcome))
for(ni in 1:length(ggderivplotsacc$allplots)){
  print(ni)
  niforfinallist<-noaorder[ni]
  finallist[[niforfinallist]]<-ggderivplotsacc$allplots[[ni]]
}

naorder<-seq(1:length(ggderivplotsacc$pairs$outcome))[!seq(1:length(ggderivplotsacc$pairs$outcome)) %in% noaorder]
###NCANDA has four extras and only there for spacing, setting to NA
for (nulli in 8:(max(length(c("Accuracycomposite",accvarsplotorder)))+1)){
  print(nulli)
  finallist[[nulli]]<-NA
}

finalorderggmatacc<-cowplot::plot_grid(plotlist=c(finallist,list(ggderivplotsacc$tilescaleplot)),ncol = 1)

save(finalorderggmatacc,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/MatRasters/NCANDAaccmatrasters.plot.Rdata")
save(ggderivplotsacc,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/MatRasters/NCANDAaccmatrasters.fulldata.Rdata")


#####lat#####
latarsfromgrouplong<-groupslong[groupslong$type=="lat",]
latarsfromgrouplong<-latarsfromgrouplong[latarsfromgrouplong$group!="Composite",]
latarsfromgrouplongorder<-latarsfromgrouplong[order(latarsfromgrouplong$group),]
latvarsplotorder<-latarsfromgrouplongorder$outcome
latvarsplotorder<-latvarsplotorder[!is.na(latvarsplotorder)]

ggderivplotslat<-multi_mgcvgam_growthrate_multiplot_rasteronly(df=coglongdata,outcomevars = c("Latencycomposite",latvarsplotorder),idvar="id",predvars='cnp_age',mformula=baseformula,totaltiles = 8,datarangegrey=TRUE,xrangetickvals=c(8,35),derivrangemanual=c(-.2,.2),derivcolourmanual=c("#2f42bde6","#c9270e"))
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

save(finalorderggmatlat,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/MatRasters/NCANDAlatmatrasters.plot.Rdata")
save(ggderivplotslat,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/MatRasters/NCANDAlatmatrasters.fulldata.Rdata")

#####scaled to 20-30###For Figure 5
source("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/GeneralScripts/metabyage.R")
coglongdatascale20<-read.csv("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/NCANDA/data/btc_R03cleaneddata_20220306.outlierremoved.compositeacclat.2030scale.csv")
coglongdatascale20$visitnum<-as.numeric(as.factor(coglongdatascale20$visit))
groupslong2030<-groupslong
groupslong2030$outcome<-paste0(groupslong2030$outcome,"_20_30")
accvarsscale20<-paste0(accvars,"_20_30")
scaledfitsacc2030<-mgcvscalefits(coglongdatascale20,outcomevars = c("Accuracycomposite_20_30",accvarsscale20),predvars = "cnp_age",idvar="id",mformula = baseformula,scale=FALSE)### already scaled to adults
scaledfitsacc2030$type<-"acc"
latvarsscale20<-paste0(latvars,"_20_30")
scaledfitslat2030<-mgcvscalefits(coglongdatascale20,outcomevars = c("Latencycomposite_20_30",latvarsscale20),predvars = "cnp_age",idvar="id",mformula = baseformula,scale=FALSE) ### already scaled to adults
scaledfitslat2030$type<-"lat"
scalefitsall2030<-rbind(scaledfitsacc2030,scaledfitslat2030)
scalefitsall2030<-merge(scalefitsall2030,groupslong2030,by=c("outcome","type"))
scalefitsall2030<-scalefitsall2030 %>% group_by(outcome) %>% mutate(adjusted.p=p.adjust(unique(modelp)))
scalefitsall2030sig<-scalefitsall2030
scalefitsall2030sig$sig<-as.character(dplyr::if_else(scalefitsall2030sig$adjusted.p < .05,1,0))
scalefitsall2030sig$dataset<-"NCANDA"
save(scalefitsall2030sig,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/AdultScaledFits/NCANDAallfitsadultscaled.Rdata")


