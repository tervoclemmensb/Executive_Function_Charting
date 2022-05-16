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
library(tidyr)
#######
source("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/GeneralScripts/growthrate_mgcvgam.R")
############
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
################
#########
baseformula<-as.formula('outcome~s(pred,k=5)+s(cnp_age)')
scaledfitsacc<-mgcvscalefits(coglongdata,outcomevars = c("Accuracycomposite",accvars),predvars = "visitnum",idvar="id",mformula = baseformula)
scaledfitsacc$type<-"acc"
scaledfitslat<-mgcvscalefits(coglongdata,outcomevars = c("Latencycomposite",latvars),predvars = "visitnum",idvar="id",mformula = baseformula)
scaledfitslat$type<-"lat"
scalefitsall<-rbind(scaledfitsacc,scaledfitslat)
scalefitsall<-merge(scalefitsall,groupslong,by=c("outcome","type"))
scalefitsall<-scalefitsall %>% group_by(outcome) %>% mutate(adjusted.p=p.adjust(unique(modelp)))
scalefitsallsig<-scalefitsall
scalefitsallsig$sig<-as.character(dplyr::if_else(scalefitsallsig$adjusted.p < .05,1,0))

###scaled points#######
scaledpointsacc<-mgcvscalefits(coglongdata,outcomevars = c("Accuracycomposite",accvars),predvars = "visitnum",idvar="id",mformula = baseformula,interval_inc = 3)
scaledpointsacc$type<-"acc"

scaledpointslat<-mgcvscalefits(coglongdata,outcomevars = c("Latencycomposite",latvars),predvars = "visitnum",idvar="id",mformula = baseformula,interval_inc = 3)
scaledpointslat$type<-"lat"

scalepointssall<-rbind(scaledpointsacc,scaledpointslat)
scalepointssallsig<-scalepointssall[scalepointssall$outcome %in% unique(scalefitsallsig$outcome),]
scalepointssallsig<-merge(scalepointssallsig,groupslong,by=c("outcome","type"))

ggscaledall<-ggplot()+
  geom_line(data=scalefitsallsig[scalefitsallsig$group=="Composite",],aes(x=pred,y=fit),colour="black",size=1.5,alpha=.5)+
  geom_line(data=scalefitsallsig[scalefitsallsig$group!="Composite",],aes(x=pred,y=fit,colour=group,linetype=sig),size=.5,alpha=1)+
  scale_x_continuous(limits = c(1, 10),breaks=(c(1,5,10)))+
  geom_point(data=scalepointssallsig[scalepointssallsig$group!="Composite",],aes(x=pred,y=fit,shape=group),size=2)+
  facet_grid(rows=vars(type),scales="free")+theme(strip.text.y= element_text(angle = 180))+
  scale_shape_manual(values=c(0,1,2,3,4,5,6,7))+
  scale_colour_manual(values=rep("black",length(unique(scalefitsallsig$group))))+
  scale_linetype_manual(values=c("dashed","solid"))
ggscaledall<-LNCDR::lunaize(ggscaledall)+xlab("Visit Number")+ylab("Executive Function\n (z-score)")+theme(legend.position = "top",legend.title = element_blank())+
  theme(strip.text.y= element_text(angle = 360))
ggscaledall
ggsave(ggscaledall,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/PracticeVis/NCANDAscaledvisitfits.pdf",height=6.5,width=8)
ggscaledallNCANDA<-ggscaledall
save(ggscaledallNCANDA,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/PracticeVis/NCANDAscaledvisitfits.Rdata")

