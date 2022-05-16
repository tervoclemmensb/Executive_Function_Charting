library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpattern)
#####
source("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/GeneralScripts/metabyage.R")
########
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
groupslong$groupbytype<-paste(groupslong$group,groupslong$type,sep="_")

groupslongcompares<-groupslong[groupslong$group!="Composite" & !is.na(groupslong$outcome),]
###set up commvar
groupslongcompares_sp<-groupslongcompares %>% dplyr::group_by(type) %>% dplyr::mutate(othervars=paste(outcome,collapse=" ")) ###remove vari from commonvar [leave one measure out]
groupslongcompares_sp$othervarsnooutcome<-unlist(lapply(1:nrow(groupslongcompares_sp),function(ri){
  print(ri)
  newothervars<-gsub(groupslongcompares_sp$outcome[ri],"",groupslongcompares_sp$othervars[ri])
  return(newothervars)
}))
groupslongcompares_sp<-data.frame(groupslongcompares_sp)
####set up commvar from other type (acc/lat) 
##select measures from other type
groupslongcompares_formergeother<-groupslongcompares_sp %>% group_by(type) %>% dplyr::summarize(othervars_signflip=paste(outcome,collapse=" "))
groupslongcompares_formergeother$type<-dplyr::if_else(as.character(groupslongcompares_formergeother$type)=="lat","acc","lat") ##switch acc/lat label (to fill in remaining othervars)
groupslongcompares_sp<-merge(groupslongcompares_sp,groupslongcompares_formergeother,by=c("type"))

groupslongcompares_sp$othervarsnooutcome_signflip<-unlist(lapply(1:nrow(groupslongcompares_sp),function(ri){ ###remove opposite measure from same task (e.g., no lat measure from same task for acc outcome)
  print(ri)
  thisgroup<-groupslongcompares_sp$group[ri]
  othervarsnooutcome_signflip<-gsub(paste(groupslongcompares_sp$outcome[groupslongcompares_sp$group==thisgroup],collapse="|"),"",groupslongcompares_sp$othervars_signflip[ri])
  return(othervarsnooutcome_signflip)
}))

######################
coglongdata_baseline<-coglongdata[coglongdata$visitnum==1,c("id",allefvars,"cnp_age")]
coglongdata_baseline[,allefvars]<-lapply(coglongdata_baseline[,allefvars],function(x){as.numeric(as.character(x))})
df<-coglongdata_baseline
agevar="cnp_age"
comparedf=groupslongcompares_sp
propbyage<-propspecbyage_gam_signflip(df,agevar,comparedf)
propbyagebyvar<-propbyage$byvardf
propbyagebyvar$dataset<-"LUNA"


propbyagebyvar<-pivot_longer(propbyagebyvar,cols=c(nonspecdev_propofage,specdev_propofage), names_to="type", values_to="prop")
propbyagebyvar$typef<-factor(ifelse(propbyagebyvar$type=="specdev_propofage","Measure Specific","Common EF"),levels=c("Measure Specific","Common EF"))
propbyagebyvar$prop[propbyagebyvar$prop<0]<-0
propbyagebyvar$prop[propbyagebyvar$prop>1]<-1

names(groupslong)[names(groupslong)=="type"]<-"vartype"
propbyagebyvar<-merge(propbyagebyvar,groupslong,by=c("outcome"))

propbyagebyvar$varname<-paste(propbyagebyvar$group,propbyagebyvar$vartype,sep="_")
propbyagebyvar$typebygroup<-paste(propbyagebyvar$vartype,propbyagebyvar$group,sep="_") ###type first for ordering in plot
propbyagebyvar$groupbytypef<-factor(propbyagebyvar$varname,levels=unique(propbyagebyvar$varname[order(propbyagebyvar$typebygroup)]))
propbyagebyvar$prop100<-propbyagebyvar$prop*100

# gpefvars<-ggplot(propbyagebyvar,aes(x=groupbytypef,y=prop100,fill=typef))+geom_bar(position="stack", stat="identity",colour="black",alpha=.95)+
#   ylab("")+xlab("")+theme(legend.title = element_blank())+
#   scale_y_continuous(limits=c(0,110),expand = c(0,0),breaks=c(0,25,50,75,100))+theme(strip.text.x = element_blank())+scale_fill_manual(values=c("grey92","grey28"))+
#   ylab("")
# gpefvars<-LNCDR::lunaize(gpefvars)+theme(legend.title = element_blank(),legend.position = "top")+ylab("% of Age-Related EF")+xlab("")+
#   theme(axis.text.x = element_text(angle = 90,size=15))####switched to hashed bar for easier reading

scalexdiscretelimits=c(levels(factor(propbyagebyvar$groupbytypef)),sprintf("Null%s",seq(12-length(levels(factor(propbyagebyvar$groupbytypef))))))###no added needed because 12 here

gpefvarspattern<-ggplot(propbyagebyvar,aes(x=groupbytypef,y=prop100,fill=typef))+geom_bar_pattern(aes(pattern=typef),position="stack", stat="identity",colour="black",alpha=.95,
                                                                                                  pattern_colour  = 'grey55',
                                                                                                  pattern_angle = 45,
                                                                                                  pattern_density = 0.005,
                                                                                                  pattern_spacing = 0.018,
                                                                                                  pattern_key_scale_factor = 0.6)+
  ylab("")+xlab("")+theme(legend.title = element_blank())+
  scale_y_continuous(limits=c(0,110),expand = c(0,0),breaks=c(0,25,50,75,100))+theme(strip.text.x = element_blank())+scale_fill_manual(values=c("grey92","grey28"))+
  ylab("")+
  scale_pattern_manual(values = c("stripe","none"))+
  scale_x_discrete(limits = scalexdiscretelimits)
gpefvarspattern<-LNCDR::lunaize(gpefvarspattern)+theme(legend.title = element_blank(),legend.position = "top")+ylab("% of Age-Related EF")+xlab("")+
  theme(axis.text.x = element_text(angle = 90,size=15))

####savelist######
savelistout<-list(propbyage=propbyage,gpefvarspattern=gpefvarspattern,groupslong=groupslong)
save(savelistout,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/Agecomspec/NCANDAagecomspecalldata.Rdata")

