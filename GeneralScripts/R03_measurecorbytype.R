library(dplyr)
library(tidyr)
library(patchwork)
library(ggplot2)
########load all data#######
load(file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/Factorfigs/allfactordataLuna.Rdata")

Lunabaseline<-data.frame(alleigssaves$corsbaseline)
Lunabaseline$type<-"baseline"

Lunalongitudinal<-data.frame(alleigssaves$corswithin)
Lunalongitudinal$type<-"longitudinal"
Lunalongitudinal[,c("var1","var2","varorder")]<-lapply(Lunalongitudinal[,c("var1","var2","varorder")],function(x){gsub(".wg","",x)})

allefvars<-c("Anti_CRLat","Anti_CRR","Mix_CRLat","Mix_CRR","SOC.Problems.solved.in.minimum.moves",
             "SOC.Overallmeaninitialthinkingtime","DMS.PC","DMS.Median.correct.latency","SSP.Span.length",
             "nfixbreak","best_acc_m_exclude","first_lat_m_exclude")

accvars<-c("Anti_CRR","Mix_CRR","DMS.PC","SSP.Span.length","nfixbreak_fl","SOC.Problems.solved.in.minimum.moves","best_acc_m_exclude_fl")

latvars<-c("Anti_CRLat","Mix_CRLat","SOC.Overallmeaninitialthinkingtime","DMS.Median.correct.latency","first_lat_m_exclude")

allefvars<-c(accvars,latvars)

groups<-data.frame(acc=c("Accuracycomposite","Anti_CRR","Mix_CRR","DMS.PC","SSP.Span.length","nfixbreak_fl","SOC.Problems.solved.in.minimum.moves","best_acc_m_exclude_fl"),
                   lat=c("Latencycomposite","Anti_CRLat","Mix_CRLat","DMS.Median.correct.latency",NA,NA,"SOC.Overallmeaninitialthinkingtime","first_lat_m_exclude"),
                   group=c("Composite","ANTI","MIX","DMS","SSP","FIX","SOC","MGS"))
groupslong<-gather(groups, type, outcome, acc:lat, factor_key=TRUE)

var1group<-groupslong
names(var1group)<-paste(names(var1group),"var1",sep="")
names(var1group)[names(var1group)=="outcomevar1"]<-"var1"

var2group<-groupslong
names(var2group)<-paste(names(var2group),"var2",sep="")
names(var2group)[names(var2group)=="outcomevar2"]<-"var2"

Lunabaseline<-merge(Lunabaseline,var1group,by=c("var1"))
Lunabaseline<-merge(Lunabaseline,var2group,by=c("var2"))

Lunalongitudinal<-merge(Lunalongitudinal,var1group,by=c("var1"))
Lunalongitudinal<-merge(Lunalongitudinal,var2group,by=c("var2"))

allluna<-plyr::rbind.fill(Lunabaseline,Lunalongitudinal)
allluna$dataset<-"Luna"

allluna$typepair<-paste(allluna$typevar1,allluna$typevar2,sep="_")
allluna$nbetween<-alleigssaves$nbetween
allluna$nwithin<-alleigssaves$nwithin

####
load(file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/Factorfigs/allfactordataNKI.Rdata")

NKI<-data.frame(alleigssaves$corslong)
NKI$dataset<-"NKI"
NKI$type<-"baseline"

accvarsNKI<-c("PCET_PCET_ACC2","SLNB2_SLNB_MCR","SPCPTNL_SCPT_TP","TOWacc","DFLacc")
latvarsNKI<-c("PCET_PCETRTCR","SLNB2_SLNB_MRTC","SPCPTNL_SCPT_TPRT","CWIlat","TMTlat")

allefvarsNKI<-c(accvarsNKI,latvarsNKI)

groupsNKI<-data.frame(acc=c("Accuracycomposite","PCET_PCET_ACC2","SPCPTNL_SCPT_TP","SLNB2_SLNB_MCR","TOWacc","DFLacc",NA,NA),
                   lat=c("Latencycomposite","PCET_PCETRTCR","SPCPTNL_SCPT_TPRT","SLNB2_SLNB_MRTC",NA,NA,"CWIlat","TMTlat"),
                   group=c("Composite","PCET","PCTP","PNBK","TOW","DFL","CWI","TMT"))
groupslongNKI<-gather(groupsNKI, type, outcome, acc:lat, factor_key=TRUE)

var1group<-groupslongNKI
names(var1group)<-paste(names(var1group),"var1",sep="")
names(var1group)[names(var1group)=="outcomevar1"]<-"var1"

var2group<-groupslongNKI
names(var2group)<-paste(names(var2group),"var2",sep="")
names(var2group)[names(var2group)=="outcomevar2"]<-"var2"

NKI<-merge(NKI,var1group,by=c("var1"))
NKI<-merge(NKI,var2group,by=c("var2"))

NKI$typepair<-paste(NKI$typevar1,NKI$typevar2,sep="_")
NKI$nbetween<-alleigssaves$nbetween

#########
load(file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/Factorfigs/allfactordataPNC.Rdata")

PNC<-data.frame(alleigssaves$corslong)
PNC$dataset<-"PNC"
PNC$type<-"baseline"

allcogtestsacc<-c("PADT_A","PFMT_IFAC_TOT","PEIT_CR","PWMT_KIWRD_TOT","PVRT_CR","PEDT_A","PMAT_CR","VOLT_SVT","LNB_MCR","PCET_ACC2","PCPT_T_TP","PLOT_TC","WRAT_CR_RAW")
allcogtestslatency<-c("PADT_T","PFMT_IFAC_RTC","PEIT_CRT","PWMT_KIWRD_RTC","PVRT_RTCR","PEDT_T","MP_MP2RTCR","PMAT_RTCR","VOLT_SVTCRT","LNB_MRTC","PCET_RTCR","PCPT_T_TPRT","PLOT_TCRT")

allcogvars<-c(allcogtestsacc,allcogtestslatency)
allfactorvars<-allcogvars

testgroups<-c("PADT","PFMT","PEIT","PWMT","PVRT","PEDT","PMAT","VOLT","LNB","PCET","PCPT","PLOT","WRAT","MP")
efficiencytestgroups<-c("PADT","PFMT","PEIT","PWMT","PVRT","PEDT","PMAT","VOLT","LNB","PCPT","PLOT")
eftestgroups<-c("PCET","PCPT","LNB") ###EXECUTIVE-CONTROL#HTTPS://WWW.NCBI.NLM.NIH.GOV/PMC/ARTICLES/PMC3295891/PDF/NIHMS348849.PDF

accvarsPNC<-unlist(lapply(eftestgroups,function(vg){grep(vg,allcogtestsacc,value=TRUE)}))
latvarsPNC<-unlist(lapply(eftestgroups,function(vg){grep(vg,allcogtestslatency,value=TRUE)}))

allefvarsPNC<-c(accvarsPNC,latvarsPNC) 
groupsPNC<-data.frame(acc=c("Accuracycomposite","PCET_ACC2","PCPT_T_TP","LNB_MCR"),
                   lat=c("Latencycomposite","PCET_RTCR","PCPT_T_TPRT","LNB_MRTC"),
                   group=c("Composite","PCET","PCTP","LNB"))
groupslongPNC<-gather(groupsPNC, type, outcome, acc:lat, factor_key=TRUE)

var1group<-groupslongPNC
names(var1group)<-paste(names(var1group),"var1",sep="")
names(var1group)[names(var1group)=="outcomevar1"]<-"var1"

var2group<-groupslongPNC
names(var2group)<-paste(names(var2group),"var2",sep="")
names(var2group)[names(var2group)=="outcomevar2"]<-"var2"

PNC<-merge(PNC,var1group,by=c("var1"))
PNC<-merge(PNC,var2group,by=c("var2"))

PNC$typepair<-paste(PNC$typevar1,PNC$typevar2,sep="_")
PNC$nbetween<-alleigssaves$nbetween

ggplot(PNC,aes(x=typepair,y=cor))+geom_boxplot(outlier.shape = NA)+geom_jitter()
################
load(file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/Factorfigs/allfactordataNCANDA.Rdata")

NCANDAbaseline<-data.frame(alleigssaves$corsbaseline)
NCANDAbaseline$type<-"baseline"

NCANDAlongitudinal<-data.frame(alleigssaves$corswithin)
NCANDAlongitudinal$type<-"longitudinal"
NCANDAlongitudinal[,c("var1","var2","varorder")]<-lapply(NCANDAlongitudinal[,c("var1","var2","varorder")],function(x){gsub(".wg","",x)})

allNCANDA<-plyr::rbind.fill(NCANDAbaseline,NCANDAlongitudinal) 
allNCANDA$dataset<-"NCANDA"

allaccvarsNCANDA<-c("cnp_cpf_ifac_tot","cnp_cpw_iwrd_tot","cnp_spcptnl_scpt_tp","cnp_sfnb2_sfnb_mcr","cnp_pmat24a_pmat24_a_cr","cnp_cpfd_dfac_tot","cnp_cpwd_dwrd_tot",
              "cnp_shortvolt_svt","cnp_er40d_er40_cr","cnp_pcet_pcet_acc2","cnp_medf36_medf36_a","cnp_pvoc_pvoccr","cnp_pvrt_pvrt_pc","cnp_svdelay_svt_ld","np_wais4_rawscore_computed","latentdd","Accuracycomposite")   

alllatvarsNCANDA<-c("cnp_cpf_ifac_rtc","cnp_cpw_iwrd_rtc","cnp_spcptnl_scpt_tprt","cnp_sfnb2_sfnb_mrtc","cnp_pmat24a_pmat24_a_rtcr","cnp_cpfd_dfac_rtc",
              "cnp_cpwd_dwrd_rtc","cnp_shortvolt_svtcrt","cnp_er40d_er40_crt","cnp_pcet_pcetrtcr","cnp_medf36_medf36_t","cnp_pvoc_pvocrtcr","cnp_pvrt_pvrtrtcr","cnp_svdelay_svtldrtc","stroop_total_mean","latentgroove","Latencycomposite")

accvarsNCANDA<-c("cnp_pcet_pcet_acc2","cnp_sfnb2_sfnb_mcr","cnp_spcptnl_scpt_tp")
latvarsNCANDA<-c("cnp_pcet_pcetrtcr","cnp_sfnb2_sfnb_mrtc","cnp_spcptnl_scpt_tprt","stroop_total_mean")

allefvarsNCANDA<-c(accvarsNCANDA,latvarsNCANDA)

groupsNCANDA<-data.frame(acc=c("Accuracycomposite","cnp_pcet_pcet_acc2","cnp_spcptnl_scpt_tp","cnp_sfnb2_sfnb_mcr",NA),
                   lat=c("Latencycomposite","cnp_pcet_pcetrtcr","cnp_spcptnl_scpt_tprt","cnp_sfnb2_sfnb_mrtc","stroop_total_mean"),
                   group=c("Composite","PCET","PCTP","PNBK","STRP"))
groupslongNCANDA<-gather(groupsNCANDA, type, outcome, acc:lat, factor_key=TRUE)

var1group<-groupslongNCANDA
names(var1group)<-paste(names(var1group),"var1",sep="")
names(var1group)[names(var1group)=="outcomevar1"]<-"var1"

var2group<-groupslongNCANDA
names(var2group)<-paste(names(var2group),"var2",sep="")
names(var2group)[names(var2group)=="outcomevar2"]<-"var2"

allNCANDA<-merge(allNCANDA,var1group,by=c("var1"))
allNCANDA<-merge(allNCANDA,var2group,by=c("var2"))

allNCANDA$typepair<-paste(allNCANDA$typevar1,allNCANDA$typevar2,sep="_")
allNCANDA$nbetween<-alleigssaves$nbetween
allNCANDA$nwithin<-alleigssaves$nwithin

###all datasets##########

alldata<-plyr::rbind.fill(allluna,NKI) %>% plyr::rbind.fill(.,PNC) %>% plyr::rbind.fill(.,allNCANDA)

alldata$label<-alldata$dataset
alldata$label[alldata$dataset %in% c("Luna","NCANDA")]<-paste(alldata$label[alldata$dataset %in% c("Luna","NCANDA")],alldata$type[alldata$dataset %in% c("Luna","NCANDA")],sep="\n")

ggplotallcors<-ggplot()+geom_hline(yintercept = 0,linetype="dotted",colour="grey66")+geom_boxplot(data=alldata,aes(x=label,y=cor),outlier.shape = NA,colour="grey66")+geom_jitter(data=alldata,aes(x=label,y=cor),width=.25,height=0)+facet_grid(rows=vars(typepair),scale="free")
ggplotallcors<-LNCDR::lunaize(ggplotallcors)+ylab("EF Measure Correlation")+xlab("")

ggplotallcors<-ggplot()+geom_hline(yintercept = 0,linetype="dotted",colour="grey66",size=.60)+geom_jitter(data=alldata,aes(x=label,y=cor,shape=label),width=.25,height=0)+facet_grid(rows=vars(typepair),scale="free")
ggplotallcors<-LNCDR::lunaize(ggplotallcors)+ylab("Executive Function\nMeasure Correlation")+xlab("")+theme(legend.position = "none")+theme(panel.border =element_rect(colour="black",fill=NA))

#######combine measures by type (three level meta)######
PNBKacc<-c("SLNB2_SLNB_MCR","LNB_MCR","cnp_sfnb2_sfnb_mcr")
PCETacc<-c("cnp_pcet_pcet_acc2","PCET_PCET_ACC2","PCET_ACC2")
PCPTacc<-c("cnp_spcptnl_scpt_tp","PCPT_T_TP","SPCPTNL_SCPT_TP")

alldata$varglobal2<-alldata$var2
alldata$varglobal1<-alldata$var1

alldata$varglobal2[alldata$var2 %in% PNBKacc]<-"PNBKacc"
alldata$varglobal2[alldata$var2 %in% PCETacc]<-"PCETacc"
alldata$varglobal2[alldata$var2 %in% PCPTacc]<-"PCPTacc"

alldata$varglobal1[alldata$var1 %in% PNBKacc]<-"PNBKacc"
alldata$varglobal1[alldata$var1 %in% PCETacc]<-"PCETacc"
alldata$varglobal1[alldata$var1 %in% PCPTacc]<-"PCPTacc"

PNBKlat<-c("cnp_sfnb2_sfnb_mrtc","LNB_MRTC","SLNB2_SLNB_MRTC")
PCETlat<-c("cnp_pcet_pcetrtcr","PCET_RTCR","PCET_PCETRTCR")
PCPTlat<-c("cnp_spcptnl_scpt_tprt","PCPT_T_TPRT","SPCPTNL_SCPT_TPRT")

alldata$varglobal2[alldata$var2 %in% PNBKlat]<-"PNBKlat"
alldata$varglobal2[alldata$var2 %in% PCETlat]<-"PCETlat"
alldata$varglobal2[alldata$var2 %in% PCPTlat]<-"PCPTlat"

alldata$varglobal1[alldata$var1 %in% PNBKlat]<-"PNBKlat"
alldata$varglobal1[alldata$var1 %in% PCETlat]<-"PCETlat"
alldata$varglobal1[alldata$var1 %in% PCPTlat]<-"PCPTlat"

alldata$varglobalpair<-paste(alldata$varglobal1,alldata$varglobal2,sep="_")

source("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/GeneralScripts/metabyage.R")
alldata$corSE<-NA
alldata$corSE[which(alldata$type=="baseline")]<-unlist(lapply(which(alldata$type=="baseline"),function(ri){
  print(ri)
  se<-corSEfromRandN(r=alldata$cor[ri],n=alldata$nbetween[ri])
  print(c(alldata$cor[ri],n=alldata$nbetween[ri],se))
  return(se)
}))

alldata$corSE[which(alldata$type=="longitudinal")]<-unlist(lapply(which(alldata$type=="longitudinal"),function(ri){
  print(ri)
  se<-corSEfromRandN(r=alldata$cor[ri],n=alldata$nwithin[ri])
  print(c(alldata$cor[ri],n=alldata$nwithin[ri],se))
  return(se)
}))

meansbytype<-lapply(unique(alldata$typepair),function(typi){
  typidat<-alldata[alldata$typepair==typi,]
  print(typi)
  typidat$corSEsqV<-typidat$corSE^2
  m_multi<-metafor::rma.mv(cor, V=corSEsqV, random = list(~ 1 | varglobalpair, ~ 1 | dataset), data = typidat,method="ML")
  m_multiout<-broom::tidy(m_multi)
  m_multiout$typepair<-typi
  return(m_multiout)
})
meansbytypedf<-do.call(rbind,meansbytype)
meansbytypedf$label<-"All\nMeasures"
meansbytypedf$labelf<-factor(meansbytypedf$label,levels=c(sort(unique(alldata$label)),"All\nMeasures"))
alldata$labelf<-factor(alldata$label,levels=c(sort(unique(alldata$label)),"All\nMeasures"))

alldatawithmean<-plyr::rbind.fill(alldata,meansbytypedf)
alldatawithmean$typepairf<-gsub("_","\nw/\n ",alldatawithmean$typepair)

ggplotallcors<-ggplot()+geom_hline(yintercept = 0,linetype="dotted",colour="black")+geom_pointrange(data=alldatawithmean,aes(x=labelf,y=estimate,ymin=estimate-(2*std.error),ymax=estimate+(2*std.error)),fatten=2)+geom_jitter(data=alldatawithmean,aes(x=labelf,y=cor,shape=labelf),width=.25,height=0)+
  facet_grid(rows=vars(typepairf),scale="free")
ggplotallcors<-LNCDR::lunaize(ggplotallcors)+ylab("Executive Function\nMeasure Correlation (r)")+xlab("")+theme(legend.position = "none")+theme(panel.border =element_rect(colour="black",fill=NA))+
  theme(axis.text.x =element_text(size=14,angle=45,vjust=1,hjust=1))+theme(strip.text.y= element_text(angle = 360))

load(file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/Factorfigs/allfactormainfigs.Rdata")


bars<-allsaveplot$allloaings[[1]]+allsaveplot$allloaings[[2]]+allsaveplot$allloaings[[3]]+allsaveplot$allloaings[[4]]+allsaveplot$allloaings[[5]]+allsaveplot$allloadingplots[[6]]+plot_layout(ncol = 6)+theme(strip.placement = NULL)


fullfigurepanel<-ggplotallcors/allsaveplot$ggallvar/allsaveplot$ggallscreeforpanel+plot_layout(heights=c(2,1,.4))
  

alldatawithmean_r<-alldatawithmean[c("labelf","typepairf","cor")]
write.csv(alldatawithmean_r,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Data/SupportingData/Figure3A.csv")


#   plot_spacer()/bars
# 
# # 
# # fullfigurepanel<-ggplotallcors/plot_spacer()/allsaveplot$ggallvar/allsaveplot$ggallscreeforpanel/plot_spacer()/(allsaveplot$allloaings[[1]]+allsaveplot$allloaings[[2]]+allsaveplot$allloaings[[3]]+allsaveplot$allloaings[[4]]+allsaveplot$allloaings[[5]]+allsaveplot$allloadingplots[[6]]+plot_layout(ncol = 6)+theme(strip.placement = NULL))+plot_layout(heights=c(2,.1,1,.33,.1,1))
# # ,widths=c(1,1,1,1,1,1))
# fullfigurepanel<-ggplotallcors/plot_spacer()/allsaveplot$ggallvar/allsaveplot$ggallscreeforpanel/plot_spacer()/(allsaveplot$allloaings[[1]]+allsaveplot$allloaings[[2]]+allsaveplot$allloaings[[3]]+allsaveplot$allloaings[[4]]+allsaveplot$allloaings[[5]]+allsaveplot$allloadingplots[[6]])+plot_layout(ncol=6)+plot_layout(ncol=6,heights=c(2,.1,1,.33,.1,2))
# 
# # ,widths=c(1,1,1,1,3))
# # +plot_layout(ncol = 6))
# # +plot_layout(widths=c(1,1,1,1,3))

ggsave(fullfigurepanel,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/Factorfigs/allfactormainfigs.figpanelminusbars.pdf",height=11,width=7)


