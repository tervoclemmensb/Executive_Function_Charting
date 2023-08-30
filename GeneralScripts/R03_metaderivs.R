library(patchwork)
library(ggplot2)
library(scales)
library(dplyr)
library(metafor)
library(tidyr)
source("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/GeneralScripts/growthrate_mgcvgam.R")
#########all deriv data####
##read#######
###luna###
load(file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/MatRasters/LUNAaccmatrasters.fulldata.Rdata")
Lunaderivacc<-ggderivplotsacc$derivdata
Lunaderivacc$dataset<-"Luna"
load(file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/MatRasters/LUNAlatmatrasters.fulldata.Rdata")
Lunaderivlat<-ggderivplotslat$derivdata
Lunaderivlat$dataset<-"Luna"
###NCANDA
load(file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/MatRasters/NCANDAaccmatrasters.fulldata.Rdata")
NCANDAderivacc<-ggderivplotsacc$derivdata
NCANDAderivacc$dataset<-"NCANDA"
load(file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/MatRasters/NCANDAlatmatrasters.fulldata.Rdata")
NCANDAderivlat<-ggderivplotslat$derivdata
NCANDAderivlat$dataset<-"NCANDA"
###NKI######
load(file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/MatRasters/NKIaccmatrasters.fulldata.Rdata")
NKIderivacc<-ggderivplotsacc$derivdata
NKIderivacc$dataset<-"NKI"
load(file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/MatRasters/NKIlatmatrasters.fulldata.Rdata")
NKIderivlat<-ggderivplotslat$derivdata
NKIderivlat$dataset<-"NKI"
###PNC######
load(file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/MatRasters/PNCaccmatrasters.fulldata.Rdata")
PNCderivacc<-ggderivplotsacc$derivdata
PNCderivacc$dataset<-"PNC"
load(file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/MatRasters/PNClatmatrasters.fulldata.Rdata")
PNCderivlat<-ggderivplotslat$derivdata
PNCderivlat$dataset<-"PNC"

######bind rows####
allderivdatalat<-plyr::rbind.fill(Lunaderivlat,NCANDAderivlat) %>% plyr::rbind.fill(.,NKIderivlat) %>% plyr::rbind.fill(.,PNCderivlat)
allderivdatalat<-data.frame(allderivdatalat %>% group_by(dataset,var) %>% mutate(mean_dff_cumsum=cumsum(mean_dff)) %>% 
                              mutate(mean_dff_cumsum_total_diff=max(mean_dff_cumsum,na.rm=TRUE)-min(mean_dff_cumsum,na.rm=TRUE),
                                     mean_dff_cumsum_precentotal=mean_dff_cumsum/mean_dff_cumsum_total_diff,numberofsigages=length(which(mean_dff_clip!=0))))

allderivdataacc<-plyr::rbind.fill(Lunaderivacc,NCANDAderivacc) %>% plyr::rbind.fill(.,NKIderivacc) %>% plyr::rbind.fill(.,PNCderivacc)
allderivdataacc<-data.frame(allderivdataacc %>% group_by(dataset,var) %>% mutate(mean_dff_cumsum=cumsum(mean_dff)) %>% 
                              mutate(mean_dff_cumsum_total_diff=max(mean_dff_cumsum,na.rm=TRUE)-min(mean_dff_cumsum,na.rm=TRUE),
                                     mean_dff_cumsum_precentotal=mean_dff_cumsum/mean_dff_cumsum_total_diff,numberofsigages=length(which(mean_dff_clip!=0))))

allderivdataacc$varbydataset<-paste(allderivdataacc$var,allderivdataacc$dataset)
allderivdatalat$varbydataset<-paste(allderivdatalat$var,allderivdatalat$dataset)

####visualize all mat points#######
matpointsacc<-allderivdataacc %>% group_by(varbydataset) %>% summarize(matpoint=unique(matpoint))
matpointslat<-allderivdatalat %>% group_by(varbydataset) %>% summarize(matpoint=unique(matpoint))

gghistmatacc<-ggplot(matpointsacc,aes(x=matpoint))+geom_histogram(bins=30,colour="black",fill="#c9270e",alpha=.99)+scale_y_continuous(expand=c(0,0))+
  geom_vline(xintercept = median(matpointsacc$matpoint,na.rm=TRUE),linetype="dashed",colour="grey40",size=1.25)+
  scale_x_continuous(limits = c(8,36),breaks=(c(10,15,20,25,30,35)))
gghistmatacc<-LNCDR::lunaize(gghistmatacc)+xlab("Final Age of Significance (years)")+ylab("Count (# of EF Measures)")+ggtitle("Accuracy\n")+
  theme(plot.title = element_text(size=24,hjust = 0.5))

gghistmatlat<-ggplot(matpointslat,aes(x=matpoint))+geom_histogram(bins=30,colour="black",fill="#253494",alpha=.99)+scale_y_continuous(expand=c(0,0))+
  geom_vline(xintercept = median(matpointslat$matpoint,na.rm=TRUE),linetype="dashed",colour="grey40",size=1.25)+
  scale_x_continuous(limits = c(8,36),breaks=(c(10,15,20,25,30,35)))
gghistmatlat<-LNCDR::lunaize(gghistmatlat)+xlab("Final Age of Significance (years)")+ylab("Count (# of EF Measures)")+ggtitle("Latency\n")+
  theme(plot.title = element_text(size=24,hjust = 0.5))

library(patchwork)
allhists<-(gghistmatacc)+(gghistmatlat)
ggsave(allhists,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/MatRasters/finalsignificantderivage.pdf",height=5,width=9)
matpointsacc$type<-"acc"
matpointslat$type<-"lat"
allmatpoints<-rbind(matpointsacc,matpointslat)
write.csv(allmatpoints,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Data/SupportingData/Sup4.csv")

####aggregate bars for all measures####
###################################
###acc####
##first have to interpolate on to a common age grid
agegrid<-data.frame(ages=as.numeric(seq(8,35,by=.1)))
#agegrid$interp<-approx(x=t$ages, y = t$mean_dff, xout=agegrid$ages, method="linear")[[2]]

agebydatabyvargridacc<-expand.grid(ages=as.numeric(seq(8,35,by=.1)),varbydataset=unique(allderivdataacc$varbydataset))
agebydatabyvargridacc$dataset<-unlist(lapply(strsplit(as.character(agebydatabyvargridacc$varbydataset)," "),"[[",2))
agebydatabyvargridacc$var<-unlist(lapply(strsplit(as.character(agebydatabyvargridacc$varbydataset)," "),"[[",1))
agebydatabyvargridacc$varglobal<-agebydatabyvargridacc$var

PNBK<-c("SLNB2_SLNB_MCR","LNB_MCR","cnp_sfnb2_sfnb_mcr")
PCET<-c("cnp_pcet_pcet_acc2","PCET_PCET_ACC2","PCET_ACC2")
PCPT<-c("cnp_spcptnl_scpt_tp","PCPT_T_TP","SPCPTNL_SCPT_TP")

agebydatabyvargridacc$varglobal[agebydatabyvargridacc$var %in% PNBK]<-"PNBK"
agebydatabyvargridacc$varglobal[agebydatabyvargridacc$var %in% PCET]<-"PCET"
agebydatabyvargridacc$varglobal[agebydatabyvargridacc$var %in% PCPT]<-"PCPT"

interpwithmatchedagegridacc<-do.call(bind_rows,lapply(as.character(unique(agebydatabyvargridacc$varbydataset)),function(varbydataset){
  print(varbydataset)
  allderivdataacc_i<-allderivdataacc[allderivdataacc$varbydataset==varbydataset,]
  varbydatasetagerange<-range(allderivdataacc_i$ages)
  print(varbydatasetagerange)
  allderivdataacc_i$sesq<-allderivdataacc_i$se^2
  agebydatabyvargridacc_i<-agebydatabyvargridacc[agebydatabyvargridacc$varbydataset==varbydataset &
                                                   agebydatabyvargridacc$ages>=varbydatasetagerange[1] &
                                                   agebydatabyvargridacc$ages<=varbydatasetagerange[2],]
  agebydatabyvargridacc_i$mean_dff_interp<-approx(x=allderivdataacc_i$ages, y = allderivdataacc_i$mean_dff, xout=agebydatabyvargridacc_i$ages, method="linear")[[2]]
  agebydatabyvargridacc_i$sesq_interp<-approx(x=allderivdataacc_i$ages, y = allderivdataacc_i$sesq, xout=agebydatabyvargridacc_i$ages, method="linear")[[2]]
  return(agebydatabyvargridacc_i)
}))
####meta#####
metabyageacc<-lapply(unique(agegrid$ages),function(ai){
  
  #https://cjvanlissa.github.io/Doing-Meta-Analysis-in-R/fitting-a-three-level-model.html
  #https://pdfs.semanticscholar.org/e12f/da2b1c4dbef76ddd1a80ad79756ba8b37fe6.pdf
  ####The variable labeled v contains the sampling variance that corresponds with the observed effect size in the variable y and can be obtained by
  #squaring the standard error
  print(ai)
  aidataacc<-interpwithmatchedagegridacc[interpwithmatchedagegridacc$ages==ai,]
  metaout<-data.frame(age=ai,ndatsets=length(unique(aidataacc$dataset)))
  if (length(unique(aidataacc$dataset))>1){
    m_multireturn<-tryCatch(
      {
        m_multi<-metafor::rma.mv(mean_dff_interp, sesq_interp, random = list(~ 1 | varglobal, ~ 1 | dataset), data = aidataacc,method="ML")
        broom::tidy(m_multi)
      },
      error = function(e){
        NA
      }
    )
    metaout<-cbind(metaout,m_multireturn)
  }
  return(metaout)
})
metabyageaccdf<-do.call(plyr::rbind.fill,metabyageacc)
metabyageaccdf$estimate[metabyageaccdf$age >=10 & metabyageaccdf$age <=15]
metabyageaccdf$estimate[metabyageaccdf$age >=10 & metabyageaccdf$age <=12]
metabyageaccdf$estimate[metabyageaccdf$age >=10 & metabyageaccdf$age <=12]
metabyageaccdf$estimate[metabyageaccdf$age >=18 & metabyageaccdf$age <=20]
metabyageaccdf$estimateclip<-metabyageaccdf$estimate
metabyageaccdf$estimateclip[metabyageaccdf$p.value > .05]<-0 ###set non significant to 0 for plotting 
metabyageaccdf$estimateclip[is.na(metabyageaccdf$p.value)]<-0

deriv_range<-range(metabyageaccdf$estimateclip)
deriv_range<-c(-.2,.2)

tilemetaacc <- ggplot(metabyageaccdf,aes(x=age,y=1,fill=estimateclip))+ geom_raster(interpolate = TRUE) +
  scale_fill_gradient2(low = "#253494", mid = "white", high = "#c9270e", midpoint = 0, space = "Lab", breaks = sort(c(0, deriv_range)), limits = deriv_range,oob=squish)+xlim(c(8,36))
tilemetaacclegend <- cowplot::get_legend(tilemetaacc)
tilemetaacc <- lunaize_geomraster(tilemetaacc) + theme(text = element_text(size = 20))+theme(axis.title.x=element_blank())
#d73f27
#D73027
##dummy scale bar for matching ages###
ciscaleplot<-metabyageaccdf
ciscaleplot$estimate<-0
tilescaleplot <- ggplot(ciscaleplot,aes(x=age,y=1,fill=estimate))+ geom_raster(interpolate = TRUE) +
  scale_fill_gradient2(low = "#2f42bde6", mid = "white", high = "#c9270e", midpoint = 0, space = "Lab", breaks = sort(c(0, deriv_range)), limits = deriv_range)+
  scale_x_continuous(limits = c(8,36),breaks=(c(10,15,20,25,30,35)))
tilescaleplot <- lunaize_geomrasterxkeep(tilescaleplot)+theme(text = element_text(size = 20))+theme(axis.title.x=element_blank())
tilescaleplot<-tilescaleplot+geom_segment(linetype = 2, colour = "black", aes(x = 12, xend =  12, y = 0.5, yend = 1.5))
tilescaleplot<-tilescaleplot+geom_segment(linetype = 2, colour = "black", aes(x = 15, xend=15, y = 0.5, yend = 1.5))
tilescaleplot<-tilescaleplot+geom_segment(linetype = 2, colour = "black", aes(x = 18, xend=18, y = 0.5, yend = 1.5))

##finallist for spacing###
finalmetaacclist<-lapply(1:8,function(ni){NULL}) ###repeating for spacing  
finalmetaacclist[[1]]<-tilemetaacc

finalggmetamatacc<-cowplot::plot_grid(plotlist=c(finalmetaacclist,list(tilescaleplot)),ncol = 1)

save(finalggmetamatacc,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/MatRasters/metaaccmatrasters.plot.Rdata")
save(tilemetaacclegend,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/MatRasters/metaaccmatrasters.legend.Rdata")

###lat######
##first have to interpolate on to a common age grid
agebydatabyvargridlat<-expand.grid(ages=as.numeric(seq(8,35,by=.1)),varbydataset=unique(allderivdatalat$varbydataset))
agebydatabyvargridlat$dataset<-unlist(lapply(strsplit(as.character(agebydatabyvargridlat$varbydataset)," "),"[[",2))
agebydatabyvargridlat$var<-unlist(lapply(strsplit(as.character(agebydatabyvargridlat$varbydataset)," "),"[[",1))

agebydatabyvargridlat$varglobal<-agebydatabyvargridlat$var

PNBK<-c("cnp_sfnb2_sfnb_mrtc","LNB_MRTC","SLNB2_SLNB_MRTC")
PCET<-c("cnp_pcet_pcetrtcr","PCET_RTCR","PCET_PCETRTCR")
PCPT<-c("cnp_spcptnl_scpt_tprt","PCPT_T_TPRT","SPCPTNL_SCPT_TPRT")

agebydatabyvargridlat$varglobal[agebydatabyvargridlat$var %in% PNBK]<-"PNBK"
agebydatabyvargridlat$varglobal[agebydatabyvargridlat$var %in% PCET]<-"PCET"
agebydatabyvargridlat$varglobal[agebydatabyvargridlat$var %in% PCPT]<-"PCPT"


interpwithmatchedagegridlat<-do.call(bind_rows,lapply(as.character(unique(agebydatabyvargridlat$varbydataset)),function(varbydataset){
  print(varbydataset)
  allderivdatalat_i<-allderivdatalat[allderivdatalat$varbydataset==varbydataset,]
  varbydatasetagerange<-range(allderivdatalat_i$ages)
  print(varbydatasetagerange)
  allderivdatalat_i$sesq<-allderivdatalat_i$se^2
  agebydatabyvargridlat_i<-agebydatabyvargridlat[agebydatabyvargridlat$varbydataset==varbydataset &
                                                   agebydatabyvargridlat$ages>=varbydatasetagerange[1] &
                                                   agebydatabyvargridlat$ages<=varbydatasetagerange[2],]
  agebydatabyvargridlat_i$mean_dff_interp<-approx(x=allderivdatalat_i$ages, y = allderivdatalat_i$mean_dff, xout=agebydatabyvargridlat_i$ages, method="linear")[[2]]
  agebydatabyvargridlat_i$sesq_interp<-approx(x=allderivdatalat_i$ages, y = allderivdatalat_i$sesq, xout=agebydatabyvargridlat_i$ages, method="linear")[[2]]
  return(agebydatabyvargridlat_i)
}))
####meta####
metabyagelat<-lapply(unique(agegrid$ages),function(ai){
  print(ai)
  aidatalat<-interpwithmatchedagegridlat[interpwithmatchedagegridlat$ages==ai,]
  metaout<-data.frame(age=ai,ndatsets=length(unique(aidatalat$dataset)))
  if (length(unique(aidatalat$dataset))>1){
    m_multireturn<-tryCatch(
      {
        m_multi<-metafor::rma.mv(mean_dff_interp, sesq_interp, random = list(~ 1 | varglobal, ~ 1 | dataset), data = aidatalat,method="ML")
        broom::tidy(m_multi)
      },
      error = function(e){
        NA
      }
    )
    metaout<-cbind(metaout,m_multireturn)
  }
  return(metaout)
})
metabyagelatdf<-do.call(plyr::rbind.fill,metabyagelat)
metabyagelatdfnona<-metabyagelatdf[!is.na(metabyagelatdf$estimate),]
metabyagelatdfnona$estimate[metabyagelatdfnona$age >=10 & metabyagelatdfnona$age <=15]
metabyagelatdfnona$estimateclip<-metabyagelatdfnona$estimate
metabyagelatdfnona$estimateclip[metabyagelatdfnona$p.value > .05]<-0 ###set non significant to 0 for plotting 

#deriv_range<-range(metabyagelatdfnona$estimateclip)
deriv_range<-c(-0.2,.2)

tilemetalat <- ggplot(metabyagelatdfnona,aes(x=age,y=1,fill=estimateclip))+ geom_raster(interpolate = TRUE) +
  scale_fill_gradient2(low = "#2f42bde6", mid = "white", high = "#c9270e", midpoint = 0, space = "Lab", breaks = rev(deriv_range),limits=deriv_range,
                       oob=squish)+xlim(c(8,36))
tilemetalatlegend <- cowplot::get_legend(tilemetalat)
tilemetalat <- lunaize_geomraster(tilemetalat) + theme(text = element_text(size = 20))+theme(axis.title.x=element_blank())

finalmetalatlist<-lapply(1:8,function(ni){NULL}) ###repeating for spacing  
finalmetalatlist[[1]]<-tilemetalat

finalggmetamatlat<-cowplot::plot_grid(plotlist=c(finalmetalatlist,list(tilescaleplot)),ncol = 1)

save(finalggmetamatlat,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/MatRasters/metalatmatrasters.plot.Rdata")
save(tilemetalatlegend,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/MatRasters/metalatmatrasters.legend.Rdata")


###################################supplemental all derivs unthresholded#####
##set up names###
lunagroups<-data.frame(acc=c("Accuracycomposite","Anti_CRR","Mix_CRR","DMS.PC","SSP.Span.length","nfixbreak_fl","SOC.Problems.solved.in.minimum.moves","best_acc_m_exclude_fl"),
                   lat=c("Latencycomposite","Anti_CRLat","Mix_CRLat","DMS.Median.correct.latency",NA,NA,"SOC.Overallmeaninitialthinkingtime","first_lat_m_exclude"),
                   group=c("COMP","ANTI","MIX","DMS","SSP","FIX","SOC","MGS"))
lunagroupslong<-gather(lunagroups, type, var, acc:lat, factor_key=TRUE)
lunagroupslong$dataset<-"Luna"

NCANDAgroups<-data.frame(acc=c("Accuracycomposite","cnp_pcet_pcet_acc2","cnp_spcptnl_scpt_tp","cnp_sfnb2_sfnb_mcr",NA),
                   lat=c("Latencycomposite","cnp_pcet_pcetrtcr","cnp_spcptnl_scpt_tprt","cnp_sfnb2_sfnb_mrtc","stroop_total_mean"),
                   group=c("COMP","PCET","PCTP","PNBK","STRP"))
NCANDAgroupslong<-gather(NCANDAgroups, type, var, acc:lat, factor_key=TRUE)
NCANDAgroupslong$dataset<-"NCANDA"

NKIgroups<-data.frame(acc=c("Accuracycomposite","PCET_PCET_ACC2","SPCPTNL_SCPT_TP","SLNB2_SLNB_MCR","TOWacc","DFLacc",NA,NA),
                   lat=c("Latencycomposite","PCET_PCETRTCR","SPCPTNL_SCPT_TPRT","SLNB2_SLNB_MRTC",NA,NA,"CWIlat","TMTlat"),
                   group=c("COMP","PCET","PCTP","PNBK","TOW","DFL","CWI","TMT"))
NKIgroupslong<-gather(NKIgroups, type, var, acc:lat, factor_key=TRUE)
NKIgroupslong$dataset<-"NKI"

PNCgroups<-data.frame(acc=c("Accuracycomposite","PCET_ACC2","PCPT_T_TP","LNB_MCR"),
                   lat=c("Latencycomposite","PCET_RTCR","PCPT_T_TPRT","LNB_MRTC"),
                   group=c("COMP","PCET","PCTP","PNBK"))
PNCgroupslong<-gather(PNCgroups, type, var, acc:lat, factor_key=TRUE)
PNCgroupslong$dataset<-"PNC"

####all labels#####
allabels<-plyr::rbind.fill(lunagroupslong,NCANDAgroupslong) %>% plyr::rbind.fill(.,NKIgroupslong) %>% plyr::rbind.fill(.,PNCgroupslong)

###accuracy
allderivdataaccwithlabels<-merge(allderivdataacc,allabels,by=c("dataset","var"))
allderivdataaccwithlabels$groupf<-factor(allderivdataaccwithlabels$group)
allderivdataaccwithlabels$groupf<-factor(allderivdataaccwithlabels$group,levels=c("COMP",levels(allderivdataaccwithlabels$groupf)[levels(allderivdataaccwithlabels$groupf)!="COMP"]))

allderivdataaccwithlabels$mean_dff_clip_bin<-dplyr::if_else(allderivdataaccwithlabels$mean_dff_clip!=0,1,0)
allderivdataaccwithlabels_diffcliponly<-allderivdataaccwithlabels[allderivdataaccwithlabels$mean_dff_clip_bin==1,]
allderivdataaccwithlabels_diffcliponly$mean_dff_clip_sign<-sign(allderivdataaccwithlabels_diffcliponly$mean_dff_clip)
##latency 
allderivdatalatwithlabels<-merge(allderivdatalat,allabels,by=c("dataset","var"))
allderivdatalatwithlabels$groupf<-factor(allderivdatalatwithlabels$group)
allderivdatalatwithlabels$groupf<-factor(allderivdatalatwithlabels$group,levels=c("COMP",levels(allderivdatalatwithlabels$groupf)[levels(allderivdatalatwithlabels$groupf)!="COMP"]))

allderivdatalatwithlabels$mean_dff_clip_bin<-dplyr::if_else(allderivdatalatwithlabels$mean_dff_clip!=0,1,0)
allderivdatalatwithlabels_diffcliponly<-allderivdatalatwithlabels[allderivdatalatwithlabels$mean_dff_clip_bin==1,]
allderivdatalatwithlabels_diffcliponly$mean_dff_clip_sign<-sign(allderivdatalatwithlabels_diffcliponly$mean_dff_clip)

########plots#######
##acc
ggallderivCIacc<-ggplot(allderivdataaccwithlabels,aes(x=ages,y=mean_dff))+geom_line(size=.5)+
  geom_ribbon(aes(x=ages,ymin=ci_low,ymax=ci_high),fill="grey55",colour="grey55",alpha=.2)+
  geom_ribbon(data=allderivdataaccwithlabels_diffcliponly,aes(x=ages,ymin=ci_low,ymax=ci_high,fill=as.factor(mean_dff_clip_sign),colour=as.factor(mean_dff_clip_sign)),alpha=.2)+
  geom_line(data=allderivdataaccwithlabels_diffcliponly,aes(x=ages,y=mean_dff,colour=as.factor(mean_dff_clip_sign)))+
  scale_colour_manual(values=c("#c9270e","grey80"))+
  scale_x_continuous(limits = c(8, 36),breaks=(c(10,15,20,25,30,35)))+
  geom_hline(yintercept = 0,linetype="dashed")+
  #facet_grid(cols=vars(groupf),rows=vars(dataset),scales="free",space="free")+
  facet_wrap(~dataset*groupf,scales = "free")+
  theme(legend.position = "none")

ggallderivCIacc<-LNCDR::lunaize(ggallderivCIacc)+xlab("Age (years)")+ylab("First Derivative of Age Trajectory\n (z-score)")+theme(legend.position = "none",legend.title = element_blank())+
  theme(axis.text.x = element_text(size=12),axis.text.y = element_text(size=12),strip.text.x=element_text(size=10))+ggtitle("Accuracy\n")+
  theme(plot.title = element_text(hjust = 0.5,size=24))

ggsave(ggallderivCIacc,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/MatRasters/derivfullCI.acc.pdf",height=9,width=10)

###lat

ggallderivCIlat<-ggplot(allderivdatalatwithlabels,aes(x=ages,y=mean_dff))+geom_line(size=.5)+
  geom_ribbon(aes(x=ages,ymin=ci_low,ymax=ci_high),fill="grey55",colour="grey55",alpha=.2)+
  geom_ribbon(data=allderivdatalatwithlabels_diffcliponly,aes(x=ages,ymin=ci_low,ymax=ci_high,fill=as.factor(mean_dff_clip_sign),colour=as.factor(mean_dff_clip_sign)),alpha=.2)+
  geom_line(data=allderivdatalatwithlabels_diffcliponly,aes(x=ages,y=mean_dff,colour=as.factor(mean_dff_clip_sign)))+
  scale_colour_manual(values=c("#2f42bde6","#c9270e"))+
  scale_fill_manual(values=c("#2f42bde6","#c9270e"))+
  scale_x_continuous(limits = c(8, 36),breaks=(c(10,15,20,25,30,35)))+
  geom_hline(yintercept = 0,linetype="dashed")+
  #facet_grid(cols=vars(groupf),rows=vars(dataset),scales="free",space="free")+
  facet_wrap(~dataset*groupf,scales = "free")+
  theme(legend.position = "none")

ggallderivCIlat<-LNCDR::lunaize(ggallderivCIlat)+xlab("Age (years)")+ylab("First Derivative of Age Trajectory\n (z-score)")+theme(legend.position = "none",legend.title = element_blank())+
  theme(axis.text.x = element_text(size=12),axis.text.y = element_text(size=12),strip.text.x=element_text(size=10))+ggtitle("Latency\n")+
  theme(plot.title = element_text(hjust = 0.5,size=24))

ggsave(ggallderivCIlat,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/MatRasters/derivfullCI.lat.pdf",height=9,width=10)

accandlatderivCI<-ggallderivCIacc/ggallderivCIlat

ggsave(accandlatderivCI,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/MatRasters/derivfullCI.pdf",height=18,width=10)

########
allderivCI<-plyr::rbind.fill(allderivdatalatwithlabels[,c("ages","type","dataset","groupf","mean_dff","ci_low","ci_high","mean_dff_clip","mean_dff_clip_bin")],
                             allderivdataaccwithlabels[,c("ages","type","dataset","groupf","mean_dff","ci_low","ci_high","mean_dff_clip","mean_dff_clip_bin")])

write.csv(allderivCI,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Data/SupportingData/Sup3.csv")

####allderivfullpanel
allderiv_individs<-plyr::rbind.fill(allderivdatalatwithlabels[,c("ages","type","dataset","groupf","mean_dff_clip")],
                             allderivdataaccwithlabels[,c("ages","type","dataset","groupf","mean_dff_clip")])

metabyagelatdfnona$type<-"lat"
metabyagelatdfnona$dataset<-"meta"
metabyagelatdfnona$groupf<-"meta"

metabyageaccdf$type<-"acc"
metabyageaccdf$dataset<-"meta"
metabyageaccdf$groupf<-"meta"


allderivmetas<-plyr::rbind.fill(metabyagelatdfnona[,c("age","type","dataset","groupf","estimateclip")],
                                metabyageaccdf[,c("age","type","dataset","groupf","estimateclip")])

names(allderivmetas)[names(allderivmetas)=="age"]<-"ages"
names(allderivmetas)[names(allderivmetas)=="estimateclip"]<-"mean_dff_clip"

allderiv_fullpanel<-plyr::rbind.fill(allderivmetas,allderiv_individs)
write.csv(allderiv_fullpanel,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Data/SupportingData/Figure2.csv")

# #########
# 
# 
# gderivpercentfit<-ggplot(allfitsacc,aes(x=pred,y=fitscale,colour=varbydataset))+geom_line(size=.5,alpha=.35)+
#   theme(legend.position = "none")+
#   geom_line(data=accmetafitsigsmooth,aes(x=pred,y=fit),colour="black",size=2)+
#   scale_x_continuous(limits = c(8, 36),breaks=(c(10,15,20,25,30,35)))+
#   scale_colour_manual(values=c("black","black",rep("black",length(unique(allfitsacc$varbydataset)))))+
#   scale_linetype_manual(values=c("solid"))
# ggderivpercentfit<-LNCDR::lunaize(ggderivpercentfit)+xlab("Age (years)")+ylab("% of Total")+theme(legend.position = "none",legend.title = element_blank())+
#   theme(strip.text.y= element_text(angle = 360))
# 
# ggsave(ggderivpercentfit,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/MatRasters/percentchangesigacc.fit.plot.pdf",height=4,width=6)
# 
# 
# allfitslat<-allfits[allfits$outcome %in% unique(allinterpolatedfitssig$outcome) & allfits$type=="lat",]
# 
# latmetafitsigsmooth<-mgcvscalefits(latmetafitsig,outcomevars = c("percentmax"),predvars = "age",mformula = as.formula("outcome~s(pred)"),scale=FALSE,interval_inc = .1) ###fit smooth
# 
# ggderivpercentfit<-ggplot(allfitslat,aes(x=pred,y=fitscale,colour=varbydataset))+geom_line(size=.5,alpha=.35)+
#   theme(legend.position = "none")+
#   geom_line(data=latmetafitsigsmooth,aes(x=pred,y=fit),colour="black",size=2)+
#   scale_x_continuous(limits = c(8, 36),breaks=(c(10,15,20,25,30,35)))+
#   scale_colour_manual(values=c("black","black",rep("black",length(unique(allfitslat$varbydataset)))))+
#   scale_linetype_manual(values=c("solid"))
# ggderivpercentfit<-LNCDR::lunaize(ggderivpercentfit)+xlab("Age (years)")+ylab("% of Total")+theme(legend.position = "none",legend.title = element_blank())+
#   theme(strip.text.y= element_text(angle = 360))
# 
# ggsave(ggderivpercentfit,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/MatRasters/percentchangesiglat.fit.plot.pdf",height=4,width=6)




