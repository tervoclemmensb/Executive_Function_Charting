library(patchwork)
library(ggplot2)
library(scales)
library(dplyr)
library(metafor)
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


###################################







