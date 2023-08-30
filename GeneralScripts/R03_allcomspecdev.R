library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpattern)
library(patchwork)
####
##LUNA###
load(file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/Agecomspec/LUNAagecomspecalldata.Rdata")
Lunaprobbyage<-savelistout$propbyage$byvardf
Lunaprobbyage<-merge(Lunaprobbyage,savelistout$groupslong,by=c("outcome"))
Lunaprobbyage$dataset<-"LUNA"

Lunaprobbyagesens<-savelistout$propbyagesens$byvardf
Lunaprobbyagesens<-merge(Lunaprobbyagesens,savelistout$groupslong,by=c("outcome"))
Lunaprobbyagesens$dataset<-"LUNA"

ggLunasens<-savelistout$gpefvarspatternsens+theme(legend.position = "none")

Lunaoutsavebyvar<-savelistout$propbyagebyvarsens
##NCANDA##
load(file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/Agecomspec/NCANDAagecomspecalldata.Rdata")
NCANDAprobbyage<-savelistout$propbyage$byvardf
NCANDAprobbyage<-merge(NCANDAprobbyage,savelistout$groupslong,by=c("outcome"))
NCANDAprobbyage$dataset<-"NCANDA"

ggNCANDA<-savelistout$gpefvarspattern+theme(legend.position = "none")

NCANDAoutsavebyvar<-savelistout$propbyagebyvar
##NKI###
load(file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/Agecomspec/NKIagecomspecalldata.Rdata")
NKIprobbyage<-savelistout$propbyage$byvardf
NKIprobbyage<-merge(NKIprobbyage,savelistout$groupslong,by=c("outcome"))
NKIprobbyage$dataset<-"NKI"

NKIprobbyagesens<-savelistout$propbyagesens$byvardf
NKIprobbyagesens<-merge(NKIprobbyagesens,savelistout$groupslong,by=c("outcome"))
NKIprobbyagesens$dataset<-"NKI"

ggNKIsens<-savelistout$gpefvarspatternsens+theme(legend.position = "none")

NKIoutsavebyvar<-savelistout$propbyagebyvar
##PNC##
load(file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/Agecomspec/PNCagecomspecalldata.Rdata")
PNCprobbyage<-savelistout$propbyage$byvardf
PNCprobbyage<-merge(PNCprobbyage,savelistout$groupslong,by=c("outcome"))
PNCprobbyage$dataset<-"PNC"
names(PNCprobbyage)[names(PNCprobbyage)=="type"]<-"vartype"
ggPNC<-savelistout$gpefvarspattern+theme(legend.position = "none")

PNCoutsavebyvar<-savelistout$propbyagebyvar
#####alldata#######
allprobbyagesens<-plyr::rbind.fill(Lunaprobbyagesens,NCANDAprobbyage) %>% plyr::rbind.fill(.,NKIprobbyagesens) %>% plyr::rbind.fill(.,PNCprobbyage) ###sens version in Luna and NKI ensures all tests are out of meausre and out of domain

###varglobal###
allprobbyagesens$varglobal<-allprobbyagesens$outcome
PNBKacc<-c("SLNB2_SLNB_MCR","LNB_MCR","cnp_sfnb2_sfnb_mcr")
PCETacc<-c("cnp_pcet_pcet_acc2","PCET_PCET_ACC2","PCET_ACC2")
PCPTacc<-c("cnp_spcptnl_scpt_tp","PCPT_T_TP","SPCPTNL_SCPT_TP")

allprobbyagesens$varglobal[allprobbyagesens$outcome %in% PNBKacc]<-"PNBKacc"
allprobbyagesens$varglobal[allprobbyagesens$outcome %in% PCETacc]<-"PCETacc"
allprobbyagesens$varglobal[allprobbyagesens$outcome %in% PCPTacc]<-"PCPTacc"
##lat###
PNBKlat<-c("cnp_sfnb2_sfnb_mrtc","LNB_MRTC","SLNB2_SLNB_MRTC")
PCETlat<-c("cnp_pcet_pcetrtcr","PCET_RTCR","PCET_PCETRTCR")
PCPTlat<-c("cnp_spcptnl_scpt_tprt","PCPT_T_TPRT","SPCPTNL_SCPT_TPRT")

allprobbyagesens$varglobal[allprobbyagesens$outcome %in% PNBKlat]<-"PNBKlat"
allprobbyagesens$varglobal[allprobbyagesens$outcome %in% PCETlat]<-"PCETlat"
allprobbyagesens$varglobal[allprobbyagesens$outcome %in% PCPTlat]<-"PCPTlat"

allprobbyagesens$sepropofage<-allprobbyagesens$sefullminusmodcomp/allprobbyagesens$devmodelbaseage
allprobbyagesens$sepropofage_sq<-allprobbyagesens$sepropofage^2
allprobbyagesens$specdev_propofage[allprobbyagesens$specdev_propofage<0]<-0 ###setting to negative to zero for aggregation
allprobbyagesens$specdev_propofage[allprobbyagesens$specdev_propofage>1]<-1 ###setting to over one to one for aggregation


specprobmetaaccandlat<-data.frame(allprobbyagesens %>% do(m_multi=(metafor::rma.mv(specdev_propofage, V=sepropofage_sq, random = list(~ 1 | varglobal, ~ 1 | dataset), data = .,method="ML") %>% 
                                                                                                   broom::tidy(.))) %>% unnest(m_multi))
specprobmetaaccandlat$nonspecdev<-1-specprobmetaaccandlat$estimate

specprobmeta<-data.frame(allprobbyagesens %>% group_by(vartype) %>% do(m_multi=(metafor::rma.mv(specdev_propofage, V=sepropofage_sq, random = list(~ 1 | varglobal, ~ 1 | dataset), data = .,method="ML") %>% 
                    broom::tidy(.))) %>% unnest(m_multi))    

specprobmeta$nonspecdev<-1-specprobmeta$estimate
specprobmeta$type<-NULL
specprobmeta_long<-pivot_longer(specprobmeta,cols=c(estimate,nonspecdev), names_to="type", values_to="prop")
specprobmeta_long$typef<-factor(ifelse(specprobmeta_long$type=="estimate","Measure Specific","Common EF"),levels=c("Measure Specific","Common EF"))
specprobmeta_long$prop100<-specprobmeta_long$prop*100

scalexdiscretelimits=c(levels(specprobmeta_long$vartype),sprintf("Null%s",seq(1:10)))###adding nulls to ensure equal bar width across datasets


gpefvarspattern_meta<-ggplot(specprobmeta_long,aes(x=vartype,y=prop100,fill=typef))+geom_bar_pattern(aes(pattern=typef),position="stack", stat="identity",colour="black",alpha=.95,
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
patternlegend<-cowplot::get_legend(gpefvarspattern_meta)
grid::grid.draw(patternlegend)
gpefvarspattern_meta<-LNCDR::lunaize(gpefvarspattern_meta)+theme(legend.title = element_blank(),legend.position = "top")+ylab("% of Age-Related EF")+xlab("")+
  theme(axis.text.x = element_text(angle = 90,size=15))+theme(legend.position = "none")

####

allspecsensgg<-(ggLunasens+ggNCANDA)/(ggNKIsens+ggPNC)+(gpefvarspattern_meta+plot_spacer())
ggsave(allspecsensgg,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/Agecomspec/allcomspec.pdf",height=12,width=10)

########save outdata
meta_save_r<-specprobmeta_long[,c("vartype","typef","prop100")]
meta_save_r$dataset<-"meta"
meta_save_r$groupbytypef<-"meta"

####individual

Lunaoutsavebyvar
NCANDAoutsavebyvar
NKIoutsavebyvar
PNCoutsavebyvar
allbyvarsaveout<-plyr::rbind.fill(Lunaoutsavebyvar,NCANDAoutsavebyvar) %>% 
  plyr::rbind.fill(.,NKIoutsavebyvar) %>% plyr::rbind.fill(PNCoutsavebyvar)

allbyvarsaveout_r<-allbyvarsaveout[,c("groupbytypef","vartype","typef","prop100","dataset")]
allbyvarsaveout_rwithmeta<-plyr::rbind.fill(allbyvarsaveout_r,meta_save_r)

write.csv(allbyvarsaveout_rwithmeta,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Data/SupportingData/Figure4.csv")






















########iter pred#######

load(file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/Agecomspec/LUNAiterpropbyage.Rdata")
propbyageitersens_LUNA<-propbyageitersens_save
propbyageitersens_LUNA$nonspecdev_propofage[propbyageitersens_LUNA$nonspecdev_propofage<0]<-0###setting to negative to zero for aggregation
propbyageitersens_LUNA$nonspecdev_propofage[propbyageitersens_LUNA$nonspecdev_propofage>1]<-1 ###setting over one to one for aggregation

propbyageitersens_LUNA$nonspecdev_propofage100<-propbyageitersens_LUNA$nonspecdev_propofage*100

ggplot(propbyageitersens_LUNA,aes(x=as.factor(thismanyvars),y=nonspecdev_propofage100))+geom_boxplot()
propbyageitersens_LUNA$dataset<-"Luna"
load(file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/Agecomspec/NKIiterpropbyage.Rdata")
propbyageitersens_NKI<-propbyageitersens_save
propbyageitersens_NKI$nonspecdev_propofage[propbyageitersens_NKI$nonspecdev_propofage<0]<-0###setting to negative to zero for aggregation
propbyageitersens_NKI$nonspecdev_propofage[propbyageitersens_NKI$nonspecdev_propofage>1]<-1 ###setting over one to one for aggregation

propbyageitersens_NKI$nonspecdev_propofage100<-propbyageitersens_NKI$nonspecdev_propofage*100

ggplot(propbyageitersens_NKI,aes(x=as.factor(thismanyvars),y=nonspecdev_propofage100))+geom_boxplot()
propbyageitersens_NKI$dataset<-"NKI"
########
allpropdata_LUNANKI<-plyr::rbind.fill(propbyageitersens_LUNA,propbyageitersens_NKI)
medians<-allpropdata_LUNANKI %>% group_by(dataset,thismanyvars) %>% dplyr::summarize(nonspecdev_propofage100med=median(nonspecdev_propofage100),se=sd(nonspecdev_propofage100)/sqrt(n()),
                                                                                     upper=nonspecdev_propofage100med+2*se,lower=nonspecdev_propofage100med-2*se,uniquevals=length(unique(nonspecdev_propofage100)),uniqueoutcomes=length(unique(outcome)))
propbyvarnumber<-ggplot(medians,aes(x=thismanyvars,y=nonspecdev_propofage100med,ymin=lower,ymax=upper,shape=dataset))+geom_point(colour="grey28")+geom_smooth(method="lm",se=FALSE,colour="grey28")+
  scale_x_continuous(breaks=c(2:max(medians$thismanyvars)))
propbyvarnumber<-LNCDR::lunaize(propbyvarnumber)+xlab("# Vars Included in Composite")+ylab("% of Age-Related EF\nby Common EF")+theme(legend.position = "top",legend.title = element_blank())


###allpropdata_LUNANKI limit to 7 vars for Luna as this is the max where all outcomes are included#####
allpropdata_LUNANKI[allpropdata_LUNANKI$dataset=="Luna" & allpropdata_LUNANKI$thismanyvars>7,]<-NA
allpropdata_LUNANKI<-allpropdata_LUNANKI[!is.na(allpropdata_LUNANKI$outcome),]

propbyvarnumber<-ggplot(allpropdata_LUNANKI,aes(x=as.factor(thismanyvars),y=nonspecdev_propofage100))+geom_boxplot(colour="grey55")+geom_smooth(method = "gam",formula=y~s(x,k=6),se=FALSE, aes(group=1),colour="grey28")+facet_grid(cols=vars(dataset),scale="free_x")
propbyvarnumber<-LNCDR::lunaize(propbyvarnumber)+xlab("# Vars Included in Composite")+ylab("% of Age-Related EF\nby Common EF")+theme(legend.position = "top",legend.title = element_blank())
ggsave(propbyvarnumber,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/Agecomspec/allcomspec.box.pdf",height=12,width=10)

supportingdata_r_s<-allpropdata_LUNANKI[,c("thismanyvars","nonspecdev_propofage100","dataset")]

write.csv(supportingdata_r_s,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Data/SupportingData/Sup10.csv")

allspecsensggwithscatter<-(ggLunasens+ggNCANDA)/(ggNKIsens+ggPNC)/(gpefvarspattern_meta+propbyvarnumber)
ggsave(propbyvarnumber,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/Agecomspec/allcomspec.withscatter.pdf",height=12,width=10)
