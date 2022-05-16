library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpattern)
#####
source("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/GeneralScripts/metabyage.R")
########
coglongdata<-read.csv("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/data/btc_R03cleaneddata_20220306.outlierremoved.compositeacclat.csv")
coglongdata$id<-as.factor(coglongdata$id)
coglongdata$visitnum<-as.numeric(coglongdata$visit)
####vars############
####################
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
coglongdata_baseline<-coglongdata[coglongdata$visitnum==1,c("id",allefvars,"Ageatvisit")]
coglongdata_baseline[,allefvars]<-lapply(coglongdata_baseline[,allefvars],function(x){as.numeric(as.character(x))})
df<-coglongdata_baseline
agevar="Ageatvisit"
comparedf=groupslongcompares_sp

#########sensivity removing same domain
domainpairs<-data.frame(domain=c("inhibition","workingmem"),measures=c(paste(c("Anti_CRR","nfixbreak_fl","Mix_CRR","Anti_CRLat","Mix_CRLat"),collapse="+"),paste(c("DMS.PC","SSP.Span.length","best_acc_m_exclude_fl","DMS.Median.correct.latency","first_lat_m_exclude"),collapse="+")))
comparedfnodomainpairs<-comparedf
for (dpi in 1:nrow(domainpairs)){
  print(domainpairs$domain[dpi])
  domainvars<-unlist(strsplit(domainpairs$measures[dpi],"[+]"))
  for (dmvi in domainvars){
    print(dmvi)
    dmvnoti<-domainvars[!grepl(dmvi,domainvars)]
    for (dmvnotithisi in dmvnoti){
      print(dmvnotithisi)
    comparedfnodomainpairs$othervarsnooutcome[comparedfnodomainpairs$outcome==dmvi]<-gsub(dmvnotithisi,"",comparedfnodomainpairs$othervarsnooutcome[comparedfnodomainpairs$outcome==dmvi])
    checkflipvarlength<-length(grep(dmvnotithisi,comparedfnodomainpairs$othervarsnooutcome_signflip[comparedfnodomainpairs$outcome==dmvi]))
    if(checkflipvarlength >0){
    comparedfnodomainpairs$othervarsnooutcome_signflip[comparedfnodomainpairs$outcome==dmvi]<-gsub(dmvnotithisi,"",comparedfnodomainpairs$othervarsnooutcome_signflip[comparedfnodomainpairs$outcome==dmvi])
    }
    }
  }
}

propbyagesens<-propspecbyage_gam_signflip(df,agevar,comparedfnodomainpairs)
propbyagebyvarsens<-propbyagesens$byvardf
propbyagebyvarsens$dataset<-"NKI"

propbyagebyvarsens<-pivot_longer(propbyagebyvarsens,cols=c(nonspecdev_propofage,specdev_propofage), names_to="type", values_to="prop")
propbyagebyvarsens$typef<-factor(ifelse(propbyagebyvarsens$type=="specdev_propofage","Measure Specific","Common EF"),levels=c("Measure Specific","Common EF"))
propbyagebyvarsens$prop[propbyagebyvarsens$prop<0]<-0
propbyagebyvarsens$prop[propbyagebyvarsens$prop>1]<-1

names(groupslong)[names(groupslong)=="type"]<-"vartype"
propbyagebyvarsens<-merge(propbyagebyvarsens,groupslong,by=c("outcome"))

propbyagebyvarsens$varname<-paste(propbyagebyvarsens$group,propbyagebyvarsens$vartype,sep="_")
propbyagebyvarsens$typebygroup<-paste(propbyagebyvarsens$vartype,propbyagebyvarsens$group,sep="_") ###type first for ordering in plot
propbyagebyvarsens$groupbytypef<-factor(propbyagebyvarsens$varname,levels=unique(propbyagebyvarsens$varname[order(propbyagebyvarsens$typebygroup)]))
propbyagebyvarsens$prop100<-propbyagebyvarsens$prop*100

# gpefvars<-ggplot(propbyagebyvarsens,aes(x=groupbytypef,y=prop100,fill=typef))+geom_bar(position="stack", stat="identity",colour="black",alpha=.95)+
#   ylab("")+xlab("")+theme(legend.title = element_blank())+
#   scale_y_continuous(limits=c(0,110),expand = c(0,0),breaks=c(0,25,50,75,100))+theme(strip.text.x = element_blank())+scale_fill_manual(values=c("grey92","grey28"))+
#   ylab("")
# gpefvars<-LNCDR::lunaize(gpefvars)+theme(legend.title = element_blank(),legend.position = "top")+ylab("% of Age-Related EF")+xlab("")+
#   theme(axis.text.x = element_text(angle = 90,size=15)) ####switched to hashed bar for easier reading
# 

gpefvarspatternsens<-ggplot(propbyagebyvarsens,aes(x=groupbytypef,y=prop100,fill=typef))+geom_bar_pattern(aes(pattern=typef),position="stack", stat="identity",colour="black",alpha=.95,
                                                                                                          pattern_colour  = 'grey55',
                                                                                                          pattern_angle = 45,
                                                                                                          pattern_density = 0.005,
                                                                                                          pattern_spacing = 0.018,
                                                                                                          pattern_key_scale_factor = 0.6)+
  ylab("")+xlab("")+theme(legend.title = element_blank())+
  scale_y_continuous(limits=c(0,110),expand = c(0,0),breaks=c(0,25,50,75,100))+theme(strip.text.x = element_blank())+scale_fill_manual(values=c("grey92","grey28"))+
  ylab("")+
  scale_pattern_manual(values = c("stripe","none"))#+
  #scale_x_discrete(limits = scalexdiscretelimits)
gpefvarspatternsens<-LNCDR::lunaize(gpefvarspatternsens)+theme(legend.title = element_blank(),legend.position = "top")+ylab("% of Age-Related EF")+xlab("")+
  theme(axis.text.x = element_text(angle = 90,size=15))

####savelist######

savelistout<-list(propbyagesens=propbyagesens,gpefvarspatternsens=gpefvarspatternsens,groupslong=groupslong)
save(savelistout,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/Agecomspec/LUNAagecomspecalldata.Rdata")

#######iter######

propbyageitersens<-propspecbyage_gam_signflip_iter(df,agevar,comparedfnodomainpairs)
propbyageitersens_save<-propbyageitersens$byvardf
save(propbyageitersens_save,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/Agecomspec/LUNAiterpropbyage.Rdata")

