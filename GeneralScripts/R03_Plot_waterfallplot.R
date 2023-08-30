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
library(RColorBrewer)
library(patchwork)
library(gridExtra)
library(grid)
########################
mycolors<-brewer.pal(n=9,"Blues")
mycolors<-c("grey72",mycolors)
mycolors[length(mycolors)]<-"#253494"
######################
###LUNA#########
coglongdata<-read.csv("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/data/btc_R03cleaneddata_20220306.outlierremoved.compositeacclat.csv")
coglongdata$id<-as.factor(coglongdata$id)
coglongdata$id<-as.factor(coglongdata$id)
coglongdata$age<-coglongdata$Ageatvisit

visitnumsbyparticipant<-coglongdata %>% dplyr::group_by(id) %>% dplyr::summarize(visitnumtotal=max(visitnum))
range(visitnumsbyparticipant$visitnumtotal)
median(visitnumsbyparticipant$visitnumtotal)
#####################################
age_rankedluna <- LNCDR::waterfall_group(coglongdata)

p <- ggplot(age_rankedluna) + aes(x = age, y = age_id, group = age_id) + geom_line(size=.2) + geom_point()
p <- lunaize(p) + ylab("Participants") + xlab("Age") + theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())+
  theme(text = element_text(size=32))

lunap<-p
ggsave(lunap,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/Waterfallplots/lunawaterfallplot.pdf",height=7,width=7)


coglongdata$visitnumplot<-factor(round(coglongdata$visitnum), levels=max(round(coglongdata$visitnum)):1)

gghist <- ggplot(coglongdata) + aes(x = age,colour=visitnumplot,fill=visitnumplot) + geom_histogram(colour="black")+scale_color_manual(values=c(mycolors))+scale_fill_manual(values=c(mycolors))+
  coord_cartesian(xlim =c(8,35))
legend <- cowplot::get_legend(gghist)
grid.newpage()
grid.draw(legend)

gghist <- lunaize(gghist) + ylab("Visit Count") + xlab("Age(years)") +theme(text = element_text(size=26))+theme(legend.position = "none")+scale_y_continuous(limits = c(0,55), expand = c(0, 0))
ggsave(gghist,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/Waterfallplots/lunawaterfallplot.hist.pdf",height=7,width=7)

# 
# p <- ggplot(age_rankedluna) + aes(x = age, y = age_id, group = age_id,colour=Gender) + geom_line() + geom_point()
# p <- lunaize(p) + ylab("") + xlab("Age") + theme(axis.title.y = element_blank(), 
#                                                  axis.ticks.y = element_blank(), axis.text.y = element_blank())+
#   scale_colour_manual(values=c("white","#F8766D","#00BFC4"))+theme(legend.title = element_blank(),legend.position="top") 


Luna_coglongdata<-coglongdata[,c("age","visitnumplot")]
Luna_coglongdata$dataset<-"Luna"
########################################
############NCANDA######################
coglongdata<-read.csv("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/NCANDA/Data/btc_NCANDAscoredmeasures_20220306.outlierremoved.compositeacclat.csv")
coglongdata$id<-as.factor(coglongdata$subject)
coglongdata$age<-coglongdata$cnp_age
coglongdata$visitnum<-as.numeric(as.factor(coglongdata$visit))
visitnumsbyparticipant<-coglongdata %>% dplyr::group_by(id) %>% dplyr::summarize(visitnumtotal=max(visitnum))
range(visitnumsbyparticipant$visitnumtotal)
median(visitnumsbyparticipant$visitnumtotal)

age_ranked<-waterfall_group(coglongdata)

ncandap <- ggplot(age_ranked) + aes(x = age, y = age_id, group = age_id) + geom_line(size=.2) + geom_point()+theme(legend.position = "none")
ncandap <- lunaize(ncandap) + ylab("Participants") + xlab("Age") + theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())+
  theme(text = element_text(size=32))+theme(legend.position = "none")

ggsave(ncandap,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/Waterfallplots/NCANDAwaterfallplot.pdf",height=7,width=7)


coglongdata$visitnumplot<-factor(round(coglongdata$visitnum), levels=max(round(coglongdata$visitnum)):1)

gghistNCANDA <- ggplot(coglongdata) + aes(x = age,colour=visitnumplot,fill=visitnumplot) +geom_histogram(colour="black",bins=25)+scale_color_manual(values=c(rev(mycolors[10:6])))+scale_fill_manual(values=c(rev(mycolors[10:6])))+
  coord_cartesian(xlim =c(8,35))

gghistNCANDA <- lunaize(gghistNCANDA) + ylab("Visit Count") + xlab("Age (years)") +theme(text = element_text(size=26))+theme(legend.position = "none")+scale_y_continuous(limits = c(0,325),expand = c(0, 0))

ggsave(gghistNCANDA,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/Waterfallplots/NCANDAwaterfallplot.hist.pdf",height=7,width=7)

NCANDA_coglongdata<-coglongdata[,c("age","visitnumplot")]
NCANDA_coglongdata$dataset<-"NCANDA"

######NKI########
coglongdata<-read.csv("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/NKI/Data/btc_NKIscoredmeasures_20220306.outlierremoved.compositeacclat.csv")
coglongdata$id<-as.factor(coglongdata$subject)

age_ranked<-waterfall_group(coglongdata)

NKIp <- ggplot(age_ranked) + aes(x = age, y = age_id, group = age_id) + geom_line(size=.2) + geom_point()+theme(legend.position = "none")
NKIp <- lunaize(NKIp) + ylab("Participants") + xlab("Age") + theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())+
  theme(text = element_text(size=32))+theme(legend.position = "none")

ggsave(NKIp,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/Waterfallplots/NCANDAwaterfallplot.pdf",height=7,width=7)

coglongdata$visitnumplot<-factor(round(coglongdata$visitnum), levels=max(round(coglongdata$visitnum)):1)

gghistNKI <- ggplot(coglongdata) + aes(x = age,colour=visitnumplot,fill=visitnumplot) +geom_histogram(colour="black",bins=25)+scale_color_manual(values=c(rev(mycolors[10:6])))+scale_fill_manual(values=c(rev(mycolors[10])))+
  coord_cartesian(xlim =c(8,35))

gghistNKI <- lunaize(gghistNKI) + ylab("Visit Count") + xlab("Age (years)") +theme(text = element_text(size=26))+theme(legend.position = "none")+scale_y_continuous(limits = c(0,45),expand = c(0, 0))
ggsave(gghistNKI,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/Waterfallplots/NKIwaterfallplot.hist.pdf",height=7,width=7)

NKI_coglongdata<-coglongdata[,c("age","visitnumplot")]
NKI_coglongdata$dataset<-"NKI"
########################################
############PNC######################
coglongdata<-read.csv("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/PNC/Data/btc_PNCscoredmeasures_20220306.outlierremoved.compositeacclat.csv")
coglongdata$id<-coglongdata$dbGaP_Subject_ID
age_rankedPNC<-waterfall_group(coglongdata)

PNCp <- ggplot(age_rankedPNC) + aes(x = age, y = age_id, group = age_id) + geom_line(size=.2) + geom_point()+theme(legend.position = "none")
PNCp <- lunaize(PNCp) + ylab("Participants") + xlab("Age") + theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())+
  theme(text = element_text(size=32))+theme(legend.position = "none")

coglongdata$visitnumplot<-factor(1)

gghistPNC <- ggplot(coglongdata) + aes(x = age,colour=visitnumplot,fill=visitnumplot) + geom_histogram(colour="black")+scale_color_manual(values=c(rev(mycolors[10])))+scale_fill_manual(values=c(rev(mycolors[10])))+
  coord_cartesian(xlim =c(8,35))

gghistPNC <- lunaize(gghistPNC) + ylab("Visit Count") + xlab("Age (years)") +theme(text = element_text(size=26))+theme(legend.position = "none")+scale_y_continuous(limits = c(0,500), expand = c(0, 0))

library(patchwork)
allhists<-(gghist)/ plot_spacer()/ (gghistNCANDA)/plot_spacer()/(gghistNKI)/plot_spacer()/(gghistPNC) + plot_layout(heights = c(6, 2.5 ,6,2.5,6,2.5,6))
ggsave(allhists,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/Waterfallplots/allsampleswatefall.hist.pdf",height=18,width=7)

allhistquads<-(gghist+gghistNCANDA)/plot_spacer()/(gghistNKI+gghistPNC)+plot_layout(heights = c(6, 1.75 ,6))
ggsave(allhistquads,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/Waterfallplots/allsampleswatefall.quadhist.pdf",height=12,width=14)

PNC_coglongdata<-coglongdata[,c("age","visitnumplot")]
PNC_coglongdata$dataset<-"PNC"
###all sup data####
Luna_coglongdata
NCANDA_coglongdata
NKI_coglongdata
PNC_coglongdata

allcoghistdata<-rbind(Luna_coglongdata,NCANDA_coglongdata) %>% rbind(.,NKI_coglongdata) %>% rbind(.,PNC_coglongdata)
write.csv(allcoghistdata,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Data/SupportingData/Sup1.csv")
