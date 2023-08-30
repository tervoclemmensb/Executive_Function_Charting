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
library(lavaan)
############load all data#######
load(file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/Factorfigs/FAbyagebin.Luna.Rdata")
Lunafactoragebin<-FAbyagebindf
Lunafactoragebin$dataset<-"Luna"
#####NCANDA########
load(file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/Factorfigs/FAbyagebin.NCANDA.Rdata")
NCANDAfactoragebin<-FAbyagebindf
NCANDAfactoragebin$dataset<-"NCANDA"
######NKI#######
load(file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/Factorfigs/FAbyagebin.NKI.Rdata")
NKIfactoragebin<-FAbyagebindf
NKIfactoragebin$dataset<-"NKI"
######PNC###########
load(file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/Factorfigs/FAbyagebin.PNC.Rdata")
PNCfactoragebin<-FAbyagebindf
PNCfactoragebin$dataset<-"PNC"

#####alldata######

allfactoragebin<-rbind(Lunafactoragebin,NCANDAfactoragebin) %>% rbind(.,NKIfactoragebin) %>% rbind(.,PNCfactoragebin)

###plot####

ggfactorbyagebin<-ggplot(allfactoragebin,aes(x=agebin,y=ML1propvarboot,group=agebin))+geom_boxplot(colour="grey55")+
  geom_smooth(method = "gam",formula=y~s(x,k=4),se=FALSE,colour="grey28",aes(group=1))+
  scale_x_continuous(limits = c(8, 36),breaks=(c(10,15,20,25,30,35)))+facet_wrap(~ dataset,nrow = 2)
ggfactorbyagebin<-LNCDR::lunaize(ggfactorbyagebin)+xlab("Age (years)")+ylab("% Total EF Variance Explained\n Domain General (Factor 1)")

write.csv(allfactoragebin,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Data/SupportingData/Sup6.csv")


  



