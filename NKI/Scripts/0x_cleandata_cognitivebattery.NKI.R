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
#####0X PNC cognitive data#######
#########NKI CNP##########
NKIcnp<-read.csv("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/NKI/Data/NKICNP_longdata-2022-02-06T02_59_23.382Z.csv",header=TRUE,fill=TRUE,sep=",")
NKIcnp$SUBJID<-unlist(lapply(strsplit(NKIcnp$Identifiers,","),"[[",1))
NKIcnp$visit<-unlist(lapply(strsplit(NKIcnp$Identifiers,","),"[[",2))

NKIdkefs<-read.csv("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/NKI/Data/NKIDKEFS_longdata-2022-02-06T21_53_47.361Z.csv",header=TRUE,fill=TRUE,sep=",")
NKIdkefs$SUBJID<-unlist(lapply(strsplit(NKIdkefs$Identifiers,","),"[[",1))
NKIdkefs$visit<-unlist(lapply(strsplit(NKIdkefs$Identifiers,","),"[[",2))
NKIdkefs$Identifiers<-NULL

demo<-read.table("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/NKI/Data/ASR_YSR_demo_age_long.csv",header=TRUE,fill=TRUE,sep=",")
demo$SUBJID<-unlist(lapply(strsplit(demo$Identifiers,","),"[[",1))
demo$visit<-unlist(lapply(strsplit(demo$Identifiers,","),"[[",2))


####NKI WASI WIAT for sensitivty
NKIwasiwiat<-read.csv("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/NKI/Data/WASI_WIAT-2022-02-12T00_28_53.609Z.csv")
NKIwasiwiat$SUBJID<-unlist(lapply(strsplit(NKIwasiwiat$Identifiers,","),"[[",1))
NKIwasiwiat$visit<-unlist(lapply(strsplit(NKIwasiwiat$Identifiers,","),"[[",2))
NKIwasiwiat$WASIBDraw<-NKIwasiwiat$int_01
NKIwasiwiat$WASIVCraw<-NKIwasiwiat$int_03
NKIwasiwiat$WASIMRraw<-NKIwasiwiat$int_05
NKIwasiwiat$WASISIMraw<-NKIwasiwiat$int_07
NKIwasiwiat$WASIFSIQ<-NKIwasiwiat$int_14
NKIwasiwiat$Identifiers<-NULL
###data dictionary#####
#https://ftp.ncbi.nlm.nih.gov/dbgap/studies/phs000607/phs000607.v3.p2/pheno_variable_summaries/phs000607.v3.pht003445.v3.p2.Neurodevelopmental_Genomics_Subject_Phenotypes.var_report.xml
##NKI var definiton https://www.nature.com/articles/s41598-021-03782-y#Sec2
datadict<-read.csv("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/NKI/Data/NKI_RS_CODEBOOK.csv")
cnpvars<-datadict[grep("Penn CNP",datadict$COINS_instrument),]
##labels####
allcnpefvars<-c("PCET_PCET_ACC2","PCET_PCETRTCR","SLNB2_SLNB_MCR","SLNB2_SLNB_MRTC","SPCPTNL_SCPT_TP","SPCPTNL_SCPT_TPRT")
loriscnpvars<-cnpvars$loris_variable[cnpvars$label %in% allcnpefvars]
for (li in loriscnpvars){
  print(li)
  varout<-cnpvars$label[cnpvars$loris_variable==li]
  print(varout)
  NKIcnp[,varout]<-NKIcnp[,li]
}

cnpaccvars<-c("PCET_PCET_ACC2","SLNB2_SLNB_MCR","SPCPTNL_SCPT_TP")
cnplatvars<-c("PCET_PCETRTCR","SLNB2_SLNB_MRTC","SPCPTNL_SCPT_TPRT")

##DKEFS
NKIdkefs[,c("dkefscwi_09","dkefscwi_12")]<-lapply(NKIdkefs[,c("dkefscwi_09","dkefscwi_12")],function(x){as.numeric(as.character(x))})
NKIdkefs$CWIlat<-rowMeans(NKIdkefs[,c('dkefscwi_09', 'dkefscwi_12')], na.rm=TRUE)
#DKEFS Color Word Interference	dkefs_color_word_interference	dkefscwi_09	dkefscwi_09	Inhibition-Total Time to Complete	Inhibition-Total Time to Complete- (seconds)			
#DKEFS Color Word Interference	dkefs_color_word_interference	dkefscwi_12	dkefscwi_12	Inhibition/Switching-Total Time to Complete	Inhibition/Switching-Total Time to Complete- (seconds)			
#no accuracy due to clinical standards + floor of errors

#DKEFS TOWER	dkefs_tower	tower_46	tow_46	Total Achievement Score Total Raw	Total Achievement Score Raw- Item Achievement Scores Summed Across All Items				#no latency due to clinical standards + floor
NKIdkefs$TOWacc<-NKIdkefs$tow_46
#DKEFS Design Fluency	dkefs_design_fluency	df_03	df_03	Switching Total Correct	Switching Total Correct- Count Correct				
#no latency released
NKIdkefs$DFLacc<-NKIdkefs$df_03
#DKEFS Sorting	dkefs_sorting	dkefssort_42	dkefssort_42	Confirmed Correct Sorts Total Raw score	Confirmed Correct Sorts Total Raw score- Number  Initial Target Sorts with at Least 1 Correct Description Score for the Sort				
#no latency released
NKIdkefs$SORTacc<-NKIdkefs$dkefssort_42

#DKEFS Trails	dkefs_trails	dkefstmt_04	dkefstmt_04	Number-Letter Switching	Number-Letter Switching- Completion Time (seconds)				
#no accuracy due to clinical standards + ceiling 
NKIdkefs$TMTlat<-NKIdkefs$dkefstmt_04

#############
allphenodemo<-merge(demo,NKIcnp,by=c("SUBJID","visit"),all.y=TRUE) %>% merge(.,NKIdkefs,by=c("SUBJID","visit"),all.y=TRUE) %>% merge(.,NKIwasiwiat,by=c("SUBJID","visit"),all.y=TRUE)
allphenodemo<-allphenodemo[!is.na(allphenodemo$SUBJID),]
allphenodemo<-allphenodemo[as.numeric(as.character(allphenodemo$age_04)) <= 35 & as.numeric(as.character(allphenodemo$age_04)) >= 8,]
allphenodemo<-allphenodemo[!is.na(allphenodemo$SUBJID),]


accvars<-c("PCET_PCET_ACC2","SLNB2_SLNB_MCR","SPCPTNL_SCPT_TP","TOWacc","DFLacc")
latvars<-c("PCET_PCETRTCR","SLNB2_SLNB_MRTC","SPCPTNL_SCPT_TPRT","CWIlat","TMTlat")

allefvars<-c(accvars,latvars)

eftestgroups<-c("PCET_PCET","SLNB2_SLNB","SPCPTNL_SCPT") 

##all numeric####
allphenodemo[,c(allefvars,"age_04")]<-lapply(allphenodemo[,c(allefvars,"age_04")],function(x){as.numeric(as.character(x))})

##########

efficiency_function<-function(df,accvar,latvar){
  print("accuracy and latency correlated at:")
  acclatcor<-cor(as.numeric(as.character(df[,accvar])),as.numeric(as.character(df[,latvar])),use="complete")
  print(acclatcor)
  corrsign<-sign(acclatcor)
  print(sprintf("efficiency based on %s sign",corrsign))
  eff<-scale(as.numeric(as.character(df[,accvar])),center=TRUE,scale=TRUE)+scale(corrsign*as.numeric(as.character(df[,latvar])),center=TRUE,scale=TRUE) ####lat sign flipped so higher means "better" on both scales
  return(eff)
} 

allefficiency<-lapply(eftestgroups,function(vg){
  print(vg)
  accvar<-grep(vg,accvars,value=TRUE)
  nacc<-length(which(!is.na(allphenodemo[,accvar])))
  latvar=grep(vg,latvars,value=TRUE)
  nlat<-length(which(!is.na(allphenodemo[,latvar])))
  print(sprintf("Latency has %s useable datapoints; Accuracy has %s usable datapoints",nlat,nacc))
  eff<-efficiency_function(allphenodemo,accvar=accvar,latvar=latvar)
  effdf<-as.data.frame(eff)
  names(effdf)<-paste0(vg,"_eff")
  return(effdf)
})

alleffs<-as.data.frame(do.call(cbind,allefficiency))
allcogtestsaccwitheffs<-cbind(allphenodemo,alleffs)

allcogtestsaccwitheffs$age<-allcogtestsaccwitheffs$age_04 
allcogtestsaccwitheffs$subject<-allcogtestsaccwitheffs$SUBJID

allcogtestsaccwitheffs<-allcogtestsaccwitheffs %>% dplyr::group_by(subject) %>% dplyr::mutate(visitnum=rank(age))
allcogtestsaccwitheffs<-allcogtestsaccwitheffs[allcogtestsaccwitheffs$visitnum==1,] ###only 10 longitudinal visits, treating as cross-sectional RE poor coverage of longitudinal effects
write.csv(allcogtestsaccwitheffs,"~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/NKI/Data/btc_NKIscoredmeasures_20220208.csv")
