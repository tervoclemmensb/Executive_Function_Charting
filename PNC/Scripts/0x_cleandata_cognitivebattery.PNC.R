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
#########CNP##########
allpheno<-read.table("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/PNC/Data/decrypt/phs000607.v3.pht003445.v3.p2.c1.Neurodevelopmental_Genomics_Subject_Phenotypes.GRU-NPU.txt",header=TRUE,fill=TRUE,sep = '\t')
demo<-read.table("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/PNC/Data/decrypt/phs000607.v3.pht007881.v1.p2.c1.Interview_Dates_Update.GRU-NPU.txt",header=TRUE,fill=TRUE)
allphenodemo<-merge(demo,allpheno,by=c("dbGaP_Subject_ID","SUBJID"))
hisubjects<-allpheno[grep("HI",allpheno$Race),]
hisubjectsCNB<-hisubjects[,c("dbGaP_Subject_ID",grep("PCET|PCPT|LNB",names(hisubjects),value=TRUE))]

write.csv(allphenodemo[,c("dbGaP_Subject_ID","Sex")],"~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/PNC/Data/PNCsex.csv")
write.csv(allphenodemo[,c("dbGaP_Subject_ID","Mother_Education","Father_Education")],"~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/PNC/Data/PNCparentaled.csv")
write.csv(allphenodemo[,c("dbGaP_Subject_ID","PVRT_CR")],"~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/PNC/Data/PNCPVRT.csv")
write.csv(allphenodemo[,c("dbGaP_Subject_ID","SCR007","SCR006")],"~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/PNC/Data/PNCPSY.csv")
#https://ftp.ncbi.nlm.nih.gov/dbgap/studies/phs000607/phs000607.v3.p2/pheno_variable_summaries/phs000607.v3.pht003445.v3.p2.Neurodevelopmental_Genomics_Subject_Phenotypes.var_report.xml
###PADT #Age differentation test
####PFMT penn face memory test
####PEIT penn emotion indentification test
####PWMT penn word memory test ##ef 
###PVRT penn verbal reasoning test
###PEDT penn emotion differntiation test
###MP motor praxis
###PMAT Penn matrix reasoning test  #EF
####TAP: missing for now#####
#####VOLT Visual object learning test 
####LNB Letter N-back  ##ef
###PCET conditional exclusion test ##ef 
###PCPT continous performance test
###PLOT penn line orientation test
###WRAT Wide range assessement test

allcogtestsacc<-c("PADT_A","PFMT_IFAC_TOT","PEIT_CR","PWMT_KIWRD_TOT","PVRT_CR","PEDT_A","PMAT_CR","VOLT_SVT","LNB_TP","PCET_ACC2","PCPT_T_TP","PLOT_TC","WRAT_CR_RAW")
allcogtestslatency<-c("PADT_T","PFMT_IFAC_RTC","PEIT_CRT","PWMT_KIWRD_RTC","PVRT_RTCR","PEDT_T","MP_MP2RTCR","PMAT_RTCR","VOLT_SVTCRT","LNB_RTC","PCET_RTCR","PCPT_T_TPRT","PLOT_TCRT")

allcogvars<-c(allcogtestsacc,allcogtestslatency)

testgroups<-c("PADT","PFMT","PEIT","PWMT","PVRT","PEDT","PMAT","VOLT","LNB","PCET","PCPT","PLOT","WRAT","MP")
efficiencytestgroups<-c("PADT","PFMT","PEIT","PWMT","PVRT","PEDT","PMAT","VOLT","LNB","PCPT","PLOT")
eftestgroups<-c("PCET","PCPT","LNB") ###EXECUTIVE-CONTROL#HTTPS://WWW.NCBI.NLM.NIH.GOV/PMC/ARTICLES/PMC3295891/PDF/NIHMS348849.PDF
#eftestgroupsplusPVRT<-c("")
varsvalidgenusclean<-function(df,subjvar="dbGaP_Subject_ID",agevar="ageAtCnb",vargroup){
  allvarsingroup<-grep(vargroup,names(df),value=TRUE)
  datavars<-allvarsingroup[!grepl("VALID|GENUS",allvarsingroup)]
  vargroupdf<-df[,c(allvarsingroup,subjvar,agevar)]
  if (vargroup=="PMAT"){
    vargroupdf$PMAT_PC<-NULL
    datavars<-datavars[datavars!="PMAT_PC"]
  }
  vargroupdf<-vargroupdf[complete.cases(vargroupdf),]
  genus<-grep("GENUS",allvarsingroup,value=TRUE)
  if (length(genus)!=0){
    gsummary<-summary(vargroupdf[,genus])
    print(gsummary)
    usegenus<- readline(prompt="limit to most common genus?[y/n]")
    if (usegenus=="y"){
      commongenus<-names(gsummary)[which.max(gsummary)]
      vargroupdf[(vargroupdf[,genus]!=commongenus),datavars]<-NA
    }
    if (usegenus=="n"){
    multiplegenus<-readline(prompt="use more than one genus?[n/genus string (e.g., PMAT24)]")
    if (multiplegenus!="n"){
      genera<-unique(grep(multiplegenus,vargroupdf[,genus],value=TRUE))
    vargroupdf[!(vargroupdf[,genus] %in% genera),datavars]<-NA
    }
    }
  }
  vargroupdf<-vargroupdf[complete.cases(vargroupdf),]
  valid<-grep("VALID",allvarsingroup,value=TRUE)
  if (length(valid)!=0){
    print(unique(vargroupdf[,valid]))
    print(summary(vargroupdf[,valid]))
    print("setting non-valid rows to NA")
    vargroupdf[(vargroupdf[,valid]!="V"),datavars]<-NA
  }
  return(vargroupdf)
}

allvarsgvclean<-lapply(eftestgroups,function(vg){
  print(vg)
  vargccleaned<-varsvalidgenusclean(allphenodemo,vargroup=vg)
})
require(purrr)
require(dplyr)
allvarsgvcleanmat<-allvarsgvclean %>% purrr::reduce(dplyr::left_join,by=c("dbGaP_Subject_ID","ageAtCnb"))

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
  accvar<-grep(vg,allcogtestsacc,value=TRUE)
  nacc<-length(which(!is.na(allvarsgvcleanmat[,accvar])))
  latvar=grep(vg,allcogtestslatency,value=TRUE)
  nlat<-length(which(!is.na(allvarsgvcleanmat[,latvar])))
  print(sprintf("Latency has %s useable datapoints; Accuracy has %s usable datapoints",nlat,nacc))
  eff<-efficiency_function(allvarsgvcleanmat,accvar=accvar,latvar=latvar)
  effdf<-as.data.frame(eff)
  names(effdf)<-paste0(vg,"_eff")
  return(effdf)
})

alleffs<-as.data.frame(do.call(cbind,allefficiency))
allcogtestsaccwitheffs<-cbind(allvarsgvcleanmat,alleffs)

allcogtestsaccwitheffs$age<-allcogtestsaccwitheffs$ageAtCnb/12
load(file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/PNC/Data/PNCRaceSexinfo.Rdata")
allcogtestsaccwitheffswithdem<-merge(allcogtestsaccwitheffs,allpheno_r,by="dbGaP_Subject_ID")

write.csv(allcogtestsaccwitheffs,"~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/PNC/Data/btc_PNCscoredmeasures_20220214.csv")
#allcogtestsaccwitheffs<-read.csv("/Users/brendenclemmens/Desktop/Projects/R03_behavioral/PNC/Data/btc_PNCscoredmeasures_20200605.csv")
