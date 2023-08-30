library(survey)
library(haven)
library(mice)
library(dplyr)
#################################
LUNA<-read.csv("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/data/btc_R03cleaneddata_20220306.outlierremoved.compositeacclat.2030scale.csv")
LUNA$metaage<-LUNA$cnp_age
LUNA$visitnum<-as.numeric(as.factor(LUNA$visit))
LUNA$dataset<-"LUNA"

LUNA<-LUNA[LUNA$visit==1,] ###baseline only for demo report 

Lunethalldemo<-read.csv("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Data/cog_demo/cog_demo.csv")
incomecode<-data.frame(codresponse=c(09,08,07,06,05,04,-8,-9),labels=c(250000,175000,87500,62500,37500,25000,NA,NA)) ###midpoint code for family income sum
incomecodelist <- list("9"=250000,"8"=175000,"7"=87500, "6"=62500,"5"=37500,"4"=25000,"-8"=NA,"-9"=NA)
incometemp<-Lunethalldemo[,grep("level.occ",names(Lunethalldemo))]
Lunethalldemoincome<-Lunethalldemo %>% mutate(value.occf=incomecodelist[as.character(level.occf)],value.occm=incomecodelist[as.character(level.occm)],
                                              value.occs=incomecodelist[as.character(level.occs)],value.occsp=incomecodelist[as.character(level.occsp)])

Lunaincome_r<-Lunethalldemoincome[,c("id","vdate",grep("level.occ|value.occ",names(Lunethalldemoincome),value=TRUE))]
Lunaincome_r<-Lunaincome_r[Lunaincome_r$id %in% LUNA$id,]
Lunaincome_r$vdate.d<-as.Date(Lunaincome_r$vdate)
Lunaincome_rorder<-Lunaincome_r %>% arrange(id,vdate.d) %>% group_by(id) %>% mutate(visitno=rank(vdate.d))
Lunaincome_first<-Lunaincome_rorder[Lunaincome_rorder$visitno==1,]

Lunaincome_firstwithrest<-merge(LUNA[,c("id","Ageatvisit")],Lunaincome_first,by="id",all.x=TRUE)
Lunaincome_firstwithrest$value.occf[sapply(Lunaincome_firstwithrest$value.occf, is.null)]<-NA
Lunaincome_firstwithrest$value.occm[sapply(Lunaincome_firstwithrest$value.occm, is.null)]<-NA
Lunaincome_firstwithrest$value.occs[sapply(Lunaincome_firstwithrest$value.occs, is.null)]<-NA
Lunaincome_firstwithrest$value.occsp[sapply(Lunaincome_firstwithrest$value.occsp, is.null)]<-NA
Lunaincome_firstwithrest[,c("value.occf","value.occm","value.occs","value.occsp")]<-lapply(Lunaincome_firstwithrest[,c("value.occf","value.occm","value.occs","value.occsp")],function(x){unlist(x)})
Lunaincome_firstwithrest$value.fam<-rowSums(Lunaincome_firstwithrest[,c("value.occf","value.occm")],na.rm=TRUE)
Lunaincome_firstwithrest$value.fam[Lunaincome_firstwithrest$value.fam==0]<-NA ###missing both set to NA

Lunethracedemo<-read.csv("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Data/cog_demo/cog_ethcollapsed_only.csv")
Lunethracedemo<-Lunethracedemo %>% group_by(id) %>% summarize(raceth=unique(eth)[1])

Lunaalledemohisp<-read.csv("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Data/cog_demo/20230123_hisp_luna/hisp_endorsed.csv",sep=",")
Lunaalledemohisp<-Lunaalledemohisp[!duplicated(Lunaalledemohisp),]
Lunaalledemohisp_r<-Lunaalledemohisp %>% group_by(id) %>% summarize(anyhsiplength=length(which(endorsed_hisp==TRUE)))
Lunaalledemohisp_r$Hispanicfac<-dplyr::if_else(Lunaalledemohisp_r$anyhsiplength>0,"Hispanic","not Hispanic")


LUNAsex<-read.csv("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Data/gender_20220214.csv")
LUNAsex$sex[LUNAsex$Gender=="Female"]<-"F" ###coded as "gender" but choices were limited to male female; framed as sex
LUNAsex$sex[LUNAsex$Gender=="Male"]<-"M"
LUNAsexadd<-read.csv("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Data/gender_20220214_manualadd.csv",sep=" ",header=FALSE)
names(LUNAsexadd)<-c("id","X","Y","sex","Z")
LUNAsex<-plyr::rbind.fill(LUNAsex,LUNAsexadd[c("id","sex")])
LUNAsex<-LUNAsex %>% dplyr::group_by(id) %>% dplyr::summarise(sex=unique(sex),sexn=length(unique(sex)))
LUNAsex<-LUNAsex[LUNAsex$id %in% unique(LUNA$id),]
LUNAsex<-LUNAsex[LUNAsex$sexn==1,] ###reported in table but left out for demographics

Lunawithalldem<-merge(LUNA[,c("id","Ageatvisit")],Lunethracedemo,by="id",all.x=TRUE) %>% merge(.,LUNAsex[,c("id","sex")],by="id",all.x=TRUE) %>%
  merge(.,Lunaalledemohisp_r[,c("id","Hispanicfac")],all.x=TRUE) %>% merge(.,Lunaincome_firstwithrest[c("id","value.fam")],all.x=TRUE)

###recode to Match ACS#
###race/eth,sex, age,income 
##race/ethnicity
Lunawithalldem$raceeth_recode[Lunawithalldem$Hispanicfac=="Hispanic"]<-"Hispanic"
Lunawithalldem$raceeth_recode[Lunawithalldem$Hispanicfac=="not Hispanic" & grepl("Black",Lunawithalldem$raceth) ]<-"Non-Hispanic Black/AA"
Lunawithalldem$raceeth_recode[Lunawithalldem$Hispanicfac=="not Hispanic" & grepl("Asian|HawaiianPacIsl|AmerIndianAlaskan",Lunawithalldem$raceth)]<-"Asian AIAN NHPI OTHER"
Lunawithalldem$raceeth_recode[Lunawithalldem$Hispanicfac=="not Hispanic" & Lunawithalldem$raceth=="White"]<-"Non-Hispanic White"
#Lunawithalldem$raceeth_recode[is.na(Lunawithalldem$raceeth_recode)]<-"Other"
###Sex
#Luna labels Male = M, Female = F
##ACS   1   Male 2 Female
Lunawithalldem$SEX<-as.factor(dplyr::if_else(Lunawithalldem$sex=="M",1,2))
###Age bins (as in https://jamanetwork.com/journals/jamanetworkopen/fullarticle/2764298)
Lunawithalldem$age_cat<-NA
Lunawithalldem$age_cat[round(Lunawithalldem$Ageatvisit) %in% seq(8,12,1)]<-"8-12"
Lunawithalldem$age_cat[round(Lunawithalldem$Ageatvisit) %in% seq(13,17,1)]<-"13-17"
Lunawithalldem$age_cat[round(Lunawithalldem$Ageatvisit) %in% seq(18,22,1)]<-"18-22"
Lunawithalldem$age_cat[round(Lunawithalldem$Ageatvisit) %in% seq(23,27,1)]<-"23-27"
Lunawithalldem$age_cat[round(Lunawithalldem$Ageatvisit) %in% seq(28,35,1)]<-"28-35"

Lunawithalldem$income_recode<-NA 
Lunawithalldem$income_recode[Lunawithalldem$value.fam < 25000 ]<-"<25k"
Lunawithalldem$income_recode[Lunawithalldem$value.fam >= 25000 & Lunawithalldem$value.fam < 50000]<-"25k-49k"
Lunawithalldem$income_recode[Lunawithalldem$value.fam >= 50000 & Lunawithalldem$value.fam < 75000]<-"50k-74k"
Lunawithalldem$income_recode[Lunawithalldem$value.fam >= 75000 & Lunawithalldem$value.fam < 100000]<-"75k-99k"
Lunawithalldem$income_recode[Lunawithalldem$value.fam >= 100000 & Lunawithalldem$value.fam < 200000]<-"100k-199k"
Lunawithalldem$income_recode[Lunawithalldem$value.fam >= 200000]<-"200k+"

Lunawithalldem$income_recode_impute<-Lunawithalldem$income_recode ###match ACS name

###income save####
Luna_r<-Lunawithalldem[,c("id","income_recode_impute","value.fam")]
write.csv(Luna_r,"~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_Behavioral/Data/lunaincome.csv")
#############
Lunawithalldemfortable<-Lunawithalldem[,c("raceeth_recode","SEX","age_cat","income_recode_impute")]
Lunawithalldemfortable$raceeth_recode<-factor(Lunawithalldemfortable$raceeth_recode,levels=c("Non-Hispanic White","Non-Hispanic Black/AA","Asian AIAN NHPI OTHER","Hispanic","Other"))
Lunawithalldemfortable$age_cat<-factor(Lunawithalldemfortable$age_cat,levels=c("8-12","13-17","18-22","23-27","28-35"))
Lunawithalldemfortable$income_recode_impute<-factor(Lunawithalldemfortable$income_recode_impute,levels=c("<25k","25k-49k","50k-74k","75k-99k","100k-199k","200k+"))
demolist<-sapply(Lunawithalldemfortable, function(x) prop.table(table(x)))
demolistdf<-do.call(rbind,lapply(demolist,function(x){data.frame(x)}))
demolistdf$Freq<-demolistdf$Freq*100
write.csv(demolistdf,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_Behavioral/GeneralScripts/RakingSens/Lunademotable.csv")




