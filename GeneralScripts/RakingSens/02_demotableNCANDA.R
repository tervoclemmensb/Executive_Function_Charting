library(survey)
library(haven)
library(mice)
library(dplyr)
#################################
NCANDA<-read.csv("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/NCANDA/data/btc_R03cleaneddata_20220306.outlierremoved.compositeacclat.2030scale.csv")
NCANDA$id<-as.factor(NCANDA$subject)
NCANDA$metaage<-NCANDA$cnp_age
NCANDA$visitnum<-as.numeric(as.factor(NCANDA$visit))
NCANDA$dataset<-"NCANDA"

NCANDA<-NCANDA[NCANDA$visit=="baseline",] ###baseline only for demo report 

NCANDAdemo<-read.csv("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/NCANDA_Behavioral/Data/NCANDA_RELEASE_4Y_REDCAP_MEASUREMENTS_V02/summaries/redcap/demographics.csv")
NCANDAdemo<-NCANDAdemo[NCANDAdemo$visit=="baseline",]

NCANDAparentreport<-read.csv("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/NCANDA_Behavioral/Data/NCANDA_RELEASE_4Y_REDCAP_MEASUREMENTS_V02/summaries/redcap/parentreport.csv")
NCANDAparentreport<-NCANDAparentreport[NCANDAparentreport$visit=="baseline",]
NCANDAparentreport$income<-NCANDAparentreport$parentreport_sesp14

NCANDAwithalldem<-merge(NCANDA,NCANDAdemo[,c("subject","hispanic","race","sex")],by="subject") %>% merge(.,NCANDAparentreport[,c("subject","income")],by="subject")

###recode to Match ACS#
###race/eth,sex, age,income 
##race/ethnicity
#NCANDA labels 1, Native American/American Indian | 2, Asian | 3, Pacific Islander | 4, African-American/Black | 5, Caucasian/White | 6, Other (Specify) | 8, Invalid Entry | 9, No Response
##[1] "Asian AIAN NHPI OTHER" "Hispanic"              "Non-Hispanic Black/AA" "Non-Hispanic White"    "Other"(ACS levels based on ABCD-style analysis)
NCANDAwithalldem$raceeth_recode[NCANDAwithalldem$hispanic=="Y"]<-"Hispanic"
NCANDAwithalldem$raceeth_recode[NCANDAwithalldem$hispanic=="N" & NCANDAwithalldem$race %in% c("1","2","3")]<-"Asian AIAN NHPI OTHER"
NCANDAwithalldem$raceeth_recode[NCANDAwithalldem$hispanic=="N" & NCANDAwithalldem$race %in% c("4")]<-"Non-Hispanic Black/AA"
NCANDAwithalldem$raceeth_recode[NCANDAwithalldem$hispanic=="N" & NCANDAwithalldem$race %in% c("5")]<-"Non-Hispanic White"
NCANDAwithalldem$raceeth_recode[is.na(NCANDAwithalldem$raceeth_recode)]<-"Other"
###Sex
#NCANDA labels Male = M, Female = F
##ACS   1   Male 2 Female
NCANDAwithalldem$SEX<-as.factor(dplyr::if_else(NCANDAwithalldem$sex=="M",1,2))
###Age bins (as in https://jamanetwork.com/journals/jamanetworkopen/fullarticle/2764298)
NCANDAwithalldem$age_cat<-NA
NCANDAwithalldem$age_cat[round(NCANDAwithalldem$cnp_age) %in% seq(12,17,1)]<-"13-17"
NCANDAwithalldem$age_cat[round(NCANDAwithalldem$cnp_age) %in% seq(18,22,1)]<-"18-22"
######income####
##NCANDA 1, Less than $5,000 | 2, $5,000 through $11,999 | 3, $12,000 through $15,999 | 4, $16,000 through $24,999 | 5, $25,000 through $34,999 | 6, $35,000 through $49,999 | 7, $50,000 through $74,999 | 8, $75,000 through $99,999 | 9, $100,000 through $199,999 | 10, $200,000 and greater | 11, Don't know
##ACS levels"<25k"      "100k-199k" "200k+"     "25k-49k"   "50k-74k"   "75k-99k"
NCANDAwithalldem$income_recode<-NA 
NCANDAwithalldem$income_recode[NCANDAwithalldem$income %in% c(1,2,3,4)]<-"<25k"
NCANDAwithalldem$income_recode[NCANDAwithalldem$income %in% c(5,6)]<-"25k-49k"
NCANDAwithalldem$income_recode[NCANDAwithalldem$income %in% c(7)]<-"50k-74k"
NCANDAwithalldem$income_recode[NCANDAwithalldem$income %in% c(8)]<-"75k-99k"
NCANDAwithalldem$income_recode[NCANDAwithalldem$income %in% c(9)]<-"100k-199k"
NCANDAwithalldem$income_recode[NCANDAwithalldem$income %in% c(10)]<-"200k+"

NCANDAwithalldem$income_recode_impute<-NCANDAwithalldem$income_recode ###match ACS name

NCANDA_r<-NCANDAwithalldem[,c("id","income_recode_impute","income")]
write.csv(NCANDA_r,"~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_Behavioral/NCANDA/Data/ncandaincome.csv")

NCANDAwithalldemfortable<-NCANDAwithalldem[,c("raceeth_recode","SEX","age_cat","income_recode_impute")]
NCANDAwithalldemfortable$raceeth_recode<-factor(NCANDAwithalldemfortable$raceeth_recode,levels=c("Non-Hispanic White","Non-Hispanic Black/AA","Asian AIAN NHPI OTHER","Hispanic","Other"))
NCANDAwithalldemfortable$age_cat<-factor(NCANDAwithalldemfortable$age_cat,levels=c("8-12","13-17","18-22","23-27","28-35"))
NCANDAwithalldemfortable$income_recode_impute<-factor(NCANDAwithalldemfortable$income_recode_impute,levels=c("<25k","25k-49k","50k-74k","75k-99k","100k-199k","200k+"))
demolist<-sapply(NCANDAwithalldemfortable, function(x) prop.table(table(x)))
demolistdf<-do.call(rbind,lapply(demolist,function(x){data.frame(x)}))
demolistdf$Freq<-demolistdf$Freq*100
write.csv(demolistdf,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_Behavioral/GeneralScripts/RakingSens/NCANDAdemotable.csv")




