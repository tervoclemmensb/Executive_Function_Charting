library(survey)
library(haven)
library(mice)
library(dplyr)
######adapted from https://github.com/ABCD-STUDY/abcd_acs_raked_propensity
load(file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_Behavioral/GeneralScripts/RakingSens/ACSdataages.Rdata")
dat.2011_15<-data.frame(dataages)
dat.2011_15<-haven::zap_labels(dat.2011_15)
nrow(dat.2011_15)

dat.2011_15<-dat.2011_15[dat.2011_15$AGE >= floor(min(NKIwithalldem$age)) & dat.2011_15$AGE <= floor(max(NKIwithalldem$age)),]
dat.2011_15$raceeth_recode[dat.2011_15$RACE==1 & dat.2011_15$HISPAN==0]<-"Non-Hispanic White"
dat.2011_15$raceeth_recode[dat.2011_15$RACE==2 & dat.2011_15$HISPAN==0]<-"Non-Hispanic Black/AA"
dat.2011_15$raceeth_recode[dat.2011_15$HISPAN==1]<-"Hispanic"
dat.2011_15$raceeth_recode[dat.2011_15$RACE %in% c(3,4,5,6)]<-"Asian AIAN NHPI OTHER"
dat.2011_15$raceeth_recode[is.na(dat.2011_15$raceeth_recode)]<-"Other"

dat.2011_15$FTOTINC[dat.2011_15$FTOTINC==9999999]<-NA
dat.2011_15$income_recode<-NA
dat.2011_15$income_recode[dat.2011_15$FTOTINC < 25000]<-"<25k"
dat.2011_15$income_recode[dat.2011_15$FTOTINC >= 25000 & dat.2011_15$FTOTINC <= 49999]<-"25k-49k"
dat.2011_15$income_recode[dat.2011_15$FTOTINC >= 50000 & dat.2011_15$FTOTINC <= 74999]<-"50k-74k"
dat.2011_15$income_recode[dat.2011_15$FTOTINC >= 75000 & dat.2011_15$FTOTINC <= 99999]<-"75k-99k"
dat.2011_15$income_recode[dat.2011_15$FTOTINC >= 100000 & dat.2011_15$FTOTINC <= 199999]<-"100k-199k"
dat.2011_15$income_recode[dat.2011_15$FTOTINC >= 200000]<-"200k+"

dat.2011_15$famtype_recode<-"Other Family"
dat.2011_15$famtype_recode[dat.2011_15$HHTYPE==1]<-"Married"
###imputation as in ABCD approach######
###convert to factor##
imputevars<-c("SEX","AGE","raceeth_recode","income_recode","REGION","famtype_recode","FAMSIZE")
dat.2011_15[,imputevars]<-lapply(dat.2011_15[,imputevars],function(x){as.factor(x)})
###imputation as in ABCD approach##
tmpimp<-dat.2011_15[,c("AGE","REGION","FAMSIZE","income_recode")]
mice::md.pattern(tmpimp)
summary(tmpimp)

init = mice(tmpimp, maxit=0) 
meth = init$method
meth
tmpimpm1 <- mice(tmpimp, method=meth, maxit=1,m=1, seed=0123) 
summary(tmpimpm1)

#extract the first imputed data set 
longm1 <- complete(tmpimpm1, 1) 
summary(longm1) 

longm1$income_recode_impute<-longm1$income_recode

allacswitimpute<-cbind(dat.2011_15,longm1$income_recode_impute)
names(allacswitimpute)[names(allacswitimpute)=="longm1$income_recode_impute"]<-"income_recode_impute"

save(allacswitimpute,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_Behavioral/GeneralScripts/RakingSens/ACSdataages.imputed.Rdata")


