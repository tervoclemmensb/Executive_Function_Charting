library(survey)
library(haven)
library(mice)
library(dplyr)
#################################
NKI<-read.csv("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/NKI/Data/btc_NKIscoredmeasures_20220306.outlierremoved.compositeacclat.2030scale.csv")
NKI$id<-as.factor(NKI$subject)
NKI$metaage<-NKI$age
NKI$dataset<-"NKI"
NKI$visitnum<-1 ##cross-sectional

NKIdemo<-read.csv("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_Behavioral/NKI/Data/data_demo_20220916.csv")
NKIdemo$subject<-unlist(lapply(strsplit(NKIdemo$Identifiers,"[,]"),"[[",1))
NKIdemo<-NKIdemo[NKIdemo$subject %in% NKI$subject,]
NKIdemo$eth<-NKIdemo$dem_003
NKIdemo$race<-NKIdemo$dem_004
NKIdemo$sex<-NKIdemo$dem_002
NKIdemoonobs<-NKIdemo %>% group_by(subject) %>% summarize(eth=unique(eth)[1],race=unique(race)[1],sex=unique(sex)[1])

NKIdemosup<-read.csv("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_Behavioral/NKI/Data/data_demo_supp_20220916.csv")
NKIdemosup$subject<-unlist(lapply(strsplit(NKIdemosup$Identifiers,"[,]"),"[[",1))
NKIdemosup<-NKIdemosup[NKIdemosup$subject %in% NKI$subject,]
NKIdemosup$income<-NKIdemosup$demos_02
NKIdemosuponeobs<-NKIdemosup %>% group_by(subject) %>% summarize(income=unique(income)[1])

NKIwithalldem<-merge(NKI,NKIdemoonobs[,c("subject","eth","race","sex")],by="subject") %>% merge(.,NKIdemosuponeobs[,c("subject","income")],by="subject")

###recode to Match ACS#
###race/eth,sex, age,income 
##race/ethnicity
#NKI labels American Indian or Native Alaskan = 1, Asian = 2, Black  or African American = 3, Native Hawiian or Other Pacific Islander = 4, White = 5, Other Race = 6, Don't Know, Missing Data
##[1] "Asian AIAN NHPI OTHER" "Hispanic"              "Non-Hispanic Black/AA" "Non-Hispanic White"    "Other"(ACS levels based on ABCD-style analysis)
NKIwithalldem$raceeth_recode<-NA
NKIwithalldem$raceeth_recode[NKIwithalldem$eth=="1"]<-"Hispanic"
NKIwithalldem$raceeth_recode[NKIwithalldem$eth=="0" & NKIwithalldem$race %in% c("1","2","4")]<-"Asian AIAN NHPI OTHER"
NKIwithalldem$raceeth_recode[NKIwithalldem$eth=="0" & NKIwithalldem$race %in% c("3")]<-"Non-Hispanic Black/AA"
NKIwithalldem$raceeth_recode[NKIwithalldem$eth=="0" & NKIwithalldem$race %in% c("5")]<-"Non-Hispanic White"
NKIwithalldem$raceeth_recode[is.na(NKIwithalldem$raceeth_recode)]<-"Other"
###Sex
#NKI labels Male = 1, Female = 2, Don't Know, Missing Data
##ACS   1   Male 2 Female , matched okay change name
NKIwithalldem$SEX<-as.factor(as.character(NKIwithalldem$sex))
###Age bins (as in https://jamanetwork.com/journals/jamanetworkopen/fullarticle/2764298)
#8-12, 13-17, 18-22,23-27,28-35
NKIwithalldem$age_cat<-NA
NKIwithalldem$age_cat[round(NKIwithalldem$age) %in% seq(8,12,1)]<-"8-12"
NKIwithalldem$age_cat[round(NKIwithalldem$age) %in% seq(13,17,1)]<-"13-17"
NKIwithalldem$age_cat[round(NKIwithalldem$age) %in% seq(18,22,1)]<-"18-22"
NKIwithalldem$age_cat[round(NKIwithalldem$age) %in% seq(23,27,1)]<-"23-27"
NKIwithalldem$age_cat[round(NKIwithalldem$age) %in% seq(28,35,1)]<-"28-35"
######income####
##NKI Less than 10k = 1, 10k to 14.9K = 2, 15k to 24.9K = 3, 25K to 34.9K = 4, 35K to 49.9K = 5, 50K to 74.9K = 6, 75K to 99.9 K = 7, 100K to 199.9K = 8, 200k or more 9, Prefer not to disclose = 10, Don't Know, Missing Data
##ACS levels"<25k"      "100k-199k" "200k+"     "25k-49k"   "50k-74k"   "75k-99k"
NKIwithalldem$income_recode<-NA 
NKIwithalldem$income_recode[NKIwithalldem$income %in% c(1,2,3)]<-"<25k"
NKIwithalldem$income_recode[NKIwithalldem$income %in% c(4,5)]<-"25k-49k"
NKIwithalldem$income_recode[NKIwithalldem$income %in% c(6)]<-"50k-74k"
NKIwithalldem$income_recode[NKIwithalldem$income %in% c(7)]<-"75k-99k"
NKIwithalldem$income_recode[NKIwithalldem$income %in% c(8)]<-"100k-199k"
NKIwithalldem$income_recode[NKIwithalldem$income %in% c(9)]<-"200k+"

NKIwithalldem$income_recode_impute<-NKIwithalldem$income_recode ###match ACS name

NKI_r<-NKIwithalldem[,c("id","income_recode_impute")]
write.csv(NKI_r,"~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_Behavioral/NKI/Data/nkiincome.csv")


NKIwithalldemfortable<-NKIwithalldem[,c("raceeth_recode","SEX","age_cat","income_recode_impute")]
NKIwithalldemfortable$raceeth_recode<-factor(NKIwithalldemfortable$raceeth_recode,levels=c("Non-Hispanic White","Non-Hispanic Black/AA","Asian AIAN NHPI OTHER","Hispanic","Other"))
NKIwithalldemfortable$age_cat<-factor(NKIwithalldemfortable$age_cat,levels=c("8-12","13-17","18-22","23-27","28-35"))
NKIwithalldemfortable$income_recode_impute<-factor(NKIwithalldemfortable$income_recode_impute,levels=c("<25k","25k-49k","50k-74k","75k-99k","100k-199k","200k+"))
demolist<-sapply(NKIwithalldemfortable, function(x) prop.table(table(x)))
demolistdf<-do.call(rbind,lapply(demolist,function(x){data.frame(x)}))
demolistdf$Freq<-demolistdf$Freq*100
write.csv(demolistdf,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_Behavioral/GeneralScripts/RakingSens/NKIdemotable.csv")

NKIwithalldemcomplete<-NKIwithalldem[complete.cases(NKIwithalldem[,c("subject","visit","raceeth_recode","SEX","age_cat","income_recode_impute")]),]
####ACS Load####
load(file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_Behavioral/GeneralScripts/RakingSens/ACSdataages.imputed.Rdata")
dat.2011_15<-allacswitimpute
dat.2011_15$AGE<-as.numeric(as.character(dat.2011_15$AGE))
dat.2011_15$age_cat<-NA
dat.2011_15$age_cat[round(dat.2011_15$AGE) %in% seq(8,12,1)]<-"8-12"
dat.2011_15$age_cat[round(dat.2011_15$AGE) %in% seq(13,17,1)]<-"13-17"
dat.2011_15$age_cat[round(dat.2011_15$AGE) %in% seq(18,22,1)]<-"18-22"
dat.2011_15$age_cat[round(dat.2011_15$AGE) %in% seq(23,27,1)]<-"23-27"
dat.2011_15$age_cat[round(dat.2011_15$AGE) %in% seq(28,35,1)]<-"28-35"

dat.2011_15$dataset<-"ACS"

ACSdemfortable<-dat.2011_15[,c("raceeth_recode","SEX","age_cat","income_recode_impute")]
ACSdemfortable$raceeth_recode<-factor(ACSdemfortable$raceeth_recode,levels=c("Non-Hispanic White","Non-Hispanic Black/AA","Asian AIAN NHPI OTHER","Hispanic","Other"))
ACSdemfortable$age_cat<-factor(ACSdemfortable$age_cat,levels=c("8-12","13-17","18-22","23-27","28-35"))
ACSdemfortable$income_recode_impute<-factor(ACSdemfortable$income_recode_impute,levels=c("<25k","25k-49k","50k-74k","75k-99k","100k-199k","200k+"))
ACSdemolist<-sapply(ACSdemfortable, function(x) prop.table(table(x)))
ACSdemolistdf<-do.call(rbind,lapply(ACSdemolist,function(x){data.frame(x)}))

ACSdemolistdf$Freq<-ACSdemolistdf$Freq*100
write.csv(ACSdemolistdf,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_Behavioral/GeneralScripts/RakingSens/ACSdemotable.csv")

##############
alldatarake<-dplyr::bind_rows(NKIwithalldemcomplete[,c("subject","visit","raceeth_recode","SEX","age_cat","income_recode_impute","dataset")],dat.2011_15)

alldatarake$meth1wgt = alldatarake$PERWT 
alldatarake$meth1wgt [alldatarake$dataset=="NKI"] <- 1
summary(alldatarake$meth1wgt, exclude = NULL)

alldatarake$acsflag<-dplyr::if_else(alldatarake$dataset=="ACS",1,0)
alldatarake[,c("raceeth_recode","SEX","age_cat","income_recode_impute")]<-lapply(alldatarake[,c("raceeth_recode","SEX","age_cat","income_recode_impute")],function(x){as.factor(as.character(x))})


m1fit1 <- glm(relevel(factor(acsflag),ref='1') ~ 1+
              + relevel(factor(income_recode_impute),ref="200k+")###reference level matches ABCD approach
              + relevel(factor(raceeth_recode),ref='Other')
              + relevel(factor(age_cat),ref='28-35')
              + relevel(factor(SEX),ref='2'), 
              data=alldatarake, family=quasibinomial(), weights=meth1wgt)
summary(m1fit1)

# create new variables 
alldatarake$propmeth1 = predict(m1fit1, type="response") 
summary(alldatarake$propmeth1) 

# assign different weight values depending on acsflag status 
alldatarake$pwgtmeth1=1/alldatarake$propmeth1
summary(alldatarake$pwgtmeth1) 

alldatarake$pwgtmeth1 <- ifelse(alldatarake$acsflag %in% (1), alldatarake$PERWT, alldatarake$pwgtmeth1) 
summary(alldatarake$pwgtmeth1)

# load psych package for use of descriptive tools 
library(psych)

describeBy(alldatarake$pwgtmeth1, alldatarake$acsflag)
quantile(alldatarake$pwgtmeth1, c(0,.01,.02,.05,.25,.50,.75, .95,.98,.99,1.0))

# by group of acsflag 
# acsflag==1
acssub <- subset(alldatarake,acsflag==1)
quantile(acssub$pwgtmeth1, c(0,.01,.02,.05,.25,.50,.75, .95,.98,.99,1.0))

# acsflag==0 
NKISub <- subset(alldatarake,acsflag==0)
quantile(NKISub$pwgtmeth1, c(0,.01,.02,.05,.25,.50,.75, .95,.98,.99,1.0))

# Trim weight at .02 and .98 percentiles for acsflag==0 only
alldatarake$pwgtmeth1 [alldatarake$acsflag==0 & alldatarake$pwgtmeth1 < quantile(NKISub$pwgtmeth1,.02)] <- quantile(NKISub$pwgtmeth1,.02)
alldatarake$pwgtmeth1 [alldatarake$acsflag==0 & alldatarake$pwgtmeth1 > quantile(NKISub$pwgtmeth1,.98)] <- quantile(NKISub$pwgtmeth1,.98)

describeBy(alldatarake$pwgtmeth1, alldatarake$acsflag)

# obtain sum of weight for NKI 
sum(NKISub$pwgtmeth1)

# sum of weight for acs 
sum(acssub$pwgtmeth1) 

# Adjust to control total, note that the control value should be changed with each change in data (denominator)  

alldatarake$pwgtmeth1 <- ifelse(alldatarake$acsflag %in% c(0), 
                                       (alldatarake$pwgtmeth1*sum(acssub$pwgtmeth1)/sum(NKISub$pwgtmeth1)), alldatarake$pwgtmeth1) 

describeBy(alldatarake$pwgtmeth1, alldatarake$acsflag) 

library(survey)

# set survey design with single psu for each respondent and weight, no strata or FPC 
svyd <- svydesign(id=~1, weights=alldatarake$pwgtmeth1, data=alldatarake)

# tables of key vars by acsflag, with labels using trimmed and adjusted weight 
# weighted two way crosstab from Survey package using svyby and svymean
#"raceeth_recode","SEX","age_cat","income_recode_impute",

svyby(~age_cat, ~acsflag, design=svyd, svymean)
svyby(~income_recode_impute, ~acsflag, design=svyd, svymean) 
svyby(~SEX, ~acsflag, design=svyd, svymean) 
svyby(~raceeth_recode, ~acsflag, design=svyd, svymean) 

# Regress weight on the propensity model covariates by acsflag

# for acs subset first 
wtfit_acs <- glm(pwgtmeth1 ~ 1+
                 + relevel(factor(income_recode_impute),ref="200k+")###reference level matches ABCD approach
                 + relevel(factor(raceeth_recode),ref='Other')
                 + relevel(factor(age_cat),ref='28-35')
                 + relevel(factor(SEX),ref='2'), 
                 data=acssub)
summary(wtfit_acs)

# for NKI subset 
wtfit_NKI <- glm(pwgtmeth1 ~ 1 
                  + relevel(factor(income_recode_impute),ref="200k+")###reference level matches ABCD approach
                  + relevel(factor(raceeth_recode),ref='Other')
                  + relevel(factor(age_cat),ref='28-35')
                  + relevel(factor(SEX),ref='2'),
                  data=NKISub)
summary(wtfit_NKI)

# create population total marginal data sets for rake of weight variable, age sex and collaped race 
# rake is done for NKI weights using ACS population marginal totals as controls for age sex collapsed race/eth

library(Hmisc)
library(questionr)
# population data sets for raking 
pop.age <- data.frame(age_cat=names(wtd.table(acssub$age_cat,weights=acssub$PERWT)),Freq=as.numeric(wtd.table(acssub$age_cat,weights=acssub$PERWT))) 
pop.sex <- data.frame(SEX=names(wtd.table(acssub$SEX,weights=acssub$PERWT)),Freq=as.numeric(wtd.table(acssub$SEX,weights=acssub$PERWT)))
pop.race <- data.frame(raceeth_recode=names(wtd.table(acssub$raceeth_recode,weights=acssub$PERWT)), Freq=as.numeric(wtd.table(acssub$raceeth_recode,weights=acssub$PERWT))) 
pop.age
pop.sex
pop.race

# survey design data for NKI only, set psu =1 and no strata, with weight 
svyNKI1 <- svydesign(id=~1, strata=NULL, weights=~pwgtmeth1, data=NKISub)

# rake using population data above 
svyNKI_test_3  <- rake(svyNKI1, list(~age_cat, ~SEX, ~raceeth_recode), list(pop.age,pop.sex,pop.race ))

# extract weight 
svyNKI_test_3$rpwgtmeth1 <-(weights(svyNKI_test_3))

# add to abcdsub data frame 
NKISub$rpwgtmeth1 <- svyNKI_test_3$rpwgtmeth1 
names(NKISub)

# check raked totals 
svytable (~age_cat, svyNKI_test_3) 
svytable (~SEX, svyNKI_test_3)
svytable (~raceeth_recode, svyNKI_test_3)

# check with raked weight rpwgtmeth1 
wtd.table(NKISub$age_cat, weights=NKISub$rpwgtmeth1) 
wtd.table(NKISub$SEX, weights=NKISub$rpwgtmeth1)
wtd.table(NKISub$raceeth_recode, weights=NKISub$rpwgtmeth1)

# evaluate raked propensity method 1 weight by regressing on propensity model covariates
wtfit_rake <- glm(rpwgtmeth1 ~ pwgtmeth1 
                  + relevel(factor(income_recode_impute),ref="200k+")###reference level matches ABCD approach
                  + relevel(factor(raceeth_recode),ref='Other')
                  + relevel(factor(age_cat),ref='28-35')
                  + relevel(factor(SEX),ref='2'),
                  data=NKISub)
summary(wtfit_rake)

# weighted estimates of ABCD variables using raked propensity weight rpwgtmeth1 
wtd.table(NKISub$age_cat,weights=NKISub$rpwgtmeth1)
wtd.table(NKISub$SEX,weights=NKISub$rpwgtmeth1) 
wtd.table(NKISub$income_recode_impute,weights=NKISub$rpwgtmeth1) 
wtd.table(NKISub$raceeth_recode,weights=NKISub$rpwgtmeth1) 

# save data for future use
# save data with revised weights including trim values and sum of weight in denominator for controls 
saveRDS(NKISub, file= "~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_Behavioral/GeneralScripts/RakingSens/NKISubwithweights.Rdata")

# Use saved data set for weighted tables after final raked weight is prepared 
# Set svydesign data 

NKIwithalldemandweights<-merge(NKIwithalldem,NKISub[,c("subject","rpwgtmeth1")],by="subject")

NKIwithalldemandweights$ageyear<-round(NKIwithalldemandweights$age)

#######unweighted#######
source("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/GeneralScripts/growthrate_mgcvgam.R")
library(mgcv)
baseformula<-as.formula('outcome~s(pred)')
unweightedGAMfits<-mgcvscalefits(NKIwithalldemandweights,outcomevars = c("Accuracycomposite","Latencycomposite"),predvars = "age",mformula = baseformula,scale=FALSE) ####no additional scaling
unweightedGAMfits$type<-dplyr::if_else(unweightedGAMfits$outcome=="Accuracycomposite","acc","lat")
unweightedGAMfits$wtype<-"unweighted"
#ymin=fit-2*se,ymax=fit+2*se
#################weighted
svyNKId <- svydesign(data=NKIwithalldemandweights, id=~1, strata=NULL, weights=NKIwithalldemandweights$rpwgtmeth1) 
varsforweightestimate<-c("Accuracycomposite","Latencycomposite")
library(splines)
library(emmeans)
weightedfits<-lapply(varsforweightestimate,function(thisoutcome){
  print(thisoutcome)
  svyNKId$variables$thisoutcome<-svyNKId$variables[,thisoutcome]
  #svyweightedmeans<-svyby(~thisoutcome, ~age, design=svyNKId,svymean)
  svymod<-svyglm("thisoutcome~ns(age,4)",svyNKId) ###basis number matched to GAM
  svypred<-data.frame(summary(emmeans(svymod, ~ age,type="response", at=list(age=seq(min(NKIwithalldemandweights$age,na.rm=TRUE),max(NKIwithalldemandweights$age,na.rm=TRUE),.1)))))
  svypred$pred<-svypred$age
  svypred$fit<-svypred$emmean
  svypred$se<-svypred$SE
  svypred$wtype<-"pop. weighted"
  svypred$outcome<-thisoutcome
  return(svypred)
})
weightedfits<-do.call(rbind,weightedfits)
weightedfits$type<-dplyr::if_else(weightedfits$outcome=="Accuracycomposite","acc","lat")
####unweighted and weighted#####
wanduwfits<-dplyr::bind_rows(unweightedGAMfits[,c("pred","fit","se","wtype","outcome","type")],weightedfits[,c("pred","fit","se","wtype","outcome","type")])
wanduwfits$outcomelabel<-unlist(strsplit(wanduwfits$outcome,"composite"))
wanduwfits$wtype<-factor(wanduwfits$wtype,levels=c("unweighted","pop. weighted"))

library(ggplot2)

ggNKIweightedunweighted<-ggplot(wanduwfits,aes(x=pred,y=fit,colour=wtype,fill=wtype))+geom_line(colour="black")+geom_ribbon(aes(ymin=fit-2*se,ymax=fit+2*se,fill=wtype),colour="black",alpha=.5)+
  theme(legend.position = "top")+
  facet_grid(cols=vars(outcomelabel),scales="free")+theme(strip.text.y= element_text(angle = 180))+
  scale_x_continuous(limits = c(8, 36),breaks=(c(10,15,20,25,30,35)))+scale_colour_manual(values=c("grey11","#9ebcda"))+scale_fill_manual(values=c("grey11","#9ebcda"))
ggNKIweightedunweighted<-LNCDR::lunaize(ggNKIweightedunweighted)+xlab("Age (years)")+ylab("Executive Function\n (z-score)")+theme(legend.position = "top",legend.title = element_blank())+
  theme(strip.text.y= element_text(angle = 360))+theme(strip.background = element_blank())

save(ggNKIweightedunweighted,file= "~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_Behavioral/GeneralScripts/RakingSens/ggNKIacclatwithweights.Rdata")
# svyfrm<-as.formula(sprintf("%s~ns(age,%s)",riskvar,round(edfagevar)))

wanduwfits$dataset<-"NKI"
save(wanduwfits,file= "~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_Behavioral/GeneralScripts/RakingSens/ggNKIacclatwithweights.suppdata.Rdata")



