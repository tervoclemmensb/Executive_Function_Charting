library(survey)
library(haven)
library(mice)
library(dplyr)
#################################
PNC<-read.csv("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/PNC/Data/btc_R03cleaneddata_20220306.outlierremoved.compositeacclat.2030scale.csv")
PNC$id<-as.factor(PNC$dbGaP_Subject_ID)
PNC$metaage<-PNC$age
PNC$dataset<-"PNC"
PNC$visitnum<-1 ##cross-sectional

allpheno<-read.table("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/PNC/Data/decrypt/phs000607.v3.pht003445.v3.p2.c1.Neurodevelopmental_Genomics_Subject_Phenotypes.GRU-NPU.txt",header=TRUE,fill=TRUE,sep = '\t')
allpheno_r<-allpheno[,c("dbGaP_Subject_ID","Race","Sex")]
allpheno_r<-allpheno_r[!duplicated(allpheno_r),]
save(allpheno_r,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/PNC/Data/PNCRaceSexinfo.Rdata")
#Race,#Sex

PNCwithalldem<-merge(PNC,allpheno_r[,c("dbGaP_Subject_ID","Race","Sex")],by="dbGaP_Subject_ID")
allpheno<-NULL
demo<-NULL
gc()
###recode to Match ACS#
###race/eth,sex, age,income 
##race/ethnicity
#PNC labels AA=Black or African American; AI=American Indian or Alaska Native; AS=Asian; EA=European American; HI=Hispanic/Latino; PI=Native Hawaiian/Pacific Islander; OT=Other; NA=Not Available/Pending Validation 
##[1] "Asian AIAN NHPI OTHER" "Hispanic"              "Non-Hispanic Black/AA" "Non-Hispanic White"    "Other"(ACS levels based on ABCD-style analysis)
PNCwithalldem$raceeth_recode<-PNCwithalldem$Race
PNCwithalldem$raceeth_recode[PNCwithalldem$Race=="EA"]<-"Non-Hispanic White"
PNCwithalldem$raceeth_recode[grepl("PI|AI|Eskimo/Alaskan",PNCwithalldem$Race)]<-"Asian AIAN NHPI OTHER"
PNCwithalldem$raceeth_recode[PNCwithalldem$Race=="OT"]<-"Other"
PNCwithalldem$raceeth_recode[grepl("AA",PNCwithalldem$Race)]<-"Non-Hispanic Black/AA"
PNCwithalldem$raceeth_recode[grepl("HI",PNCwithalldem$Race)]<-"Hispanic"
PNCwithalldem$raceeth_recode[grepl("EA",PNCwithalldem$raceeth_recode)]<-"Non-Hispanic White" ###remaining EA with OT set to NHW

###Sex
#PNC labels M, F, Recode to 1, 2 to match ACS
##ACS   1   Male 2 Female , matched okay change name
PNCwithalldem$SEX[PNCwithalldem$Sex=="M"]<-1
PNCwithalldem$SEX[PNCwithalldem$Sex=="F"]<-2
PNCwithalldem$SEX<-as.factor(PNCwithalldem$SEX)
###Age bins (as in https://jamanetwork.com/journals/jamanetworkopen/fullarticle/2764298)
#8-12, 13-17, 18-22,23-27,28-35
PNCwithalldem$age_cat<-NA
PNCwithalldem$age_cat[round(PNCwithalldem$age) %in% seq(8,12,1)]<-"8-12"
PNCwithalldem$age_cat[round(PNCwithalldem$age) %in% seq(13,17,1)]<-"13-17"
PNCwithalldem$age_cat[round(PNCwithalldem$age) %in% seq(18,22,1)]<-"18-22"
PNCwithalldem$age_cat[round(PNCwithalldem$age) %in% seq(23,27,1)]<-"23-27"
PNCwithalldem$age_cat[round(PNCwithalldem$age) %in% seq(28,35,1)]<-"28-35"

PNCwithalldemcomplete<-PNCwithalldem[complete.cases(PNCwithalldem[,c("dbGaP_Subject_ID","raceeth_recode","SEX","age_cat")]),]

PNCwithalldemfortable<-PNCwithalldem[,c("raceeth_recode","SEX","age_cat")]
PNCwithalldemfortable$raceeth_recode<-factor(PNCwithalldemfortable$raceeth_recode,levels=c("Non-Hispanic White","Non-Hispanic Black/AA","Asian AIAN NHPI OTHER","Hispanic","Other"))
PNCwithalldemfortable$age_cat<-factor(PNCwithalldemfortable$age_cat,levels=c("8-12","13-17","18-22","23-27","28-35"))
demolist<-sapply(PNCwithalldemfortable, function(x) prop.table(table(x)))
demolistdf<-do.call(rbind,lapply(demolist,function(x){data.frame(x)}))
demolistdf$Freq<-demolistdf$Freq*100
write.csv(demolistdf,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_Behavioral/GeneralScripts/RakingSens/PNCdemotable.csv")

####ACS Load####
load(file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_Behavioral/GeneralScripts/RakingSens/ACSdataages.imputed.Rdata")
dat.2011_15<-allacswitimpute
dat.2011_15$AGE<-as.numeric(as.character(dat.2011_15$AGE))
dat.2011_15<-dat.2011_15[dat.2011_15$AGE>=8 & dat.2011_15$AGE<=22,] ####limit to PNC AGES######
dat.2011_15$age_cat<-NA
dat.2011_15$age_cat[round(dat.2011_15$AGE) %in% seq(8,12,1)]<-"8-12"
dat.2011_15$age_cat[round(dat.2011_15$AGE) %in% seq(13,17,1)]<-"13-17"
dat.2011_15$age_cat[round(dat.2011_15$AGE) %in% seq(18,22,1)]<-"18-22"

dat.2011_15$dataset<-"ACS"


##############
alldatarake<-dplyr::bind_rows(PNCwithalldemcomplete[,c("dbGaP_Subject_ID","raceeth_recode","SEX","age_cat","dataset")],dat.2011_15)

alldatarake$meth1wgt = alldatarake$PERWT 
alldatarake$meth1wgt [alldatarake$dataset=="PNC"] <- 1
summary(alldatarake$meth1wgt, exclude = NULL)

alldatarake$acsflag<-dplyr::if_else(alldatarake$dataset=="ACS",1,0)
alldatarake[,c("raceeth_recode","SEX","age_cat","income_recode_impute")]<-lapply(alldatarake[,c("raceeth_recode","SEX","age_cat","income_recode_impute")],function(x){as.factor(as.character(x))})


m1fit1 <- glm(relevel(factor(acsflag),ref='1') ~ 1+
              + relevel(factor(raceeth_recode),ref='Other')
              + relevel(factor(age_cat),ref='18-22')
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
PNCSub <- subset(alldatarake,acsflag==0)
quantile(PNCSub$pwgtmeth1, c(0,.01,.02,.05,.25,.50,.75, .95,.98,.99,1.0))

# Trim weight at .02 and .98 percentiles for acsflag==0 only
alldatarake$pwgtmeth1 [alldatarake$acsflag==0 & alldatarake$pwgtmeth1 < quantile(PNCSub$pwgtmeth1,.02)] <- quantile(PNCSub$pwgtmeth1,.02)
alldatarake$pwgtmeth1 [alldatarake$acsflag==0 & alldatarake$pwgtmeth1 > quantile(PNCSub$pwgtmeth1,.98)] <- quantile(PNCSub$pwgtmeth1,.98)

describeBy(alldatarake$pwgtmeth1, alldatarake$acsflag)

# obtain sum of weight for PNC 
sum(PNCSub$pwgtmeth1)

# sum of weight for acs 
sum(acssub$pwgtmeth1) 
alldatarake %>% group_by(acsflag) %>% summarize(sum(alldatarake$pwgtmeth1))
# Adjust to control total, note that the control value should be changed with each change in data (denominator)  

alldatarake$pwgtmeth1 <- ifelse(alldatarake$acsflag %in% c(0), 
                                       (alldatarake$pwgtmeth1*sum(acssub$pwgtmeth1)/sum(PNCSub$pwgtmeth1)), alldatarake$pwgtmeth1) 

describeBy(alldatarake$pwgtmeth1, alldatarake$acsflag) 
alldatarake %>% group_by(acsflag) %>% summarize(sum(alldatarake$pwgtmeth1))
library(survey)

# set survey design with single psu for each respondent and weight, no strata or FPC 
svyd <- svydesign(id=~1, weights=alldatarake$pwgtmeth1, data=alldatarake)

# tables of key vars by acsflag, with labels using trimmed and adjusted weight 
# weighted two way crosstab from Survey package using svyby and svymean
#"raceeth_recode","SEX","age_cat","income_recode_impute",

svyby(~age_cat, ~acsflag, design=svyd, svymean)
svyby(~SEX, ~acsflag, design=svyd, svymean) 
svyby(~raceeth_recode, ~acsflag, design=svyd, svymean) 

# Regress weight on the propensity model covariates by acsflag

# for acs subset first 
wtfit_acs <- glm(pwgtmeth1 ~ 1+
                 + relevel(factor(raceeth_recode),ref='Other')###reference level matches ABCD approach
                 + relevel(factor(age_cat),ref='18-22')
                 + relevel(factor(SEX),ref='2'), 
                 data=acssub)
summary(wtfit_acs)

# for PNC subset 
wtfit_PNC <- glm(pwgtmeth1 ~ 1 
                 + relevel(factor(raceeth_recode),ref='Other')###reference level matches ABCD approach
                 + relevel(factor(age_cat),ref='18-22')
                  + relevel(factor(SEX),ref='2'),
                  data=PNCSub)
summary(wtfit_PNC)

# create population total marginal data sets for rake of weight variable, age sex and collaped race 
# rake is done for PNC weights using ACS population marginal totals as controls for age sex collapsed race/eth

library(Hmisc) 
library(questionr)
# population data sets for raking 
pop.age <- data.frame(age_cat=names(wtd.table(acssub$age_cat,weights=acssub$PERWT)),Freq=as.numeric(wtd.table(acssub$age_cat,weights=acssub$PERWT))) 
pop.sex <- data.frame(SEX=names(wtd.table(acssub$SEX,weights=acssub$PERWT)),Freq=as.numeric(wtd.table(acssub$SEX,weights=acssub$PERWT)))
pop.race <- data.frame(raceeth_recode=names(wtd.table(acssub$raceeth_recode,weights=acssub$PERWT)), Freq=as.numeric(wtd.table(acssub$raceeth_recode,weights=acssub$PERWT))) 
pop.age
pop.sex
pop.race

# survey design data for PNC only, set psu =1 and no strata, with weight 
svyPNC1 <- svydesign(id=~1, strata=NULL, weights=~pwgtmeth1, data=PNCSub)

# rake using population data above 
svyPNC_test_3  <- rake(svyPNC1, list(~age_cat, ~SEX, ~raceeth_recode), list(pop.age,pop.sex,pop.race ))

# extract weight 
svyPNC_test_3$rpwgtmeth1 <-(weights(svyPNC_test_3))

# add to abcdsub data frame 
PNCSub$rpwgtmeth1 <- svyPNC_test_3$rpwgtmeth1 
names(PNCSub)

# check raked totals 
svytable (~age_cat, svyPNC_test_3) 
svytable (~SEX, svyPNC_test_3)
svytable (~raceeth_recode, svyPNC_test_3)

# check with raked weight rpwgtmeth1 
wtd.table(PNCSub$age_cat, weights=PNCSub$rpwgtmeth1) 
wtd.table(PNCSub$SEX, weights=PNCSub$rpwgtmeth1)
wtd.table(PNCSub$raceeth_recode, weights=PNCSub$rpwgtmeth1)

# evaluate raked propensity method 1 weight by regressing on propensity model covariates
wtfit_rake <- glm(rpwgtmeth1 ~ pwgtmeth1 
                  + relevel(factor(raceeth_recode),ref='Other')###reference level matches ABCD approach
                  + relevel(factor(age_cat),ref='18-22')
                  + relevel(factor(SEX),ref='2'),
                  data=PNCSub)
summary(wtfit_rake)

# weighted estimates of ABCD variables using raked propensity weight rpwgtmeth1 
wtd.table(PNCSub$age_cat,weights=PNCSub$rpwgtmeth1)
wtd.table(PNCSub$SEX,weights=PNCSub$rpwgtmeth1) 
wtd.table(PNCSub$income_recode_impute,weights=PNCSub$rpwgtmeth1) 
wtd.table(PNCSub$raceeth_recode,weights=PNCSub$rpwgtmeth1) 

# save data for future use
# save data with revised weights including trim values and sum of weight in denominator for controls 
saveRDS(PNCSub, file= "~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_Behavioral/GeneralScripts/RakingSens/PNCSubwithweights.Rdata")

# Use saved data set for weighted tables after final raked weight is prepared 
# Set svydesign data 

PNCwithalldemandweights<-merge(PNCwithalldem,PNCSub[,c("dbGaP_Subject_ID","rpwgtmeth1")],by="dbGaP_Subject_ID")

PNCwithalldemandweights$ageyear<-round(PNCwithalldemandweights$age)
#######unweighted#######
source("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/GeneralScripts/growthrate_mgcvgam.R")
library(mgcv)
baseformula<-as.formula('outcome~s(pred)')
unweightedGAMfits<-mgcvscalefits(PNCwithalldemandweights,outcomevars = c("Accuracycomposite","Latencycomposite"),predvars = "age",mformula = baseformula,scale=FALSE) ####no additional scaling
unweightedGAMfits$type<-dplyr::if_else(unweightedGAMfits$outcome=="Accuracycomposite","acc","lat")
unweightedGAMfits$wtype<-"unweighted"
#ymin=fit-2*se,ymax=fit+2*se
#################weighted
svyPNCd <- svydesign(data=PNCwithalldemandweights, id=~1, strata=NULL, weights=PNCwithalldemandweights$rpwgtmeth1) 
varsforweightestimate<-c("Accuracycomposite","Latencycomposite")
library(splines)
library(emmeans)
weightedfits<-lapply(varsforweightestimate,function(thisoutcome){
  print(thisoutcome)
  svyPNCd$variables$thisoutcome<-svyPNCd$variables[,thisoutcome]
  #svyweightedmeans<-svyby(~thisoutcome, ~age, design=svyPNCd,svymean)
  svymod<-svyglm("thisoutcome~ns(age,4)",svyPNCd) ###basis number matched to GAM
  svypred<-data.frame(summary(emmeans(svymod, ~ age,type="response", at=list(age=seq(min(PNCwithalldemandweights$age,na.rm=TRUE),max(PNCwithalldemandweights$age,na.rm=TRUE),.1)))))
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

ggPNCweightedunweighted<-ggplot(wanduwfits,aes(x=pred,y=fit,colour=wtype,fill=wtype))+geom_line(colour="black")+geom_ribbon(aes(ymin=fit-2*se,ymax=fit+2*se,fill=wtype),colour="black",alpha=.5)+
  theme(legend.position = "top")+
  facet_grid(cols=vars(outcomelabel),scales="free")+theme(strip.text.y= element_text(angle = 180))+
  scale_x_continuous(limits = c(8, 36),breaks=(c(10,15,20,25,30,35)))+scale_colour_manual(values=c("grey11","#9ebcda"))+scale_fill_manual(values=c("grey11","#9ebcda"))
ggPNCweightedunweighted<-LNCDR::lunaize(ggPNCweightedunweighted)+xlab("Age (years)")+ylab("Executive Function\n (z-score)")+theme(legend.position = "top",legend.title = element_blank())+
  theme(strip.text.y= element_text(angle = 360))+theme(strip.background = element_blank())

save(ggPNCweightedunweighted,file= "~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_Behavioral/GeneralScripts/RakingSens/ggPNCacclatwithweights.Rdata")


wanduwfits$dataset<-"PNC"
save(wanduwfits,file= "~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_Behavioral/GeneralScripts/RakingSens/ggPNCacclatwithweights.suppdata.Rdata")



