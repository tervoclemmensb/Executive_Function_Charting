ages<-data.frame(linearage=seq(8,35,by=.1))
ages$age<-ages$linearage
ages$invage<-1/ages$linearage

ages_long<-ages %>% tidyr::pivot_longer(cols=c("linearage","invage"),names_to="agetype",values_to="value")
ages_long<-ages_long %>% group_by(agetype) %>% mutate(diff = value - lag(value))



derivexample<-ggplot(ages_long,aes(x=age,y=diff))+geom_smooth()+facet_wrap(~agetype,scales="free")+ylab("local age-related change (first deriv)")+
  scale_x_continuous(breaks=seq(10,100,by=5),limits=c(7,100))+xlab("Age")
derivexample

ages<-data.frame(linearage=seq(8,50,by=.1))
ages$age<-ages$linearage
ages$invage<-1/ages$linearage

ages_long<-ages %>% tidyr::pivot_longer(cols=c("linearage","invage"),names_to="agetype",values_to="value")
ages_long<-ages_long %>% group_by(agetype) %>% mutate(diff = value - lag(value))



derivexample<-ggplot(ages_long,aes(x=age,y=diff))+geom_smooth()+facet_wrap(~agetype,scales="free")+ylab("local age-related change (first deriv)")+
  scale_x_continuous(breaks=seq(10,100,by=5),limits=c(7,100))+xlab("Age")
derivexample


ages<-data.frame(linearage=seq(8,100,by=.1))
ages$age<-ages$linearage
ages$invage<-1/ages$linearage

ages_long<-ages %>% tidyr::pivot_longer(cols=c("linearage","invage"),names_to="agetype",values_to="value")
ages_long<-ages_long %>% group_by(agetype) %>% mutate(diff = value - lag(value))



derivexample<-ggplot(ages_long,aes(x=age,y=diff))+geom_smooth()+facet_wrap(~agetype,scales="free")+ylab("local age-related change (first deriv)")+
  scale_x_continuous(breaks=seq(10,100,by=5),limits=c(7,100))+xlab("Age")
derivexample

ggclosedfit<-ggplot(ages,aes(x=age,y=invage))+geom_line(size=1)+ylab("Inverse Age Fit")+xlab("Age")
ggclosedfit<-LNCDR::lunaize(ggclosedfit)+scale_x_continuous(limits=c(8,35))

ggclosedformderiv<-ggplot(ages,aes(x=age,y=-1/age^2))+geom_line(size=1)+ylab("Derivative of 1/Age")+xlab("Age")
ggclosedformderiv<-LNCDR::lunaize(ggclosedformderiv)+scale_x_continuous(limits=c(8,35))

require(patchwork)

fitwithderivativedev<-ggclosedfit/ggclosedformderiv

ggclosedfit<-ggplot(ages,aes(x=age,y=invage))+geom_line(size=1)+ylab("Inverse Age Fit")+xlab("Age")
ggclosedfit<-LNCDR::lunaize(ggclosedfit)

ggclosedformderiv<-ggplot(ages,aes(x=age,y=-1/age^2))+geom_line(size=1)+ylab("Derivative of 1/Age")+xlab("Age")
ggclosedformderiv<-LNCDR::lunaize(ggclosedformderiv)

fitwithderivativedevfull<-ggclosedfit/ggclosedformderiv

#########
fitwithderivativedev<-ggclosedfit/ggclosedformderiv

ggclosedfit<-ggplot(ages,aes(x=age,y=invage))+geom_line(size=1)+ylab("Inverse Age Fit")+xlab("Age")
ggclosedfit<-LNCDR::lunaize(ggclosedfit)+scale_y_continuous(limits=c(0,.2))

ggclosedfit10<-ggplot(ages,aes(x=age,y=invage))+geom_line(size=1)+ylab("Inverse Age Fit")+xlab("Age")
ggclosedfit10<-LNCDR::lunaize(ggclosedfit10)+scale_x_continuous(limits=c(10,25))

ggclosedfit14<-ggplot(ages,aes(x=age,y=invage))+geom_line(size=1)+ylab("Inverse Age Fit")+xlab("Age")
ggclosedfit14<-LNCDR::lunaize(ggclosedfit14)+scale_x_continuous(limits=c(14,25))

ggclosedfit18<-ggplot(ages,aes(x=age,y=invage))+geom_line(size=1)+ylab("Inverse Age Fit")+xlab("Age")
ggclosedfit18<-LNCDR::lunaize(ggclosedfit18)+scale_x_continuous(limits=c(18,25))

fitallageszoom<-ggclosedfit10+ggclosedfit14+ggclosedfit18




ggclosedformderiv<-ggplot(ages,aes(x=age,y=-1/age^2))+geom_line(size=1)+ylab("Derivative of 1/Age")+xlab("Age")
ggclosedformderiv<-LNCDR::lunaize(ggclosedformderiv)





################
ages$invagebigeffect<-ages$invage*100
ggclosedfit<-ggplot(ages,aes(x=age,y=invagebigeffect))+geom_line(size=1)+ylab("Inverse Age Fit")+xlab("Age")
ggclosedfit<-LNCDR::lunaize(ggclosedfit)

ggclosedformderiv<-ggplot(ages,aes(x=age,y=-1/(age*100)^2))+geom_line(size=1)+ylab("Derivative of 1/Age")+xlab("Age")
ggclosedformderiv<-LNCDR::lunaize(ggclosedformderiv)

fitwithderivativedevfull<-ggclosedfit/ggclosedformderiv




