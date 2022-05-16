require(dplyr)
require(metafor)
require(broom)

interpolatebyage<-function(df,agevar,vars,extrainterp=NULL,agegrid=NULL,agegridvar=NULL,interp_age_range=NULL,
                           interval_inc=.1,longformat=FALSE,datanamefromlongformat=NULL,varnameforlongformat=NULL,
                           limitinterptodatarange=TRUE,namemodifier=NULL,returnlongformat=FALSE){
  if (!is.null(agegrid) && limitinterptodatarange){
    agerangefromdf<-range(df[,agevar])
    agegrid<-data.frame(agegrid[agegrid[,agegridvar] >= agerangefromdf[1] & agegrid[,agegridvar] <= agerangefromdf[2],])
    names(agegrid)<-agegridvar
  }
  if (longformat){
    longformatdf<-df %>%
      pivot_wider(names_from = {{varnameforlongformat}}, values_from = {{datanamefromlongformat}}) %>%
      distinct()
    df<-data.frame(longformatdf)
  }
  if(is.null(interp_age_range) && is.null(agegrid)){
    interp_age_range<-range(df[,agevar],na.rm=TRUE)
    agegrid<-data.frame(ages=as.numeric(seq(interp_age_range[1],interp_age_range[2],by=interval_inc)))
    agegridvar<-names(agegrid)
  }
  
  agegrid[,vars]<-lapply(vars,function(vari){
    print(vari)
    dfinterp<-df[!is.na(df[,vari]),]
    approx(x=df[,agevar], y = df[,vari], xout=agegrid[,agegridvar], method="linear")[[2]]
  })
  
  if(returnlongformat){
    agegrid<-agegrid %>% pivot_longer(cols=vars,names_to = "var",values_to="value")
  }
  if (!is.null(namemodifier)){
    if(returnlongformat){
    names(agegrid)[names(agegrid)=="value"]<-paste(namemodifier,"value",sep="_")
    }else{
    names(agegrid)[names(agegrid) %in% vars]<-paste(namemodifier,vars,sep="_")
  }}
  
  return(agegrid)
}
metabyage<-function(interpolatedagedfs,agevar,valcol,secol,datasetvar,verbose=FALSE){
  interpolatedagedfs$valcol<-interpolatedagedfs[,valcol]
  interpolatedagedfs$se<-interpolatedagedfs[,secol]
  interpolatedagedfs$ages<-interpolatedagedfs[,agevar]
  metabyageout<-lapply(unique(interpolatedagedfs[,agevar]),function(ai){
    if(verbose){print(ai)}
    aidata<-interpolatedagedfs[round(interpolatedagedfs$ages,1)==ai,]
    metaout<-data.frame(age=ai,ndatsets=length(unique(aidata$dataset)),valcol=valcol,secol=secol)
    if (length(unique(aidata$dataset))>1){
      m_metareturn<-tryCatch(
        {
          m_metareturnm<-metafor::rma(yi=valcol,sei=se, data = aidata,method="ML")
          broom::tidy(m_metareturnm)
        },
        error = function(e){
          NA
        }
      )
      metaout<-cbind(metaout,m_metareturn)
    }
    return(metaout)
  })
  metabyageoutdf<-do.call(bind_rows,metabyageout)
}

metabyagethreelevel<-function(interpolatedagedfs,agevar,valcol,secol,variablenestvar,datasetvar,verbose=FALSE){
  #https://bookdown.org/MathiasHarrer/Doing_Meta_Analysis_in_R/multilevel-ma.html
  #https://www.tqmp.org/RegularArticles/vol12-3/p154/p154.pdf
  interpolatedagedfs$valcol<-interpolatedagedfs[,valcol]
  interpolatedagedfs$se<-interpolatedagedfs[,secol]
  interpolatedagedfs$sesqv<-interpolatedagedfs[,secol]^2
  interpolatedagedfs$ages<-interpolatedagedfs[,agevar]
  interpolatedagedfs$variablenestvar<-interpolatedagedfs[,variablenestvar]
  metabyageout<-lapply(unique(interpolatedagedfs[,agevar]),function(ai){
    if(verbose){print(ai)}
    aidata<-interpolatedagedfs[round(interpolatedagedfs$ages,1)==round(ai,1),]
    metaout<-data.frame(age=ai,ndatsets=length(unique(aidata$dataset)),valcol=valcol,secol=secol,variablenestvar=variablenestvar)
    if (length(unique(aidata$dataset))>1){
      m_multi<-NULL
      m_metareturn<-tryCatch(
        {
          m_multi<-metafor::rma.mv(valcol, V=sesqv, random = list(~ 1 | variablenestvar, ~ 1 | dataset), data = aidata,method="ML")
          broom::tidy(m_multi)
        },
        error = function(e){
          "model non converge"
        }
      )
      metaout<-cbind(metaout,m_metareturn)
    }
    return(metaout)
  })
  metabyageoutdf<-do.call(bind_rows,metabyageout)
}

scalebyidx<-function(var,idx){
  varscale<-(var-mean(var[idx],na.rm=TRUE))/sd(var[idx],na.rm=TRUE)
}

scalebyidxrenamereturndf<-function(df,vars,idx,appendname=NULL){
  allvarout<-lapply(vars,function(vari){
    print(vari)
    variscale<-data.frame(scalebyidx(df[,vari],idx))
    names(variscale)<-paste(vari,appendname,sep="_")
    return(variscale)
  })
  allvaroutdf<-do.call(cbind,allvarout)
  df<-cbind(df,allvarout)
  return(df)
}

compositecols<-function(cols,data,mincols=2){
  #return composite score (z scored) of data[,cols]
  data_z<-scale(data[,cols])
  compositeout<-scale(base::rowMeans(data_z,na.rm=TRUE))
  NAfilter<-which(rowSums(!is.na(data[,cols])) < mincols) ###minimum no na cols for composite
  compositeout[NAfilter]<-NA
  return(compositeout)
} 

corSEfromRandN<-function(r,n){
  SE=sqrt((1-r^2)/(n-2))
  return(SE)
}

propspecbyage_gam<-function(df,agevar,comparedf,byeachvar=FALSE){
  ####non sign flip version 
  ###comparedf expects primarvar (measure xi) and othervarsnooutcome (meausures [!=xi]) for common var estimation
  byvar<-lapply(1:nrow(comparedf),function(ci){
    print(comparedf$outcome[ci])
    othervars<-unlist(lapply(strsplit(comparedf$othervarsnooutcome[ci]," "), function(x){x[!x ==""]}))
    df$compositei<-compositecols(c(othervars),df)
    if(length(which(!is.na(unique(df[,comparedf$outcome[ci]])))) < 10){
      print("reducing edf to match max dim of outcome")
      edf<-length(which(!is.na(unique(df[,comparedf$outcome[ci]])))) 
      modbaseage<-mgcv::gam(as.formula(sprintf("%s~s(%s,k=%s)",agevar,comparedf$outcome[ci],edf)),data=df) 
      devmodelbaseage<-summary(modbaseage)$dev.expl
      
      modcomp<-mgcv::gam(as.formula(sprintf("%s~s(%s)",agevar,"compositei")),data=df)
      devmodelcomp<-summary(modcomp)$dev.expl
      
      fullmod<-mgcv::gam(as.formula(sprintf("%s~s(%s,k=%s)+s(%s)",agevar,comparedf$outcome[ci],edf,"compositei")),data=df)
      devmodelfull<-summary(fullmod)$dev.expl
    }else{
    modbaseage<-mgcv::gam(as.formula(sprintf("%s~s(%s)",agevar,comparedf$outcome[ci])),data=df) 
    devmodelbaseage<-summary(modbaseage)$dev.expl
    
    modcomp<-mgcv::gam(as.formula(sprintf("%s~s(%s)",agevar,"compositei")),data=df)
    devmodelcomp<-summary(modcomp)$dev.expl
    
    fullmod<-mgcv::gam(as.formula(sprintf("%s~s(%s)+s(%s)",agevar,comparedf$outcome[ci],"compositei")),data=df)
    devmodelfull<-summary(fullmod)$dev.expl
    }
    
    devianceCIboot<-mgcvgam_bootfromformsdeviance(mods=list(fullmod,modcomp,modbaseage),dataforboot=df)
    
    
    formtosinglchar<-function(x){paste(as.character(formula(x)),collapse=",")}
    outdata<-data.frame(outcome=comparedf$outcome[ci],devmodelbaseage=devmodelbaseage,ageform=formtosinglchar(modbaseage),
                        devmodelcomp=devmodelcomp,compform=formtosinglchar(modcomp),othervars=paste(othervars,collapse=","),
                        devmodelfull=devmodelfull,fullform=formtosinglchar(fullmod),
                        sefullminusmodcomp=devianceCIboot$bootalldf_CIwide$form1minusform2_se_fromqnt,
                        secomp=devianceCIboot$bootalldf_CIwide$form2_se_fromqnt,
                        seagebase=devianceCIboot$bootalldf_CIwide$form3_se_fromqnt)
    return(outdata)
  })
  byvardf<-do.call(rbind,byvar)
  byvardf$specdev<-byvardf$devmodelfull-byvardf$devmodelcomp
  byvardf$allprop<-byvardf$devmodelfull+byvardf$specdev
  byvardf$devmodelfull_prop<-byvardf$devmodelfull/byvardf$allprop
  byvardf$specdev_prop<-byvardf$specdev/byvardf$allprop
  
  byvardf$specdev_propofage<-byvardf$specdev/byvardf$devmodelbaseage
  byvardf$nonspecdev_propofage<-1-byvardf$specdev_propofage
  byvardf$nonspecdev<-byvardf$devmodelbaseage-byvardf$specdev

  byeachothervardf<-NA
  if(byeachvar){byeachothervar<-lapply(1:nrow(comparedf),function(ci){
    othervars<-unlist(lapply(strsplit(comparedf$othervarsnooutcome[ci]," "), function(x){x[!x ==""]}))
    othervarloop<-lapply(unique(othervars),function(othervari){
      modbaseage<-mgcv::gam(as.formula(sprintf("%s~s(%s)",agevar,comparedf$outcome[ci])),data=df)
      devmodelbaseage<-summary(modbaseage)$dev.expl
      
      
      modothervari<-mgcv::gam(as.formula(sprintf("%s~s(%s)",agevar,othervari)),data=df)
      devmodelothervari<-summary(modothervari)$dev.expl

      fullmodwithothervar<-mgcv::gam(as.formula(sprintf("%s~s(%s)+s(%s)",agevar,comparedf$outcome[ci],othervari)),data=df)
      devmodelfullwithothervar<-summary(fullmodwithothervar)$dev.expl
      
      formtosinglchar<-function(x){paste(as.character(formula(x)),collapse=",")}
      outdataothervar<-data.frame(outcome=comparedf$outcome[ci],devmodelbaseage=devmodelbaseage,ageform=formtosinglchar(modbaseage),
                          devmodelothervari=devmodelothervari,othervars=othervari,othervarform=formtosinglchar(modothervari),
                          devmodelfullwithothervar=devmodelfullwithothervar,fullform=formtosinglchar(fullmodwithothervar))
      return(outdataothervar)
    })
    othervarloopdf<-do.call(rbind,othervarloop)
  })
  byeachothervardf<-do.call(rbind,byeachothervar)
  
  byeachothervardf$specdev<-byeachothervardf$devmodelfullwithothervar-byeachothervardf$devmodelothervari
  byeachothervardf$allprop<-byeachothervardf$devmodelfullwithothervar+byeachothervardf$specdev
  byeachothervardf$devmodelfull_prop<-byeachothervardf$devmodelfullwithothervar/byeachothervardf$allprop
  byeachothervardf$specdev_prop<-byeachothervardf$specdev/byeachothervardf$allprop
  
  byeachothervardf$specdev_propofage<-byeachothervardf$specdev/byeachothervardf$devmodelbaseage
  byeachothervardf$nonspecdev_propofage<-1-byeachothervardf$specdev_propofage
  byeachothervardf$nonspecdev<-byeachothervardf$devmodelbaseage-byeachothervardf$specdev
  }
  allreturn<-list(byvardf=byvardf,byeachothervardf=byeachothervardf)
  return(allreturn)
}


propspecbyage_gam_signflip<-function(df,agevar,comparedf,byeachvar=FALSE){
  ###comparedf expects primarvar (measure xi) and othervarsnooutcome (meausures [!=xi]) for common var estimation
  ###comparedf expects primarvar (measure xi) and othervarsnooutcome_signflip (meausures [!=xi], acc/lat for opposite signflipped) for common var estimation
  byvar<-lapply(1:nrow(comparedf),function(ci){
    print(comparedf$outcome[ci])
    othervars<-unlist(lapply(strsplit(comparedf$othervarsnooutcome[ci]," "), function(x){x[!x ==""]}))
    othervarssignflip<-unlist(lapply(strsplit(comparedf$othervarsnooutcome_signflip[ci]," "), function(x){x[!x ==""]}))
    dfforcomposite<-df
    dfforcomposite[othervarssignflip]<-lapply(dfforcomposite[othervarssignflip],function(x){x*-1}) ###flip sign for sign flips
    df$compositei<-compositecols(c(othervars,othervarssignflip),dfforcomposite)
    edf<-10####not edf limited for certain cx versions
    if(length(which(!is.na(unique(df[,comparedf$outcome[ci]])))) < 10){
      print("reducing edf to match max dim of outcome")
      edf<-length(which(!is.na(unique(df[,comparedf$outcome[ci]])))) 
      modbaseage<-mgcv::gam(as.formula(sprintf("%s~s(%s,k=%s)",agevar,comparedf$outcome[ci],edf)),data=df) 
      devmodelbaseage<-summary(modbaseage)$dev.expl
      
      modcomp<-mgcv::gam(as.formula(sprintf("%s~s(%s)",agevar,"compositei")),data=df)
      devmodelcomp<-summary(modcomp)$dev.expl
      
      fullmod<-mgcv::gam(as.formula(sprintf("%s~s(%s,k=%s)+s(%s)",agevar,comparedf$outcome[ci],edf,"compositei")),data=df)
      devmodelfull<-summary(fullmod)$dev.expl
    }else{
      modbaseage<-mgcv::gam(as.formula(sprintf("%s~s(%s)",agevar,comparedf$outcome[ci])),data=df) 
      devmodelbaseage<-summary(modbaseage)$dev.expl
      
      modcomp<-mgcv::gam(as.formula(sprintf("%s~s(%s)",agevar,"compositei")),data=df)
      devmodelcomp<-summary(modcomp)$dev.expl
      
      fullmod<-mgcv::gam(as.formula(sprintf("%s~s(%s)+s(%s)",agevar,comparedf$outcome[ci],"compositei")),data=df)
      devmodelfull<-summary(fullmod)$dev.expl
    }
    
    devianceCIboot<-mgcvgam_bootfromformsdeviance(mods=list(fullmod,modcomp,modbaseage),dataforboot=df)
    
    
    formtosinglchar<-function(x){paste(as.character(formula(x)),collapse=",")}
    outdata<-data.frame(outcome=comparedf$outcome[ci],devmodelbaseage=devmodelbaseage,ageform=formtosinglchar(modbaseage),
                        devmodelcomp=devmodelcomp,compform=formtosinglchar(modcomp),othervars=paste(othervars,collapse=","),othervarssignflip=paste(othervarssignflip,collapse=","),
                        devmodelfull=devmodelfull,fullform=formtosinglchar(fullmod),
                        sefullminusmodcomp=devianceCIboot$bootalldf_CIwide$form1minusform2_se_fromqnt,
                        secomp=devianceCIboot$bootalldf_CIwide$form2_se_fromqnt,
                        seagebase=devianceCIboot$bootalldf_CIwide$form3_se_fromqnt)
    return(outdata)
  })
  byvardf<-do.call(rbind,byvar)
  byvardf$specdev<-byvardf$devmodelfull-byvardf$devmodelcomp
  byvardf$allprop<-byvardf$devmodelfull+byvardf$specdev
  byvardf$devmodelfull_prop<-byvardf$devmodelfull/byvardf$allprop
  byvardf$specdev_prop<-byvardf$specdev/byvardf$allprop
  
  byvardf$specdev_propofage<-byvardf$specdev/byvardf$devmodelbaseage
  byvardf$nonspecdev_propofage<-1-byvardf$specdev_propofage
  byvardf$nonspecdev<-byvardf$devmodelbaseage-byvardf$specdev
  
  allreturn<-list(byvardf=byvardf)
  return(allreturn)
}

propspecbyage_gam_signflip_iter<-function(df,agevar,comparedf,iter=1000,cores=5){
  ###comparedf expects primarvar (measure xi) and othervarsnooutcome (meausures [!=xi]) for common var estimation
  ###comparedf expects primarvar (measure xi) and othervarsnooutcome_signflip (meausures [!=xi], acc/lat for opposite signflipped) for common var estimation
  message_parallel <- function(...){
    system(sprintf('echo "\n%s\n"', paste0(..., collapse="")))
  }
  byvar<-lapply(1:nrow(comparedf),function(ci){
    print(comparedf$outcome[ci])
    othervars<-unlist(lapply(strsplit(comparedf$othervarsnooutcome[ci]," "), function(x){x[!x ==""]}))
    othervarssignflip<-unlist(lapply(strsplit(comparedf$othervarsnooutcome_signflip[ci]," "), function(x){x[!x ==""]}))
    dfforcomposite<-df
    dfforcomposite[othervarssignflip]<-lapply(dfforcomposite[othervarssignflip],function(x){x*-1}) ###flip sign for sign flips
    allcompositevars<-c(othervars,othervarssignflip)
    acrossvarnumbers<-lapply(2:length(allcompositevars),function(thismanyvars){
      message_parallel(thismanyvars)
      acrossiters<-mclapply(1:iter,mc.cores=cores,function(ii){
        message_parallel(ii)
        thesevarsforcomposite<-sample(allcompositevars,thismanyvars,replace=FALSE)
        df$compositei<-compositecols(thesevarsforcomposite,dfforcomposite)
        
        if(length(which(!is.na(unique(df[,comparedf$outcome[ci]])))) < 10){
          print("reducing edf to match max dim of outcome")
          edf<-length(which(!is.na(unique(df[,comparedf$outcome[ci]])))) 
          modbaseage<-mgcv::gam(as.formula(sprintf("%s~s(%s,k=%s)",agevar,comparedf$outcome[ci],edf)),data=df) 
          devmodelbaseage<-summary(modbaseage)$dev.expl
          
          modcomp<-mgcv::gam(as.formula(sprintf("%s~s(%s)",agevar,"compositei")),data=df)
          devmodelcomp<-summary(modcomp)$dev.expl
          
          fullmod<-mgcv::gam(as.formula(sprintf("%s~s(%s,k=%s)+s(%s)",agevar,comparedf$outcome[ci],edf,"compositei")),data=df)
          devmodelfull<-summary(fullmod)$dev.expl
        }else{
          modbaseage<-mgcv::gam(as.formula(sprintf("%s~s(%s)",agevar,comparedf$outcome[ci])),data=df) 
          devmodelbaseage<-summary(modbaseage)$dev.expl
          
          modcomp<-mgcv::gam(as.formula(sprintf("%s~s(%s)",agevar,"compositei")),data=df)
          devmodelcomp<-summary(modcomp)$dev.expl
          
          fullmod<-mgcv::gam(as.formula(sprintf("%s~s(%s)+s(%s)",agevar,comparedf$outcome[ci],"compositei")),data=df)
          devmodelfull<-summary(fullmod)$dev.expl
        }
        
        
        formtosinglchar<-function(x){paste(as.character(formula(x)),collapse=",")}
        outdata<-data.frame(outcome=comparedf$outcome[ci],devmodelbaseage=devmodelbaseage,ageform=formtosinglchar(modbaseage),
                            devmodelcomp=devmodelcomp,compform=formtosinglchar(modcomp),thesevarsforcomposite=paste(thesevarsforcomposite,collapse=","),
                            devmodelfull=devmodelfull,fullform=formtosinglchar(fullmod),iter=ii)
        return(outdata)
      })
      acrossitersdf<-do.call(rbind,acrossiters)
      acrossitersdf$thismanyvars<-thismanyvars
      return(acrossitersdf)
    })    
    acrossvarnumbersdf<-do.call(rbind,acrossvarnumbers)
  })
  
  byvardf<-do.call(rbind,byvar)
  byvardf$specdev<-byvardf$devmodelfull-byvardf$devmodelcomp
  byvardf$allprop<-byvardf$devmodelfull+byvardf$specdev
  byvardf$devmodelfull_prop<-byvardf$devmodelfull/byvardf$allprop
  byvardf$specdev_prop<-byvardf$specdev/byvardf$allprop
  
  byvardf$specdev_propofage<-byvardf$specdev/byvardf$devmodelbaseage
  byvardf$nonspecdev_propofage<-1-byvardf$specdev_propofage
  byvardf$nonspecdev<-byvardf$devmodelbaseage-byvardf$specdev
  
  allreturn<-list(byvardf=byvardf)
  return(allreturn)
}


lunaize_propvarplots<-function(x){
  x+
    theme_bw()+theme(text = element_text(size = 36))+
    theme(
      axis.line = element_line(colour = "black"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      axis.title.y     = element_blank(),
      #axis.ticks.y     = element_blank(),
      axis.text.y      = element_blank(),
      panel.border = element_blank(),
      legend.position  = "none")
}

message_parallel <- function(...){
  system(sprintf('echo "\n%s\n"', paste0(..., collapse="")))
}
mgcvgam_bootfromformsdeviance<-function(mods,dataforboot,iter=1000,cores=5,qnt=c(.025,.975)){
  ###note this is only for cases when a posterior simulation (preferred gam resampling) isn't specified 
  #use with consideration####
  require(parallel)
  bootalldf<-NULL
  bootall<-mclapply(1:iter,mc.cores=cores,function(x){
    message_parallel(x)
    bootdata<-dataforboot[sample(1:nrow(dataforboot),replace=TRUE),]
    acrossmods<-lapply(1:length(mods),function(fi){
      message_parallel(fi)
      modi=mods[[fi]]
      formi=modi$formula
      formtosinglchar<-function(x){paste(as.character(formula(x)),collapse=",")}
      bootmodformreturn <- tryCatch(
        {
          bootmodformi<-mgcv::gam(formi,data=bootdata)
          bootmodformreturn<-data.frame(deviance=summary(bootmodformi)$dev.expl,form=formtosinglchar(bootmodformi),
                                        formorder=sprintf("form%s",fi))
        },
        error = function(e){
          print("pass NA")
          bootmodformreturn<-data.frame(deviance=NA,form=NA,
                                        formorder=sprintf("form%s",fi))
        }
      )
      
      return(bootmodformreturn)
    })
    acrossmodsdf<-do.call(rbind,acrossmods)
    acrossmodsdf$iter<-x
    return(acrossmodsdf)
  })
  bootalldf<-data.frame(do.call(rbind,bootall))
  rownames(bootalldf)<-seq(1:nrow(bootalldf))
  bootalldf_wide<-bootalldf[,c("deviance","formorder","iter")] %>% pivot_wider(names_from = formorder, values_from = deviance) ##convert to wide format (by formulae)
  bootalldf_wide<-bootalldf_wide[,grep("form",names(bootalldf_wide))] ###limit to formulae columns
  colpairs<-expand.grid(names(bootalldf_wide),names(bootalldf_wide)) ### all possible differences
  myfunc<-function(x){length(unique(x))}
  colpairs<-colpairs[which(apply(colpairs,1,myfunc)>1),] ##limit to only differences with other variables

  bypairdiffs<-do.call(cbind,lapply(1:nrow(colpairs),function(colpairi){
    colpairi_diffname<-paste(colpairs[colpairi,"Var1"],colpairs[colpairi,"Var2"],sep="minus")
    diff<-data.frame(bootalldf_wide[,unlist(colpairs[colpairi,"Var1"])]-bootalldf_wide[,unlist(colpairs[colpairi,"Var2"])])
    names(diff)<-colpairi_diffname
    return(diff)
  }))

  bootalldf_wide<-cbind(bootalldf_wide,bypairdiffs)
  SEdevisor<-abs(qnorm(qnt[1]))###assumes qnt are symmetrical with respect to mean e.g., defaults qnt=c(.025,.975)
  bootalldf_CI<-bootalldf_wide %>% pivot_longer(values_to='v',names_to='m',cols=where(is.numeric)) %>% dplyr::group_by(m) %>%
    dplyr::summarise(mean=mean(v,na.rm=TRUE),lower=quantile(v,na.rm=TRUE,probs=qnt[1]),upper=quantile(v,na.rm=TRUE,probs=qnt[2]),se_fromqnt=(upper-lower)/SEdevisor)

  bootalldf_CIwide<-data.frame(bootalldf_CI %>% pivot_longer(cols=where(is.numeric)) %>% pivot_wider(names_from = c(m,name), values_from = value))

  bootfullreturn<-list(bootalldf_wide=bootalldf_wide,bootalldf_CIwide=bootalldf_CIwide)

  return(bootfullreturn)
}


